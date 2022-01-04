###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import logging
import subprocess
import multiprocessing as mp
import inspect

from unitem.external.tophits import TopHitTigrFile
from unitem.external.prodigal import AA_GENE_FILE_SUFFIX

TIGRFAM_SUFFIX = "_tigrfam.tsv"
TIGRFAM_OUT = '_tigrfam.out'


class TigrfamException(Exception):
    """Exception for errors associated with running HMMER with TIGRfam marker genes."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class TigrfamSearch(object):
    """Runs TIGRFAM HMMs over a set of genomes."""

    def __init__(self, marker_dir, cpus):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

        # file with TIGRfam HMMs
        if marker_dir is None:
            hmm_file = os.path.join(TigrfamSearch.hmm_dir(), 'tigrfam.hmm')
            self.marker_dir = TigrfamSearch.hmm_dir()
            self.tigrfam_hmms = os.path.abspath(hmm_file)
        else:
            self.marker_dir = marker_dir
            self.tigrfam_hmms = os.path.abspath(
                os.path.join(marker_dir, 'tigrfam.hmm'))

        self.logger.info(f' - using TIGRfam HMMs in: {self.marker_dir}')

    @staticmethod
    def version():
        """ get HMMER version."""
        try:
            proc = subprocess.Popen(['hmmsearch', '-h'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    encoding='utf-8')

            output, _ = proc.communicate()
            for line in output.split('\n'):
                if line.startswith('# HMMER'):
                    version = line.split(';')[0].replace('# HMMER', '').strip()
                    return version
            return "(version unavailable)"
        except:
            return "(version unavailable)"

    @staticmethod
    def hmm_dir():
        """Get directory with TIGRAM HMM model files."""

        src_file_path = os.path.split(inspect.getfile(lambda: None))[0]
        hmm_dir = os.path.abspath(
            os.path.join(src_file_path, '..', '..', 'markers'))

        return hmm_dir

    def get_marker_genes(self):
        """Get bacterial and archaeal PFAM marker genes."""

        bac_ms = set()
        with open(os.path.join(self.marker_dir, 'tigr_bac.lst')) as f:
            for line in f:
                bac_ms.add(line.strip())

        ar_ms = set()
        with open(os.path.join(self.marker_dir, 'tigr_ar.lst')) as f:
            for line in f:
                ar_ms.add(line.strip())

        return bac_ms, ar_ms

    def _top_hit(self, tigrfam_file, output_dir):
        """Determine top hits to TIGRFAMs.

        A gene is assigned to a single TIGRFAM
        family. This will be the top hit among
        all TIGRFAM HMMs and pass the threshold
        for the HMM.

        Parameters
        ----------
        tigrfam_file : str
            Name of file containing hits to TIGRFAM HMMs.
        """

        filename = os.path.split(tigrfam_file)[1]
        gid = filename.replace(TIGRFAM_SUFFIX, '')
        tophit_file = TopHitTigrFile(output_dir, gid)

        # populate the top hit file
        with open(tigrfam_file, 'r') as fh_tigrfam:
            for line in fh_tigrfam:
                if line[0] == '#':
                    continue

                line_split = line.split()
                gene_id = line_split[0]
                hmm_id = line_split[3]
                evalue = float(line_split[4])
                bitscore = float(line_split[5])
                tophit_file.add_hit(gene_id, hmm_id, evalue, bitscore)

        # write the top-hit file to disk
        tophit_file.write()

    def _worker(self, gene_dir, output_dir, queue_in, queue_out):
        """Process each genome in parallel."""

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file is None:
                break

            filename = os.path.split(gene_file)[1]
            output_hit_file = os.path.join(gene_dir, filename.replace(AA_GENE_FILE_SUFFIX,
                                                                      TIGRFAM_SUFFIX))

            hmmsearch_out = os.path.join(gene_dir, filename.replace(AA_GENE_FILE_SUFFIX,
                                                                    TIGRFAM_OUT))

            if os.path.exists(output_hit_file) and os.path.exists(hmmsearch_out):
                # use previously calculated results
                queue_out.put(gene_file)
                continue

            # identify TIGRfam marker genes
            cmd = ['hmmsearch', '--noali', '--notextw', '--cut_nc',
                   '-o', hmmsearch_out,
                   '--tblout', output_hit_file,
                   '--cpu', str(self.cpus_per_genome),
                   self.tigrfam_hmms, gene_file]
            p = subprocess.Popen(cmd, stderr=subprocess.PIPE, encoding='utf-8')
            _stdout, stderr = p.communicate()

            if p.returncode != 0:
                self.logger.error(
                    f'Non-zero exit code returned when running hmmsearch: {stderr}')
                self.logger.error(f'Skipping {filename}.')

                return None

            # identify top hit for each gene
            self._top_hit(output_hit_file, output_dir)

            # allow results to be processed or written to file
            queue_out.put(gene_file)

    def _writer(self, num_genomes, writer_queue):
        """Track progress."""

        num_processed = 0
        while True:
            data = writer_queue.get(block=True, timeout=None)
            if data is None:
                break

            num_processed += 1
            status = f' -> processed {num_processed:,} of {num_genomes:,} genomes'
            sys.stdout.write(f'{status}\r')
            sys.stdout.flush()

    def run(self, gene_files, gene_dir, output_dir):
        """Annotate genes with TIGRFAM HMMs.

        Parameters
        ----------
        gene_files : iterable
            Gene files in FASTA format to process.
        """

        if len(gene_files) == 0:
            self.logger.error('No genomes to process.')
            sys.exit(1)

        self.cpus_per_genome = max(1, int(self.cpus / len(gene_files)))

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for f in gene_files:
            worker_queue.put(f)

        for _ in range(self.cpus):
            worker_queue.put(None)

        try:
            worker_proc = [mp.Process(target=self._worker, args=(gene_dir,
                                                                 output_dir,
                                                                 worker_queue,
                                                                 writer_queue)) for _ in range(self.cpus)]

            write_proc = mp.Process(target=self._writer, args=(len(gene_files),
                                                               writer_queue))

            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

                if p.exitcode != 0:
                    raise TigrfamException(
                        'HMMER returned non-zero exit code.')

            writer_queue.put(None)
            write_proc.join()
        except Exception:
            for p in worker_proc:
                p.terminate()

            write_proc.terminate()
            raise
