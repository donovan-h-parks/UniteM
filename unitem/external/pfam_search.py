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
import inspect
import multiprocessing as mp

from unitem.external.tophits import TopHitPfamFile
from unitem.external.prodigal import AA_GENE_FILE_SUFFIX
from unitem.external.pypfam.Scan.PfamScan import PfamScan

PFAM_SUFFIX = "_pfam.tsv"
PFAM_OUT = '_pfam.out'


class PfamException(Exception):
    """Exception for errors associated with running HMMER with Pfam marker genes."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class PfamSearch(object):
    """Runs pfam_search.pl over a set of genomes."""

    def __init__(self, marker_dir, cpus):
        """Initialization."""

        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')

        if marker_dir is None:
            self.marker_dir = PfamSearch.hmm_dir()
        else:
            self.marker_dir = marker_dir

        self.logger.info(f' - using Pfam HMMs in: {self.marker_dir}')

    @staticmethod
    def hmm_dir():
        """Get directory with PFAM HMM model files."""

        src_file_path = os.path.split(inspect.getfile(lambda: None))[0]
        hmm_dir = os.path.abspath(
            os.path.join(src_file_path, '..', '..', 'markers'))

        return hmm_dir

    def get_marker_genes(self):
        """Get bacterial and archaeal PFAM marker genes."""

        bac_ms = set()
        with open(os.path.join(self.marker_dir, 'pfam_bac.lst')) as f:
            for line in f:
                bac_ms.add(line.strip())

        ar_ms = set()
        with open(os.path.join(self.marker_dir, 'pfam_ar.lst')) as f:
            for line in f:
                ar_ms.add(line.strip())

        return bac_ms, ar_ms

    def _top_hit(self, pfam_file, output_dir):
        """Determine top hits to PFAMs.

        A gene may be assigned to multiple
        PFAM families from the same clan. The
        search_pfam.pl script takes care of
        most of these issues and here the results
        are simply parsed.

        Parameters
        ----------
        pfam_file : str
            Name of file containing hits to PFAM HMMs.
        """

        filename = os.path.basename(pfam_file)
        gid = filename.replace(PFAM_SUFFIX, '')
        tophit_file = TopHitPfamFile(output_dir, gid)

        with open(pfam_file, 'r') as fh_pfam:
            for line in fh_pfam:
                if line[0] == '#' or not line.strip():
                    continue

                line_split = line.split()
                gene_id = line_split[0]
                hmm_id = line_split[5]
                evalue = float(line_split[12])
                bitscore = float(line_split[11])
                tophit_file.add_hit(gene_id, hmm_id, evalue, bitscore)

        tophit_file.write()

    def _worker(self, gene_dir, output_dir, queue_in, queue_out):
        """Process each genome in parallel."""

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file is None:
                break

            filename = os.path.basename(gene_file)
            output_hit_file = os.path.join(gene_dir, filename.replace(AA_GENE_FILE_SUFFIX,
                                                                      PFAM_SUFFIX))

            if os.path.exists(output_hit_file):
                # use previously calculated results
                queue_out.put(gene_file)
                continue

            # identify Pfam marker genes
            pfam_scan = PfamScan(cpu=self.cpus_per_genome,
                                 fasta=gene_file,
                                 dir=self.marker_dir)
            pfam_scan.search()
            pfam_scan.write_results(output_hit_file, None, None, None, None)

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
        """Annotate genes with Pfam HMMs.

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
            worker_proc = [mp.Process(target=self._worker, args=(
                gene_dir,
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
                    raise PfamException(
                        'HMMER returned non-zero exit code.')

            writer_queue.put(None)
            write_proc.join()
        except Exception:
            for p in worker_proc:
                p.terminate()

            write_proc.terminate()
            raise
