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
import shutil
import subprocess
import tempfile
import multiprocessing as mp

from unitem.seq_io import read_fasta, write_fasta
from unitem.utils import make_sure_path_exists

from numpy import (sum as np_sum, zeros as np_zeros, bool as np_bool)


AA_GENE_FILE_SUFFIX = '_genes.faa'
NT_GENE_FILE_SUFFIX = '_genes.fna'


class ProdigalException(Exception):
    """Exception for errors associated with running Prodigal."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class Prodigal(object):
    """Perform gene prediction using Prodigal."""

    def __init__(self, cpus):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    @staticmethod
    def version():
        try:
            proc = subprocess.Popen(['prodigal', '-v'],
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')

            _output, error = proc.communicate()

            return error.split('\n')[1].split()[1][1:-1]
        except:
            return "(version unavailable)"

    def _run_prodigal(self, gid, genome_file):
        """Run Prodigal."""

        # check if genome has previously been processed
        trans_table_file = os.path.join(
            self.output_dir,
            f'{gid}_translation_table.tsv')
        aa_gene_file = os.path.join(
            self.output_dir, f'{gid}{AA_GENE_FILE_SUFFIX}')
        nt_gene_file = os.path.join(
            self.output_dir, f'{gid}{NT_GENE_FILE_SUFFIX}')
        gff_file = os.path.join(self.output_dir, f'{gid}.gff')

        if (os.path.exists(trans_table_file)
            and os.path.exists(aa_gene_file)
            and os.path.exists(nt_gene_file)
                and os.path.exists(gff_file)):
            # use previously calculated results
            return aa_gene_file, nt_gene_file, gff_file, trans_table_file

        seqs = read_fasta(genome_file)

        if len(seqs) == 0:
            self.logger.warning(
                f'Cannot call Prodigal on empty genome: {genome_file}')
            return None

        with tempfile.TemporaryDirectory() as tmp_dir:

            # determine number of bases
            total_bases = sum(len(seq) for seq in seqs.values())

            # Prodigal requires an uncompressed FASTA file
            prodigal_input = genome_file
            if genome_file.endswith('.gz'):
                prodigal_input = os.path.join(tmp_dir, f'{gid}.fna')
                write_fasta(seqs, prodigal_input)

            # there may be ^M character in the input file,
            # the following code is similar to dos2unix command to remove
            # those special characters.
            with open(prodigal_input, 'r') as fh:
                text = fh.read().replace('\r\n', '\n')

            processed_prodigal_input = os.path.join(
                tmp_dir,
                os.path.basename(prodigal_input))
            with open(processed_prodigal_input, 'w') as fh:
                fh.write(text)

            # call genes under different translation tables
            table_coding_density = {4: -1, 11: -1}
            for trans_table in table_coding_density:
                tt_tmp_dir = os.path.join(tmp_dir, str(trans_table))
                os.makedirs(tt_tmp_dir)
                aa_gene_file_tmp = os.path.join(
                    tt_tmp_dir,
                    f'{gid}{AA_GENE_FILE_SUFFIX}')
                nt_gene_file_tmp = os.path.join(
                    tt_tmp_dir,
                    f'{gid}{NT_GENE_FILE_SUFFIX}')
                gff_file_tmp = os.path.join(tt_tmp_dir, gid + '.gff')

                # check if there is sufficient bases to calculate prodigal parameters
                cmd = ['prodigal', '-m', '-q', '-f', 'gff']
                if total_bases < 100000:
                    cmd += ['-p', 'meta']  # use best pre-calculated parameters
                else:
                    cmd += ['-p', 'single']  # estimate parameters from data

                cmd += ['-g', str(trans_table)]
                cmd += ['-a', str(aa_gene_file_tmp)]
                cmd += ['-d', str(nt_gene_file_tmp)]
                cmd += ['-i', str(processed_prodigal_input)]

                with open(gff_file_tmp, 'wb') as stdout:
                    proc = subprocess.Popen(cmd,
                                            stdout=stdout,
                                            encoding='utf-8')

                    proc.communicate()
                    exit_code = proc.wait()
                    if exit_code != 0:
                        self.logger.error(
                            'Prodigal returned non-zero exit code.')
                        self.logger.error(f'Skipping genome {gid}.')
                        return None

                prodigalParser = ProdigalGeneFeatureParser(gff_file_tmp)

                codingBases = 0
                for seq_id, _seq in seqs.items():
                    codingBases += prodigalParser.coding_bases(seq_id)

                codingDensity = float(codingBases) / total_bases
                table_coding_density[trans_table] = codingDensity

            # determine best translation table
            trans_table_out = open(trans_table_file, 'w')
            trans_table_out.write('Table\tCoding density\n')
            trans_table_out.write(f'4\t{table_coding_density[4]:.3f}\n')
            trans_table_out.write(f'11\t{table_coding_density[11]:.3f}\n')

            best_tt = 11
            if (table_coding_density[4] - table_coding_density[11] > 0.05) and table_coding_density[4] > 0.7:
                best_tt = 4

            trans_table_out.write(f'\nSelected table\t{best_tt}\n')

            # copy results for best translation table to desired output files
            best_tt_dir = os.path.join(tmp_dir, str(best_tt))
            shutil.copyfile(os.path.join(
                best_tt_dir, f'{gid}{AA_GENE_FILE_SUFFIX}'), aa_gene_file)
            shutil.copyfile(os.path.join(
                best_tt_dir, f'{gid}{NT_GENE_FILE_SUFFIX}'), nt_gene_file)
            shutil.copyfile(os.path.join(best_tt_dir, f'{gid}.gff'), gff_file)

        return aa_gene_file, nt_gene_file, gff_file, trans_table_file

    def _worker(self, out_dict, worker_queue, writer_queue):
        """Run Prodigal in parallel across genomes."""

        while True:
            data = worker_queue.get(block=True, timeout=None)
            if data is None:
                break

            gid, genome_file = data

            out_files = self._run_prodigal(gid, genome_file)

            # Only proceed if an error didn't occur when running Prodigal
            if out_files:
                prodigal_output = {"aa_gene_path": out_files[0],
                                   "nt_gene_path": out_files[1],
                                   "gff_path": out_files[2],
                                   "translation_table_path": out_files[3]}

                out_dict[gid] = prodigal_output

            writer_queue.put(gid)

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

    def run(self, genomic_files, output_dir):
        """Run Prodigal across a set of genomes."""

        self.output_dir = output_dir
        make_sure_path_exists(output_dir)

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for gid, genome_file in genomic_files.items():
            worker_queue.put((gid, genome_file))

        for _ in range(self.cpus):
            worker_queue.put(None)

        worker_proc = []
        writer_proc = None
        try:
            manager = mp.Manager()
            out_dict = manager.dict()

            worker_proc = [mp.Process(target=self._worker, args=(out_dict,
                                                                 worker_queue,
                                                                 writer_queue))
                           for _ in range(self.cpus)]
            writer_proc = mp.Process(target=self._writer, args=(len(genomic_files),
                                                                writer_queue))

            writer_proc.start()
            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

                # Gracefully terminate the program.
                if p.exitcode != 0:
                    raise ProdigalException(
                        'Prodigal returned non-zero exit code.')

            writer_queue.put(None)
            writer_proc.join()
        except Exception as e:
            for p in worker_proc:
                p.terminate()

            if writer_proc:
                writer_proc.terminate()

            raise ProdigalException(
                f'Exception caught while running Prodigal: {e}')

        # report genomes which failed to have any genes called
        result_dict = dict()
        failed_gids = list()
        for gid, gid_dict in out_dict.items():
            if os.path.getsize(gid_dict['aa_gene_path']) <= 1:
                failed_gids.append(gid)
            else:
                result_dict[gid] = gid_dict

        if len(failed_gids) > 0:
            self.logger.warning(f'Skipping {len(failed_gids)} of {len(genomic_files)} '
                                'genomes as no genes were called by Prodigal.')

        return result_dict


class ProdigalGeneFeatureParser(object):
    """Parses prodigal gene feature files (GFF) output."""

    def __init__(self, filename):
        """Initialization.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """

        self.genes = {}
        self.last_coding_base = {}

        self._parseGFF(filename)

        self.coding_base_masks = {}
        for seq_id in self.genes:
            self.coding_base_masks[seq_id] = self._build_coding_base_mask(
                seq_id)

    def _parseGFF(self, filename):
        """Parse genes from GFF file.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        bGetTranslationTable = True
        with open(filename, 'r') as fh:
            for line in fh:
                if bGetTranslationTable and line.startswith('# Model Data'):
                    data_model_info = line.split(':')[1].strip().split(';')
                    dict_data_model = {}
                    for item in data_model_info:
                        k = item.split('=')[0]
                        v = item.split('=')[1]
                        dict_data_model[k] = v

                    self.translationTable = int(
                        dict_data_model.get('transl_table'))
                    bGetTranslationTable = False

                if line[0] == '#':
                    continue

                line_split = line.split('\t')
                seq_id = line_split[0]
                if seq_id not in self.genes:
                    geneCounter = 0
                    self.genes[seq_id] = {}
                    self.last_coding_base[seq_id] = 0

                geneId = seq_id + '_' + str(geneCounter)
                geneCounter += 1
                start = int(line_split[3])
                end = int(line_split[4])

                self.genes[seq_id][geneId] = [start, end]
                self.last_coding_base[seq_id] = max(
                    self.last_coding_base[seq_id], end)

    def _build_coding_base_mask(self, seq_id):
        """Build mask indicating which bases in a sequences are coding.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        """

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        coding_base_mask = np_zeros(
            self.last_coding_base[seq_id], dtype=np_bool)
        for pos in self.genes[seq_id].values():
            coding_base_mask[pos[0]:pos[1] + 1] = True

        return coding_base_mask

    def coding_bases(self, seq_id, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end).

        To process the entire sequence set start to 0, and
        end to None.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        start : int
            Start calculation at this position in sequence.
        end : int
            End calculation just before this position in the sequence.
        """

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end is None:
            end = self.last_coding_base[seq_id]

        return np_sum(self.coding_base_masks[seq_id][start:end])
