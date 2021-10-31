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
import gzip
from collections import defaultdict

from unitem.seq_io import read_fasta


class BinTools():
    """Functions for exploring and modifying bins."""

    def __init__(self, threads=1):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

    def bin_id_from_filename(self, filename):
        """Extract bin id from bin filename."""

        bin_id = os.path.basename(filename)
        bin_id = os.path.splitext(bin_id)[0]

        return bin_id

    def _read_seq_ids(self, bin_files):
        """Read sequence IDs of all bin files."""

        bins = {}
        for bin_file in bin_files:
            bin_id = self.bin_id_from_filename(bin_file)
            bins[bin_id] = set(read_fasta(bin_file).keys())

        return bins

    def bin_files(self, bin_dir, bin_extension):
        """Get bins in directory based on extension."""

        bin_files = []
        if bin_dir is not None:
            all_files = os.listdir(bin_dir)
            for f in all_files:
                if f.endswith(bin_extension):
                    bin_file = os.path.join(bin_dir, f)
                    if os.stat(bin_file).st_size == 0:
                        self.logger.warning(
                            "Skipping bin %s as it has a size of 0 bytes." % f)
                    else:
                        bin_files.append(bin_file)

        if not bin_files:
            self.logger.error(
                "No bins found. Check the extension (-x) used to identify bins.")
            sys.exit(1)

        return sorted(bin_files)

    def unique(self, bin_files):
        """Check if sequences are assigned to multiple bins."""

        # read sequence IDs from all bins,
        # while checking for duplicate sequences within a bin
        bin_seqs = {}
        for f in bin_files:
            bin_id = self.bin_id_from_filename(f)

            if f.endswith('.gz'):
                open_file = gzip.open
            else:
                open_file = open

            seq_ids = set()
            for line in open_file(f):
                if line[0] == '>':
                    seq_id = line[1:].split(None, 1)[0]

                    if seq_id in seq_ids:
                        print('Sequence %s found multiple times in bin %s.' %
                              (seq_id, bin_id))
                    seq_ids.add(seq_id)

            bin_seqs[bin_id] = seq_ids

        # check for sequences assigned to multiple bins
        print()
        bDuplicates = False
        bin_ids = bin_seqs.keys()
        for i in range(0, len(bin_ids)):
            for j in range(i + 1, len(bin_ids)):
                seq_inter = set(bin_seqs[bin_ids[i]]).intersection(
                    set(bin_seqs[bin_ids[j]]))

                if len(seq_inter) > 0:
                    bDuplicates = True
                    print('Sequences shared between %s and %s: ' %
                          (bin_ids[i], bin_ids[j]))
                    for seq_id in seq_inter:
                        print('  ' + seq_id)
                    print()

        if not bDuplicates:
            print('No sequences assigned to multiple bins.')

    def _binning_stats(self, bins, seq_lens):
        """Calculate binning stats for a set of bins."""

        total_uniq_binned_seqs = 0
        tota_uniq_binned_bases = 0

        bin_stats = {}
        processed_seqs = set()
        repeats = set()
        for bin_id, seqs in bins.items():
            num_binned_bases = 0
            for seq_id in seqs:
                num_binned_bases += seq_lens[seq_id]
                if seq_id not in processed_seqs:
                    processed_seqs.add(seq_id)
                    tota_uniq_binned_bases += seq_lens[seq_id]
                    total_uniq_binned_seqs += 1
                else:
                    repeats.add(seq_id)

            bin_stats[bin_id] = [len(seqs), num_binned_bases]

        return bin_stats, total_uniq_binned_seqs, tota_uniq_binned_bases, len(repeats)

    def compare(self, bin_files1, bin_files2, assembly_file, output_file):
        """Compare bins from two different binning methods."""

        # determine total number of sequences
        self.logger.info('Reading bins.')
        seqs = read_fasta(assembly_file)

        seq_lens = {}
        total_bases = 0
        num_seq1K = 0
        total_bases1K = 0
        num_seq5K = 0
        total_bases5K = 0
        for seq_id, seq in seqs.items():
            seq_len = len(seq)
            seq_lens[seq_id] = seq_len
            total_bases += seq_len
            if seq_len >= 1000:
                num_seq1K += 1
                total_bases1K += seq_len
            if seq_len >= 5000:
                num_seq5K += 1
                total_bases5K += seq_len

        # determine sequences in each bin
        bins1 = self._read_seq_ids(bin_files1)
        bins2 = self._read_seq_ids(bin_files2)

        # determine bin stats
        bin_stats1, total_uniq_binned_seqs1, tota_uniq_binned_bases1, num_repeats1 = self._binning_stats(
            bins1, seq_lens)
        bin_stats2, total_uniq_binned_seqs2, tota_uniq_binned_bases2, num_repeats2 = self._binning_stats(
            bins2, seq_lens)

        # sort bins by size
        bin_stats1 = sorted(bin_stats1.items(),
                            key=lambda x: x[1][1], reverse=True)
        bin_stats2 = sorted(bin_stats2.items(),
                            key=lambda x: x[1][1], reverse=True)

        # report summary results
        print()
        print('Assembled sequences = %d (%.2f Mbp)' %
              (len(seqs), total_bases / 1e6))
        print('  No. seqs > 1 kbp = %d (%.2f Mbp)' %
              (num_seq1K, total_bases1K / 1e6))
        print('  No. seqs > 5 kbp = %d (%.2f Mbp)' %
              (num_seq5K, total_bases5K / 1e6))
        print()
        print('Binning statistics:')
        print('  1) No. bins: %s, No. binned seqs: %d (%.2f%%), No. binned bases: %.2f Mbp (%.2f%%), No. seqs in multiple bins: %d'
              % (len(bins1),
                 total_uniq_binned_seqs1,
                 total_uniq_binned_seqs1 * 100 / len(seqs),
                 tota_uniq_binned_bases1 / 1e6,
                 tota_uniq_binned_bases1 * 100 / total_bases,
                 num_repeats1))
        print('  2) No. bins: %s, No. binned seqs: %d (%.2f%%), No. binned bases: %.2f Mbp (%.2f%%), No. seqs in multiple bins: %d'
              % (len(bins2),
                 total_uniq_binned_seqs2,
                 total_uniq_binned_seqs2 * 100 / len(seqs),
                 tota_uniq_binned_bases2 / 1e6,
                 tota_uniq_binned_bases2 * 100 / total_bases,
                 num_repeats2))
        print()

        # output report
        fout = open(output_file, 'w')
        for data in bin_stats2:
            fout.write('\t' + data[0])
        fout.write(
            '\tUnbinned\tNo. Sequences\tNo. Bases (Mbp)\tBest Match\tBases in Common (%)\tSequences in Common (%)\n')

        max_bases_in_common2 = defaultdict(int)
        max_seqs_in_common2 = defaultdict(int)
        best_matching_bins2 = {}
        binned_seqs2 = defaultdict(set)
        for data1 in bin_stats1:
            bin_id1 = data1[0]
            fout.write(bin_id1)

            seqs1 = bins1[bin_id1]

            max_bases_in_common = 0
            max_seqs_in_common = 0
            best_matching_bin = 'n/a'
            binned_seqs = set()
            for data2 in bin_stats2:
                bin_id2 = data2[0]
                seqs2 = bins2[bin_id2]

                seqs_in_common = seqs1.intersection(seqs2)
                binned_seqs.update(seqs_in_common)
                num_seqs_in_common = len(seqs_in_common)
                fout.write('\t' + str(num_seqs_in_common))

                bases_in_common = 0
                for seq_id in seqs_in_common:
                    bases_in_common += seq_lens[seq_id]

                if bases_in_common > max_bases_in_common:
                    max_bases_in_common = bases_in_common
                    max_seqs_in_common = num_seqs_in_common
                    best_matching_bin = bin_id2

                if bases_in_common > max_bases_in_common2[bin_id2]:
                    max_bases_in_common2[bin_id2] = bases_in_common
                    max_seqs_in_common2[bin_id2] = num_seqs_in_common
                    best_matching_bins2[bin_id2] = bin_id1

                binned_seqs2[bin_id2].update(seqs_in_common)
            fout.write('\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{:.2f}\n'.format(
                len(seqs1) - len(binned_seqs),
                data1[1][0],
                data1[1][1] /
                1e6,
                best_matching_bin,
                max_bases_in_common *
                100 /
                data1[1][1],
                max_seqs_in_common *
                100 /
                data1[1][0],
            ))

        fout.write('Unbinned')
        for data in bin_stats2:
            binId = data[0]
            fout.write('\t%d' % (len(bins2[binId]) - len(binned_seqs2[binId])))
        fout.write('\n')

        fout.write('No. Sequences')
        for data in bin_stats2:
            fout.write('\t%d' % data[1][0])
        fout.write('\n')

        fout.write('No. Bases (Mbp)')
        for data in bin_stats2:
            fout.write('\t%.2f' % (data[1][1] / 1e6))
        fout.write('\n')

        fout.write('Best Match')
        for data in bin_stats2:
            binId = data[0]
            fout.write('\t%s' % best_matching_bins2.get(binId, 'n/a'))
        fout.write('\n')

        fout.write('Bases in Common (%)')
        for data in bin_stats2:
            binId = data[0]
            fout.write('\t%.2f' %
                       (max_bases_in_common2[binId] * 100 / data[1][1]))
        fout.write('\n')

        fout.write('Sequences in Common (%)')
        for data in bin_stats2:
            binId = data[0]
            fout.write('\t%.2f' %
                       (max_seqs_in_common2[binId] * 100 / data[1][0]))
        fout.write('\n')

        fout.close()
