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
import copy
import random
import gzip
from shutil import copyfile
from collections import defaultdict

from unitem.defaults import *
from unitem.common import (read_bins,
                           calculateN50L50M50)
from unitem.utils import make_sure_path_exists
from unitem.markers import Markers
from unitem.plot_common_bases import PlotCommonBases
from unitem.tree_common_bases import TreeCommonBases


class Ensemble():
    """Perform ensemble binning of genomes across multiple binning methods."""

    def __init__(self, bin_prefix):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.bin_prefix = bin_prefix

    def _update_bins(self, bins, cids_to_remove):
        """Remove specified contigs from all bins.

        Parameters
        ----------
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        cids_to_remove : iterable
          Contigs to remove from bins.
        """

        cids_to_remove = set(cids_to_remove)

        # remove scaffolds in highest quality bin from all other bins
        for method_id in bins:
            for bin_id in bins[method_id]:
                bins[method_id][bin_id] -= cids_to_remove

    def _update_gene_tables(self, gene_tables, cids_to_remove):
        """Remove specified contigs from marker gene tables.

        Parameters
        ----------
        gene_tables : d[binning method] -> gene table
          Marker gene tables for bins in each binning method.
        cids_to_remove : iterable
          Contigs to remove from marker gene tables.
        """

        cids_to_remove = set(cids_to_remove)

        # remove contigs from other bins
        for bm_gene_tables in gene_tables.values():
            for bin_id in bm_gene_tables:
                updated_table = defaultdict(list)
                for marker_id, contig_ids in bm_gene_tables[bin_id].items():
                    for cid in contig_ids:
                        if cid not in cids_to_remove:
                            updated_table[marker_id].append(cid)

                bm_gene_tables[bin_id] = updated_table

    def _bin_quality(self, bins, contigs, gene_tables, quality_weight, markers):
        """Determine estimated completeness, contamination, and quality of bins.

        Parameters
        ----------
        bins : d[binning method][bin ID] -> set(cid1, cid2, ... cidN)
          Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
          Contigs across all bins.
        gene_tables : d[binning method] -> gene table
          Marker gene tables for bins in each binning method.
        quality_weight : float
          Weight given to contamination when assessing genome quality.
        markers : Markers
          Marker genes used to evaluate quality of genomes.

        Return
        ------
          List with bin metadata sorted by quality, then N50, the genome size.
        """

        q = []
        for binning_method in gene_tables:
            for bid, gene_table in gene_tables[binning_method].items():
                domain, comp, cont = markers.evaluate(gene_table)

                bin_seqs = [contigs[cid] for cid in bins[binning_method][bid]]
                n50, _l50, _m50 = calculateN50L50M50(bin_seqs)
                genome_size = sum([len(s) for s in bin_seqs])

                q.append((binning_method,
                          bid,
                          domain,
                          comp,
                          cont,
                          comp - quality_weight*cont,
                          n50,
                          genome_size))

        # sort bins by quality follwed by N50 followed by genome size, and
        # break remaining ties randomly
        q.sort(key=lambda x: (x[5], x[6], x[7], random.random()),
               reverse=True)

        return q

    def _matched_bin_sets(self, bins, contig_lens, bin_quality, min_quality, min_comp, max_cont, no_bin_matching):
        """Determine all sets of matched bins.

        Parameters
        ----------
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contig_lens : d[cid] -> length of contig
          Length of contigs.
        bin_quality : list of tuples
          Bin quality information. Must be sorted by quality!
        min_quality : float
          Minimum quality of bin to consider during bin matching.
        no_bin_matching : boolean
          Flag indicating bin matching should be skipped and all bins treated independently.

        Return
        ------
          Matched bin sets.
        """

        matched_sets = []
        processed_bins = defaultdict(set)
        for cur_bm, cur_bid, _domain, comp, cont, quality, N50, gs in bin_quality:
            if cur_bid in processed_bins[cur_bm]:
                continue  # bin has already been considered

            matched_bins = []
            matched_bins.append((cur_bm, cur_bid, quality, N50, gs))
            if not no_bin_matching and (quality >= min_quality and comp >= min_comp and cont <= max_cont):
                for test_bm, test_bid, _domain, comp, cont, quality, N50, gs in bin_quality:
                    if test_bm == cur_bm:
                        # can't group with a bin from the same binning method
                        continue

                    if test_bid in processed_bins[test_bm]:
                        continue  # bin has already been considered

                    if quality < min_quality or comp < min_comp or cont > max_cont:
                        continue

                    bp_in_common, total_bp1, total_bp2 = self._bases_in_common(bins[cur_bm][cur_bid],
                                                                               bins[test_bm][test_bid],
                                                                               contig_lens)

                    per_bp_in_common = bp_in_common / max(total_bp1, total_bp2)
                    if per_bp_in_common > 0.5:
                        matched_bins.append(
                            (test_bm, test_bid, quality, N50, gs))

            # removed matched bins
            for bm, bid, _q, _n50, _gs in matched_bins:
                processed_bins[bm].add(bid)

            matched_sets.append(tuple(matched_bins))

        # sort by total quality of bins in matched sets
        matched_sets.sort(key=lambda ms: (sum([x[2] for x in ms]),
                                          sum([x[3] for x in ms]),
                                          sum([x[4] for x in ms]),
                                          random.random()),
                          reverse=True)

        return matched_sets

    def _bases_in_common(self, contig_ids1, contig_ids2, contig_lens):
        """Calculate number of base pairs in common between two bins.

        Parameters
        ----------
        contig_ids1 : iterable
          Contigs in bin.
        contig_ids2 : iterable
          Contigs in bin.
        contig_lens : d[seq ID] -> length
          Length of contigs.

        Returns
        -------
          Base pairs in common, total bases in bin 1, total bases in bin 2
        """

        bp_in_common = sum([contig_lens[seq_id]
                            for seq_id in contig_ids2
                            if seq_id in contig_ids1])

        total_bp1 = sum([contig_lens[seq_id] for seq_id in contig_ids1])
        total_bp2 = sum([contig_lens[seq_id] for seq_id in contig_ids2])

        return bp_in_common, total_bp1, total_bp2

    def _resolve_matched_set(self,
                             matched_set,
                             bins,
                             contigs,
                             bin_quality,
                             remove_perc,
                             add_perc,
                             add_matches,
                             sel_min_quality,
                             sel_min_comp,
                             sel_max_cont,
                             greedy,
                             unanimous):
        """Select contigs to cluster from matched bin set.

        Parameters
        ----------
        matched_set : iterable with bin metadata (binning method, bin ID, quality)
          Matched set of bins. Highest quality bin must be first.
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
          Contigs across all bins.
        bin_quality : list of tuples
          Bin quality information. Must be sorted by quality!
        remove_perc : float
          Minimum percentage of bins from other binning methods require to remove contig in highest quality bin.
        add_perc : float
          Minimum percentage of matched bins required to add contig to highest quality bin.
        add_matches : float
          Minimum number of matched bins required to 'add' contigs.
        sel_min_quality : float
          Minimum quality of bins to consider for filtering contigs.
        greedy : boolean
          Perform greedy filtering.
        unanimous : boolean
          Perform unanimous filtering.
        """

        # identify contig count for bins in matched in set
        primary_bm, primary_bid, _q, _n50, _gs = matched_set[0]
        primary_contigs = bins[primary_bm][primary_bid]

        matched_contigs = defaultdict(int)
        matched_bin_ids = set()
        matched_methods = set()
        for bm, bid, _q, _n50, _gs in matched_set:
            matched_bin_ids.add(bm + bid)
            matched_methods.add(bm)
            for cid in bins[bm][bid]:
                matched_contigs[cid] += 1

        # add or remove contigs based on consensus filtering
        new_bin = {}
        removed_contigs = set()
        added_contigs = set()
        if unanimous:
            # remove all contigs from highest quality bin that are
            # not present in all other matched bins
            num_matched_bins = len(matched_set)
            for cid in primary_contigs:
                match_count = matched_contigs.get(cid, 0)
                if match_count == num_matched_bins:
                    new_bin[cid] = contigs[cid]
        elif greedy:
            for cid in primary_contigs:
                new_bin[cid] = contigs[cid]
        else:
            # identify contig count for bins not in matched set
            unmatched_contigs = defaultdict(int)
            binned_non_matched_method = defaultdict(int)
            for bm, bid, _domain, comp, cont, quality, _N50, _gs in bin_quality:
                if bm+bid in matched_bin_ids:
                    continue

                if quality < sel_min_quality or comp < sel_min_comp or cont > sel_max_cont:
                    continue

                for cid in bins[bm][bid]:
                    unmatched_contigs[cid] += 1

                    if bm not in matched_methods:
                        binned_non_matched_method[cid] += 1

            # remove contigs from highest quality bin iff it appears in a
            # consensus of non-matched bins
            # (highest quality bin with contig is also counted)
            num_binning_methods = len(bins)
            for cid in primary_contigs:
                unmatched_count = unmatched_contigs.get(cid, 0)

                per_other_bm = unmatched_count * 100 / num_binning_methods
                if per_other_bm <= remove_perc:
                    new_bin[cid] = contigs[cid]
                else:
                    removed_contigs.add(cid)

            # add new contigs to highest quality bin iff there is a suitable
            # number of matched bins and the contig is present in a consensus
            # of these matched bins
            num_matched_bins = len(matched_set)
            if num_matched_bins >= add_matches:
                for cid, match_count in matched_contigs.items():
                    if cid in primary_contigs:
                        continue

                    # number of bins in denominator is number of matched
                    # binning methods plus number of unmatched binning methods
                    # where contig is in a bin
                    bin_count = num_matched_bins + \
                        binned_non_matched_method.get(cid, 0)
                    p = match_count * 100 / bin_count
                    if p >= add_perc:
                        new_bin[cid] = contigs[cid]
                        added_contigs.add(cid)

        return new_bin, removed_contigs, added_contigs

    def _report_selection(self,
                          bin_name,
                          bin,
                          markers,
                          quality_weight,
                          binning_method,
                          original_bid,
                          num_matched_bins,
                          num_removed,
                          num_added,
                          greedy,
                          unanimous):
        """Create string with information about selected bin."""

        domain, comp, cont = markers.bin_quality(bin)
        quality = comp - quality_weight*cont

        genome_size = sum([len(seq) for seq in bin.values()])
        num_scaffolds = len(bin)
        n50, l50, m50 = calculateN50L50M50(bin.values())

        row = '%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f' % (bin_name,
                                                    binning_method,
                                                    original_bid,
                                                    domain,
                                                    comp,
                                                    cont,
                                                    quality)

        row += '\t%d\t%d\t%d\t%d' % (genome_size,
                                     n50,
                                     l50,
                                     num_scaffolds)

        if not greedy:
            row += '\t%s\t%s' % (num_matched_bins,
                                 num_removed)
            if not unanimous:
                row += '\t%s' % num_added

        row += '\n'

        return row

    def _bin_summary(self,
                     sel_bins,
                     merged_bins,
                     quality_weight,
                     markers,
                     summary_file):
        """Summarize quality of reported bins."""

        bins = sel_bins.copy()
        bins.update(merged_bins)

        comp_cont_count = defaultdict(lambda: defaultdict(int))
        quality_count = defaultdict(int)
        total_comp = 0
        total_cont = 0
        total_q = 0
        for bid, bin in bins.items():
            domain, comp, cont = markers.bin_quality(bin)
            quality = comp - quality_weight*cont

            for test_comp in [90, 80, 70]:
                for test_cont in [5, 10]:
                    if comp >= test_comp and cont <= test_cont:
                        comp_cont_count[test_comp][test_cont] += 1

            for test_q in [90, 70, 50]:
                if quality >= test_q:
                    quality_count[test_q] += 1

            total_comp += comp
            total_cont += cont
            total_q += quality

        fout = open(summary_file, 'w')
        fout.write('Criteria\tNo. Bins\n')
        for comp in [90, 80, 70]:
            for cont in [5, 10]:
                fout.write('Comp. >= %d, cont. <= %d\t%d\n' %
                           (comp, cont, comp_cont_count[comp][cont]))

        for q in [90, 70, 50]:
            fout.write('Quality >= %d\t%d\n' % (q, quality_count[q]))

        fout.write('Total completeness\t%.1f\n' % total_comp)
        fout.write('Total contamination\t%.1f\n' % total_cont)
        fout.write('Total quality\t%.1f\n' % total_q)

        fout.close()

    def _perc_bases_in_common(self, query_bin, bins, contig_lens):
        """Determine percentage of base pairs in common between a bin and a set of bins.

        Parameters
        ----------
        query_bin : d[cid] -> contig
          Bin to match.
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contig_lens : d[cid] -> contig length
          Length of all contigs.

        Returns
        -------
        list of tuples
          Each list indicates the binning method, bin id, and percentage of common bases.

        Only cases where the percentage of common bases is >= 50% are returned.
        """

        matches = []
        for bm in bins:
            for bid in bins[bm]:
                bp_in_common, total_bp1, total_bp2 = self._bases_in_common(query_bin.keys(),
                                                                           bins[bm][bid],
                                                                           contig_lens)

                perc_common = bp_in_common * 100 / max(total_bp1, total_bp2)
                if perc_common >= 50:
                    matches.append((bm, bid, perc_common))

        matches.sort(key=lambda x: (x[2], random.random()), reverse=True)

        return matches

    def _read_init_bin_quality(self, profile_dir):
        """Read quality of each bin as determined by the profile command."""

        bq = defaultdict(lambda: {})
        for f in os.listdir(profile_dir):
            if not f.endswith('_quality.tsv'):
                continue

            method_id = f.replace('_quality.tsv', '')
            with open(os.path.join(profile_dir, f)) as f:
                f.readline()

                for line in f:
                    tokens = line.strip().split('\t')
                    bid, domain, comp, cont, _ = tokens
                    bq[method_id][bid] = (domain, float(comp), float(cont))

        return bq

    def _write_initial_contig_state(self,
                                    methods_sorted,
                                    contigs_in_bins,
                                    init_bin_quality,
                                    output_dir):
        """Record initial state of each contig."""

        fout = open(os.path.join(output_dir, 'contig_info_initial.tsv'), 'w')
        fout.write('Contig ID')
        for bm in methods_sorted:
            fout.write('\t' + bm)
        fout.write('\n')

        for cid in contigs_in_bins:
            fout.write(cid)
            for bm in methods_sorted:
                if bm not in contigs_in_bins[cid]:
                    fout.write('\tunbinned')
                else:
                    bid = contigs_in_bins[cid][bm]
                    domain, comp, cont = init_bin_quality[bm][bid]
                    fout.write('\t%s,%.1f,%.1f' % (bid, comp, cont))
            fout.write('\n')
        fout.close()

    def run(self,
            profile_dir,
            bin_dirs,
            marker_dir,
            quality_weight,
            sel_min_quality,
            sel_min_comp,
            sel_max_cont,
            remove_perc,
            add_perc,
            add_matches,
            greedy,
            unanimous,
            report_min_quality,
            simple_headers,
            output_dir):
        """Perform ensemble binning of genomes across multiple binning methods.

        Parameters
        ----------
        profile_dir : str
          Directory with bin profiles (output of 'profile' command).
        bin_dirs : list of str
            Directories containing bins from different binning methods.
        quality_weight : float
          Weight given to contamination when assessing quality of bins.
        sel_min_quality : float
          Minimum quality of bins to select.
        sel_min_comp : float
          Minimum completeness of bins to select.
        sel_max_cont : float
          Maximum contamination of bins to select.
        remove_perc : float
          Minimum percentage of bins from other binning methods require to remove contig in highest quality bin.
        add_perc : float
          Minimum percentage of matched bins required to add contig to highest quality bin.
        add_matches : float
          Minimum number of matched bins required to 'add' contigs.
        greedy : boolean
          Perform greedy bin selection.
        unanimous : boolean
          Perform unanimous bin selection.
        report_min_quality : float
          Minimum quality of bins to report.
        simple_headers : boolean
          Flag indicating that information should not be appended to the headers of bin FASTA files.
        output_dir : str
          Output directory.
        """

        markers = Markers(marker_dir)

        if greedy:
            self.logger.info("Greedy selection with quality_weight = {:.1f}, sel_min_quality = {:.1f}, sel_min_comp = {:.1f}, and sel_max_cont = {:.1f}.".format(
                quality_weight, sel_min_quality, sel_min_comp, sel_max_cont))
            remove_perc = 101
            add_perc = 101
        elif unanimous:
            self.logger.info("Unanimous selection with quality_weight = {:.1f}, sel_min_quality = {:.1f}, sel_min_comp = {:.1f}, and sel_max_cont = {:.1f}.".format(
                quality_weight, sel_min_quality, sel_min_comp, sel_max_cont))
        else:
            self.logger.info("Consensus selection with quality_weight ={:.1f}, sel_min_quality = {:.1f}, sel_min_comp = {:.1f}, and sel_max_cont = {:.1f}.".format(
                quality_weight, sel_min_quality, sel_min_comp, sel_max_cont))
            self.logger.info('Removing and adding contigs by consensus with remove_perc = {:.1f}, add_perc = {:.1f}, add_matches = {:,}.'.format(
                remove_perc, add_perc, add_matches))

        self.logger.info(
            'Reporting bins with a quality >= {:.1f}.'.format(report_min_quality))

        # get scaffold IDs in bins across all binning methods
        self.logger.info('Reading all bins.')
        bins, contigs, contigs_in_bins = read_bins(bin_dirs)
        methods_sorted = sorted(bins.keys())
        contig_lens = {cid: len(contigs[cid]) for cid in contigs}
        orig_bins = copy.deepcopy(bins)

        # get marker genes for bins across all binning methods
        self.logger.info('Identifying marker genes across all bins.')
        gene_tables = markers.marker_gene_tables(
            profile_dir,
            binning_methods=bins.keys()
        )

        # create output directories
        bin_dir = os.path.join(output_dir, 'bins')
        make_sure_path_exists(bin_dir)

        # get initial quality of bins
        self.logger.info('Reading initial quality of bins.')
        init_bin_quality = self._read_init_bin_quality(profile_dir)

        # write out initial state of each contig
        self.logger.info('Recording initial state of contigs.')
        self._write_initial_contig_state(methods_sorted,
                                         contigs_in_bins,
                                         init_bin_quality,
                                         output_dir)

        # perform consensus selection of bins
        header = 'UniteM Bin ID\tBinning Method\tBin ID'
        header += '\tMarker Domain\tCompleteness (%)\tContamination (%)\tQuality (%)'
        header += '\tGenome Size\tN50\tL50\tNo. Contigs'
        if not greedy:
            header += '\tNo. Matched Bins\tNo. Removed Contigs'
            if not unanimous:
                header += '\tNo. Added Contigs'
        header += '\n'

        fout = open(os.path.join(output_dir, 'bin_info.tsv'), 'w')
        fout.write(header)

        fout_matched = open(os.path.join(
            output_dir, 'matched_set_info.tsv'), 'w')
        if greedy:
            fout_matched.write(
                'No. Matched Bins\tBin ID\tCompleteness (%)\tContamination (%)\tQuality (%)\n')
        else:
            fout_matched.write(
                'No. Matched Bins\tUniteM Bin ID\tCompleteness (%)\tContamination (%)\tQuality (%)')
            fout_matched.write(
                '\tMatched Bin ID\tCompleteness (%)\tContamination (%)\tQuality (%)\n')

        fout_bin_info = open(os.path.join(
            output_dir, 'matched_bin_info.tsv'), 'w')
        fout_bin_info.write(
            'UniteM Bin ID\tCompleteness (%)\tContamination (%)\tQuality (%)')
        fout_bin_info.write(
            '\tNo. Matched Bins\tMatched Bin ID\tPercent Common Bases\tCompleteness (%)\tContamination (%)\tQuality (%)\n')

        fout_contigs = open(os.path.join(output_dir, 'contig_info.tsv'), 'w')
        fout_contigs.write(
            'UniteM Bin ID\tContig ID\tNo. Matched Bins\tNo. Unmatched Bins\tNo. Degenerate Bins\tNo. Unbinned')
        for method in methods_sorted:
            fout_contigs.write(f'\t{method}')
        fout_contigs.write('\n')

        plot_dir = os.path.join(output_dir, 'common_bases')
        make_sure_path_exists(plot_dir)

        # perform ensemble binning
        self.logger.info('Performing ensemble binning.')
        bin_num = 0
        total_comp = 0
        total_cont = 0
        total_quality = 0
        sel_bins = {}
        selected_rows = {}
        unitem_common_bases = defaultdict(lambda: defaultdict(int))
        unitem_bin_quality = {}
        sanity_check_contigs = set()
        tree_common_bases = TreeCommonBases()
        while True:
            # determine highest quality match bin set
            bin_quality = self._bin_quality(bins,
                                            contigs,
                                            gene_tables,
                                            quality_weight,
                                            markers)

            matched_sets = self._matched_bin_sets(bins,
                                                  contig_lens,
                                                  bin_quality,
                                                  sel_min_quality,
                                                  sel_min_comp,
                                                  sel_max_cont,
                                                  greedy)

            if len(matched_sets) == 0:
                break  # no bins meeting selection criteria

            new_bin, removed_contigs, added_contigs = self._resolve_matched_set(
                matched_sets[0],
                bins,
                contigs,
                bin_quality,
                remove_perc,
                add_perc,
                add_matches,
                sel_min_quality,
                sel_min_comp,
                sel_max_cont,
                greedy,
                unanimous)

            _domain, comp, cont = markers.bin_quality(new_bin)

            quality = comp - quality_weight*cont
            if quality < report_min_quality:
                break

            total_comp += comp
            total_cont += cont
            total_quality += quality

            # report selection
            bin_num += 1
            unitem_bin_id = f'{self.bin_prefix}_{bin_num}'
            unitem_bin_quality[unitem_bin_id] = (comp, cont)
            primary_bm, primary_bid, _q, _n50, _gs = matched_sets[0][0]
            self.logger.info("Selected {} from {} with quality = {:.1f} (comp. = {:.1f}%, cont. = {:.1f}%).".format(
                primary_bid,
                primary_bm,
                quality,
                comp,
                cont))

            if not greedy and not unanimous:
                # performing consensus binning
                self.logger.info("-> Identified {} matched bins, removed {} contigs, added {} contigs.".format(
                    len(matched_sets[0]),
                    len(removed_contigs),
                    len(added_contigs)))
            elif not greedy:
                # performing unanimous binning
                self.logger.info("-> Identified {} matched bins and removed {} contigs.".format(
                    len(matched_sets[0]),
                    len(removed_contigs)))

            # write out matched set
            fout_matched.write(f'{len(matched_sets[0])}')
            if not greedy:
                fout_matched.write('\t{}\t{:.1f}\t{:.1f}\t{:.1f}'.format(
                    unitem_bin_id,
                    comp,
                    cont,
                    quality))

            for bm, bid, q, _n50, _gs in matched_sets[0]:
                _domain, comp, cont = markers.bin_quality(bins[bm][bid])
                quality = comp - quality_weight*cont
                fout_matched.write('\t{}\t{:.1f}\t{:.1f}\t{:.1f}'.format(
                    bm + '~' + bid,
                    comp,
                    cont,
                    quality))
            fout_matched.write('\n')

            # write out common base pairs
            matches = self._perc_bases_in_common(
                new_bin, orig_bins, contig_lens)
            fout_bin_info.write('{}\t{:.1f}\t{:.1f}\t{:.1f}\t{}'.format(
                unitem_bin_id,
                comp,
                cont,
                quality,
                len(matches)))

            for bm, bid, perc_common in matches:
                unitem_common_bases[unitem_bin_id][bm] = perc_common
                _domain, comp, cont = markers.bin_quality(orig_bins[bm][bid])
                quality = comp - quality_weight*cont
                fout_bin_info.write('\t{}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}'.format(
                    bm + '~' + bid,
                    perc_common,
                    comp,
                    cont,
                    quality))

            fout_bin_info.write('\n')

            output_file = os.path.join(plot_dir, unitem_bin_id + '.svg')
            tree_common_bases.plot(unitem_bin_id, matches, output_file)

            # write out summary information about selected bins
            row = self._report_selection(unitem_bin_id,
                                         new_bin,
                                         markers,
                                         quality_weight,
                                         primary_bm,
                                         primary_bid,
                                         str(len(matched_sets[0])),
                                         str(len(removed_contigs)),
                                         str(len(added_contigs)),
                                         greedy,
                                         unanimous)
            fout.write(row)
            selected_rows[unitem_bin_id] = row

            # write out contig info
            matched_bins = {}
            for m in matched_sets[0]:
                matched_bins[m[0]] = m[1]

            new_bin_file = os.path.join(bin_dir, unitem_bin_id + '.fna.gz')
            fout_bin = gzip.open(new_bin_file, 'wt')
            for cid, seq in new_bin.items():
                fout_contigs.write(f'{unitem_bin_id}\t{cid}')

                if cid in sanity_check_contigs:
                    self.logger.error(f'Contig selected twice: {cid}')
                    sys.exit(1)
                sanity_check_contigs.add(cid)

                matched = 0
                unmatched = 0
                unbinned = 0
                degenerate = 0
                row = ''
                for m in methods_sorted:
                    if m in contigs_in_bins[cid]:
                        bid = contigs_in_bins[cid][m]
                        domain, comp, cont = markers.bin_quality(bins[m][bid])
                        if m in matched_bins and bid == matched_bins[m]:
                            row += '\t{},{},{:.1f},{:.1f}'.format(
                                'matched', bid, comp, cont)
                            matched += 1
                        else:
                            q = comp - quality_weight*cont
                            if q < sel_min_quality or comp < sel_min_comp or cont > sel_max_cont:
                                row += '\t{},{},{:.1f},{:.1f}'.format(
                                    'degenerate', bid, comp, cont)
                                degenerate += 1
                            else:
                                row += '\t{},{},{:.1f},{:.1f}'.format(
                                    'unmatched', bid, comp, cont)
                                unmatched += 1
                    else:
                        row += '\tunbinned'
                        unbinned += 1

                fout_contigs.write('\t{}\t{}\t{}\t{}'.format(
                    matched,
                    unmatched,
                    degenerate,
                    unbinned))
                fout_contigs.write(row + '\n')

                if simple_headers:
                    fout_bin.write(f'>{cid}\n')
                else:
                    cid_info = '[matched={}] [unmatched={}] [degenerate={}] [unbinned={}]'.format(
                        matched,
                        unmatched,
                        degenerate,
                        unbinned)

                    if cid in added_contigs:
                        cid_info += ' [added by consensus]'

                    fout_bin.write(f'>{cid} {cid_info}\n')
                fout_bin.write(seq + '\n')
            fout_bin.close()

            # remove contigs in highest quality bin from marker gene tables and all other bins
            self._update_gene_tables(gene_tables, new_bin.keys())
            self._update_bins(bins, new_bin.keys())

            sel_bins[unitem_bin_id] = new_bin

        self.logger.info(f'Selected {bin_num} bins.')
        self.logger.info('-> total comp. = {:.1f}, total cont. = {:.1f}, total quality = {:.1f}'.format(
            total_comp,
            total_cont,
            total_quality))

        fout.close()
        fout_matched.close()
        fout_contigs.close()
        fout_bin_info.close()

        plot = PlotCommonBases()
        output_plot = os.path.join(output_dir, 'percent_common_bases.svg')
        plot.plot(unitem_common_bases, unitem_bin_quality, output_plot)

        # summarize results of reported bins
        summary_file = os.path.join(output_dir, 'bin_quality_summary.tsv')
        self._bin_summary(sel_bins,
                          {},
                          quality_weight,
                          markers,
                          summary_file)
