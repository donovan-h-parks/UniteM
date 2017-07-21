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
import random
from shutil import copyfile
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.common import make_sure_path_exists

from unitem.defaults import *
from unitem.common import (read_bins, 
                            parse_bin_stats,                         
                            calculateN50L50M50)
from unitem.markers import Markers
from unitem.merge import Merge

class Consensus():
    """Perform consensus selection of genomes across multiple binning methods."""

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
        gene_tables : d[binning method] -> (bac gene table, ar gene table)
          Bacterial and archaeal marker gene tables for bins in each binning method.
        cids_to_remove : iterable
          Contigs to remove from marker gene tables.
        """

        cids_to_remove = set(cids_to_remove)

        # remove contigs from other bins
        for binning_method, (bac_gene_tables, ar_gene_tables) in gene_tables.iteritems():
            for bin_id in bac_gene_tables:
                for marker_id, contig_ids in bac_gene_tables[bin_id].iteritems():
                    bac_gene_tables[bin_id][marker_id] = []
                    for cid in contig_ids:
                        if cid not in cids_to_remove:
                            bac_gene_tables[bin_id][marker_id].append(cid)
                            
                for marker_id, contig_ids in ar_gene_tables[bin_id].iteritems():
                    ar_gene_tables[bin_id][marker_id] = []
                    for cid in contig_ids:
                        if cid not in cids_to_remove:
                            ar_gene_tables[bin_id][marker_id].append(cid)
                            
    def _bin_quality(self, bins, contigs, gene_tables, quality_weight):
        """Determine estimated completeness, contamination, and quality of bins.
        
        Parameters
        ----------
        bins : d[binning method][bin ID] -> set(cid1, cid2, ... cidN)
          Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
          Contigs across all bins.
        gene_tables : d[binning method] -> (bac gene table, ar gene table)
          Bacterial and archaeal marker gene tables for bins in each binning method.
        quality_weight : float
          Weight given to contamination when assessing genome quality.
        
        Return
        ------
          List with bin metadata sorted by quality, then N50, the genome size.
        """
        
        markers = Markers()
        
        q = []
        for bm, (bac_gene_tables, ar_gene_tables) in gene_tables.iteritems():
            for bid in bac_gene_tables:
                domain, comp, cont = markers.evaluate(bac_gene_tables[bid], 
                                                            ar_gene_tables[bid])
                
                bin_seqs = [contigs[cid] for cid in bins[bm][bid]]
                n50, l50, m50 = calculateN50L50M50(bin_seqs)
                genome_size = sum([len(s) for s in bin_seqs])
                                                        
                q.append((bm, bid, domain, comp, cont, comp - quality_weight*cont, n50, genome_size))
                    
        # sort bins by quality follwed by N50 followed by genome size, and 
        # break remaining ties randomly
        q.sort(key=lambda x: (x[5], x[6], x[7], random.random()), reverse=True)
        
        return q
        
    def _seq_lens(self, bins, contigs):
        """Calculate length of sequences in bins.
        
        Parameters
        ----------
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
          Contigs across all bins.
          
        Return
        ------
          Dictionary indicating length of each sequence in bin.
        """
        
        seq_lens = defaultdict(dict)
        for bm in bins:
            for bid in bins[bm]:
                d = {}
                for cid in bins[bm][bid]:
                    d[cid] = len(contigs[cid])
                    
                seq_lens[bm][bid] = d
            
        return seq_lens

    def _matched_bin_sets(self, bins, contigs, bin_quality, min_quality, no_bin_matching):
        """Determine all sets of matched bins.
        
        Parameters
        ----------
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
          Contigs across all bins.
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
        
        seq_lens = self._seq_lens(bins, contigs)
        
        matched_sets = []
        for cur_bm, cur_bid, _domain, _comp, _cont, quality, N50, gs in bin_quality:
            if cur_bid not in seq_lens[cur_bm]:
                continue # bin has already been considered
                                                         
            matched_bins = []
            matched_bins.append((cur_bm, cur_bid, quality, N50, gs))
            if not no_bin_matching and quality >= min_quality:
                for test_bm, test_bid, _domain, _comp, _cont, quality, N50, gs in bin_quality:
                    if test_bm == cur_bm:
                        # can't group with a bin from the same binning method
                        continue 
                        
                    if test_bid not in seq_lens[test_bm]:
                        continue # bin has already been considered
                                                    
                    if quality < min_quality:
                        continue

                    bp_in_common, total_bp1, total_bp2 = self._bases_in_common(seq_lens[cur_bm][cur_bid], 
                                                                                seq_lens[test_bm][test_bid])
        
                    per_bp_in_common = float(bp_in_common) / total_bp1
                    if per_bp_in_common > 0.5:
                        matched_bins.append((test_bm, test_bid, quality, N50, gs))
                        
            # removed matched bins
            for bm, bid, _q, _n50, _gs in matched_bins:
                del seq_lens[bm][bid]
            
            matched_sets.append(tuple(matched_bins))
                
        # sort by total quality of bins in matched sets
        matched_sets.sort(key=lambda ms: (sum([x[2] for x in ms]),
                                            sum([x[3] for x in ms]),
                                            sum([x[4] for x in ms]),
                                            random.random()), 
                                            reverse=True)

        return matched_sets

    def _bases_in_common(self, seq_lens1, seq_lens2):
        """Calculate number of base pairs in common between two bins.
        
        Parameters
        ----------
        seq_lens1 : d[seq ID] -> seq length
          Length of sequences for first bin.
        seq_lens2 : d[seq ID] -> seq length
          Length of sequences for second bin.
        
        Returns
        -------
          Base pairs in common, total bases in bin 1, total bases in bin 2
        """
        
        seq_ids1 = set(seq_lens1.keys())
        bp_in_common = sum([seq_len
                            for seq_id, seq_len in seq_lens2.iteritems()
                            if seq_id in seq_ids1])
        
        total_bp1 = sum([seq_len for seq_len in seq_lens1.values()])
        total_bp2 = sum([seq_len for seq_len in seq_lens2.values()])
        
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
        num_removed = 0
        num_added = 0
        if unanimous:
            # remove all contigs from highest quality bin that are
            # not present in all other matched bins
            num_matched_bins = len(matched_set)
            for cid in primary_contigs:
                match_count = matched_contigs.get(cid, 0)
                if match_count == num_matched_bins:
                    new_bin[cid] = contigs[cid]
                else:
                    num_removed += 1
        elif greedy:
            for cid in primary_contigs:
                new_bin[cid] = contigs[cid]
        else:
            # identify contig count for bins not in matched set
            unmatched_contigs = defaultdict(int)
            binned_non_matched_method = defaultdict(int)
            for bm, bid, _domain, _comp, _cont, quality, _N50, _gs in bin_quality:
                if bm+bid in matched_bin_ids:
                    continue
                    
                if quality < sel_min_quality:
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

                per_other_bm = unmatched_count * 100.0 / num_binning_methods
                if per_other_bm <= remove_perc:
                    new_bin[cid] = contigs[cid]
                else:
                    num_removed += 1

            # add new contigs to highest quality bin iff there is a suitable
            # number of matched bins and the contig is present in a consensus
            # of these matched bins
            num_matched_bins = len(matched_set)
            if num_matched_bins >= add_matches:
                for cid, match_count in matched_contigs.iteritems():
                    if cid in primary_contigs:
                        continue
                        
                    # number of bins in denominator is number of matched
                    # binning methods plus number of unmatched binning methods
                    # where contig is in a bin
                    bin_count = num_matched_bins + binned_non_matched_method.get(cid, 0)
                    p = match_count * 100.0 / bin_count
                    if p >= add_perc:
                        new_bin[cid] = contigs[cid]
                        num_added += 1
                
        return new_bin, num_removed, num_added
    
    def _report_selection(self, 
                            bin_name, 
                            bin, 
                            markers,
                            report_weight,
                            binning_method,
                            original_bid,
                            num_matched_bins,
                            num_removed,
                            num_added,
                            greedy,
                            unanimous):
        """Create string with information about selected bin."""
        
        domain, comp, cont = markers.bin_quality(bin)
        quality = comp - report_weight*cont
        
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
                        report_weight, 
                        report_min_quality,
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
        for bid, bin in bins.iteritems():
            domain, comp, cont = markers.bin_quality(bin)
            quality = comp - report_weight*cont
            
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
                fout.write('Comp. >= %d, cont. <= %d\t%d\n' % (comp, cont, comp_cont_count[comp][cont]))

        for q in [90, 70, 50]:
            fout.write('Quality >= %d\t%d\n' % (q, quality_count[q]))
            
        fout.write('Total compl.\t%.1f\n' % total_comp)
        fout.write('Total cont.\t%.1f\n' % total_cont)
        fout.write('Total quality\t%.1f\n' % total_q)
        
        fout.close()
        
    def run(self, 
            profile_dir, 
            bin_dirs,
            sel_weight,
            sel_min_quality,
            remove_perc,
            add_perc,
            add_matches,
            greedy,
            unanimous,
            merge,
            report_weight, 
            report_min_quality, 
            output_dir):
                        
        """Perform consensus selection of genomes across multiple binning methods.

        Parameters
        ----------
        profile_dir : str
          Directory with bin profiles (output of 'profile' command).
        bin_dirs : list of str
            Directories containing bins from different binning methods.
        self_weight : float
          Weight given to contamination when assessing quality of bins to select.
        sel_min_quality : float
          Minimum quality of bins to select.
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
        merge : boolean
          Attempt to merge partial bins.
        report_weight : float
          Weight given to contamination when assessing quality of bins to report.
        report_min_quality : float
          Minimum quality of bins to report.
        output_dir : str
          Output directory.
        """
        
        markers = Markers()
        
        if greedy:
            self.logger.info("Performing greedy bin selection with sel_weight = %.1f and sel_min_quality = %.1f." % (sel_weight, sel_min_quality))
            remove_perc = 101
            add_perc = 101
        elif unanimous:
            self.logger.info("Performing unanimous bin selection with sel_weight = %.1f and sel_min_quality = %.1f." % (sel_weight, sel_min_quality))
        else:
            self.logger.info("Performing consensus bin selection with sel_weight = %.1f and sel_min_quality = %.1f." % (sel_weight, sel_min_quality))
            self.logger.info('Removing and adding contigs by consensus with remove_perc = %.1f, add_perc = %.1f, add_matches = %d.' % (remove_perc, add_perc, add_matches))
        
        self.logger.info('Reporting bins with weight = %.1f and min_quality = %.1f.' % (report_weight, report_min_quality))
        
        # get scaffold IDs in bins across all binning methods
        self.logger.info('Reading all bins.')
        bins, contigs = read_bins(bin_dirs)

        # get marker genes for bins across all binning methods
        self.logger.info('Identifying marker genes across all bins.')
        gene_tables = markers.marker_gene_tables(profile_dir, binning_methods=bins.keys())

        # create output directories
        sel_bin_dir = os.path.join(output_dir, 'bins_selected')
        make_sure_path_exists(sel_bin_dir)
        
        final_bin_dir = os.path.join(output_dir, 'bins_final')
        make_sure_path_exists(final_bin_dir)
        
        # perform consensus selection of bins
        header = 'UniteM Bin ID\tBinning Method\tBin ID'
        header += '\tMarker Domain\tCompleteness (%)\tContamination (%)\tQuality (%)'
        header += '\tGenome Size\tN50\tL50\tNo. Contigs'
        if not greedy:
            header += '\tNo. Matched Bins\tNo. Removed Contigs'
        if not unanimous:
            header += '\tNo. Added Contigs'
        header += '\n'
        
        fout = open(os.path.join(output_dir, 'bins_selected.tsv'), 'w')
        fout.write(header)
        
        fout_matched = open(os.path.join(output_dir, 'matched_sets.tsv'), 'w')
        fout_matched.write('Selected Bin ID\tCompleteness (%)\tContamination (%)\tQuality (%)\n')
        
        bin_num = 0
        total_comp = 0
        total_cont = 0
        total_quality = 0
        sel_gene_tables = {}
        sel_bins = {}
        selected_rows = {}
        while True:
            # determine highest quality match bin set
            bin_quality = self._bin_quality(bins, contigs, gene_tables, sel_weight)
            matched_sets = self._matched_bin_sets(bins, contigs, bin_quality, sel_min_quality, greedy)
            if len(matched_sets) == 0:
                break # no bins meeting selection criteria
                
            new_bin, num_removed, num_added = self._resolve_matched_set(matched_sets[0], 
                                                        bins,
                                                        contigs,
                                                        bin_quality,
                                                        remove_perc,
                                                        add_perc,
                                                        add_matches,
                                                        sel_min_quality,
                                                        greedy,
                                                        unanimous)
                
            domain, comp, cont = markers.bin_quality(new_bin)
            if comp - cont < 10:
                # skip extremely poor-quality bins
                break
                
            report_quality = comp - report_weight*cont
            total_comp += comp
            total_cont += cont
            total_quality += report_quality

            # report selection
            bin_num += 1
            unitem_bin_id = '%s_%d' % (self.bin_prefix, bin_num)
            primary_bm, primary_bid, _q, _n50, _gs = matched_sets[0][0]
            self.logger.info("Selected %s from %s with quality = %.1f (comp. = %.1f%%, cont. = %.1f%%)." % (
                                    primary_bid,
                                    primary_bm,
                                    report_quality,
                                    comp,
                                    cont))
                                    
            if not greedy and not unanimous:
                # performing consensus binning
                self.logger.info("-> Identified %d matched bins, removed %d contigs, added %d contigs." % (
                                        len(matched_sets[0]),
                                        num_removed, 
                                        num_added))
            elif not greedy:
                # performing unanimous binning
                self.logger.info("-> Identified %d matched bins and removed %d contigs." % (
                                        len(matched_sets[0]),
                                        num_removed))
                
            # write out matched set
            for i, (bm, bid, q, n50, gs) in enumerate(matched_sets[0]):
                domain, comp, cont = markers.bin_quality(bins[bm][bid])
                quality = comp - report_weight*cont
                if i != 0:
                    fout_matched.write('\t')
                fout_matched.write('%s\t%.1f\t%.1f\t%.1f' % (bm + '_' + bid, comp, cont, quality))
            fout_matched.write('\n')
            
            # write out highest quality bin
            new_bin_file = os.path.join(sel_bin_dir, unitem_bin_id + '.fna')
            seq_io.write_fasta(new_bin, new_bin_file)

            # write out summary information about selected bins
            row = self._report_selection(unitem_bin_id, 
                                            new_bin, 
                                            markers,
                                            report_weight,
                                            primary_bm,
                                            primary_bid,
                                            str(len(matched_sets[0])),
                                            str(num_removed),
                                            str(num_added),
                                            greedy,
                                            unanimous)
            fout.write(row)
            selected_rows[unitem_bin_id] = row

            # remove contigs in highest quality bin from marker gene tables and all other bins
            self._update_gene_tables(gene_tables, new_bin.keys())
            self._update_bins(bins, new_bin.keys())
            
            # cache selected gene tables and bin
            sel_gene_tables[unitem_bin_id] = markers.create_gene_table(new_bin.keys())  
            sel_bins[unitem_bin_id] = new_bin

        self.logger.info('Selected %d bins.' % bin_num)
        self.logger.info('-> total comp. = %.1f, total cont. = %.1f, total quality = %.1f' % (total_comp,
                                                                                            total_cont,
                                                                                            total_quality))

        fout.close()
        fout_matched.close()
        
        # merge partial bins
        merged_bids = set()
        merged_bins = {}
        if merge:
            merge = Merge(self.bin_prefix)
            merged_bids, merged_bins = merge.run(sel_bins, 
                                                    sel_gene_tables,
                                                    markers,
                                                    report_weight, 
                                                    report_min_quality, 
                                                    output_dir)
                                                    
            for bid in merged_bids:
                del sel_bins[bid]
                            
        # report final set of bins
        fout = open(os.path.join(output_dir, 'bins_final.tsv'), 'w')
        fout.write(header)
        num_reported_bins = 0
        total_comp = 0
        total_cont = 0
        total_quality = 0
        for bid, bin in sel_bins.iteritems():
            domain, comp, cont = markers.bin_quality(bin)
            quality = comp - report_weight*cont
            if quality < report_min_quality:
                continue
                
            num_reported_bins += 1
            total_comp += comp
            total_cont += cont
            total_quality += quality
                
            # copy over selected bin
            sel_bin = os.path.join(sel_bin_dir, bid + '.fna')
            final_bin = os.path.join(final_bin_dir, bid + '.fna')
            copyfile(sel_bin, final_bin)
            
            # report bin
            fout.write(selected_rows[bid])
            
        for merged_id, merged_bin in merged_bins.iteritems():
            domain, comp, cont = markers.bin_quality(merged_bin)
            quality = comp - report_weight*cont
            if quality < report_min_quality:
                continue
                
            num_reported_bins += 1
            total_comp += comp
            total_cont += cont
            total_quality += quality
            
            fout.write(self._report_selection(merged_id, 
                                                merged_bin, 
                                                markers,
                                                report_weight,
                                                'merged',
                                                'N/A',
                                                'N/A',
                                                'N/A',
                                                'N/A',
                                                greedy,
                                                unanimous))
        
        fout.close()
        
        self.logger.info('Identified %d bins passing reporting criteria.' % num_reported_bins)
        self.logger.info('-> total comp. = %.1f, total cont. = %.1f, total quality = %.1f' % (total_comp,
                                                                                            total_cont,
                                                                                            total_quality))
                                                                                            
        # summarize results of reported bins
        summary_file = os.path.join(output_dir, 'bins_final_quality_summary.tsv')
        self._bin_summary(sel_bins, 
                            merged_bins, 
                            report_weight, 
                            report_min_quality,
                            markers,
                            summary_file)
