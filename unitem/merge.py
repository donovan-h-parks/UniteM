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
import ast
import logging
import itertools

from biolib.common import find_nearest, alphanumeric_sort
import biolib.seq_io as seq_io

from unitem.common import bin_gc


class Merge():
    """Identify bins with complementary marker sets."""
    
    def __init__(self, bin_prefix):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.bin_prefix = bin_prefix
        
        self.min_delta_comp = 20.0
        self.max_delta_cont = 0.0
        self.min_merged_comp = 70.0
        self.max_merged_cont = 10.0
        
        self.gc_dist = None
        
    def _read_distribution(self, prefix):
        """Read distribution file.

        Parameters
        ----------
        prefix : str
            Prefix of distibution to read (gc_dist or td_dist)

        Returns
        -------
        dict : d[percentile] -> critical value
            Critical values at integer percentiles.
        """

        dist_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                    'distributions', prefix + '.txt')

        with open(dist_file, 'r') as f:
            s = f.read()
            d = ast.literal_eval(s)

        return d
        
    def check_gc(self, bin1, bin2):
        """Check if bins have compatible GC content."""
        
        # make sure distributions have been loaded
        if not self.gc_dist:
            self.gc_dist = self._read_distribution('gc_dist')
            
        gc1 = bin_gc(bin1)
        gc2 = bin_gc(bin2)
        
        bin_size1 = sum([len(s) for s in bin1.values()])
        bin_size2 = sum([len(s) for s in bin2.values()])
        
        # perform test based on bin sizes
        if bin_size1 < bin_size2:
            bin_size = bin_size1
            test_gc = gc2
            delta_gc = gc1 - gc2
        else:
            bin_size = bin_size2
            test_gc = gc1
            delta_gc = gc2 - gc1
        
        # find keys into GC distributions
        # gc -> [mean GC][scaffold length][percentile]
        closest_gc = find_nearest(self.gc_dist.keys(), test_gc)
        closest_len = find_nearest(self.gc_dist[closest_gc].keys(), bin_size)
        
        gc_per = 95
        gc_lower_bound_key = find_nearest(self.gc_dist[closest_gc][closest_len].keys(), 
                                            (100 - gc_per) / 2.0)
        gc_upper_bound_key = find_nearest(self.gc_dist[closest_gc][closest_len].keys(), 
                                            (100 + gc_per) / 2.0)

        gc_lower_bound = self.gc_dist[closest_gc][closest_len][gc_lower_bound_key]
        gc_upper_bound = self.gc_dist[closest_gc][closest_len][gc_upper_bound_key]
        if delta_gc < gc_lower_bound or delta_gc > gc_upper_bound:
            return False, gc1, gc2
        
        return True, gc1, gc2

    def run(self, 
            bins, 
            gene_tables,
            markers,
            report_weight, 
            report_min_quality,            
            output_dir):
        """Perform consensus selection of genomes across multiple binning methods.

        Parameters
        ----------
        bins : d[bin ID][contig Id] -> seq
          Contigs for bins across all binning methods.
        gene_tables : d[bin Id] -> (bac table, ar table)
          Bacterial and archaeal marker gene tables for bin.
        output_dir : str
          Output directory.
        """
        
        self.logger.info('Identifying pairs of bins with complementary marker sets.')

        # determine union and intersection of marker sets for each pair of bins
        output_file = os.path.join(output_dir, "bins_merged.tsv")
        fout = open(output_file, 'w')
        fout.write('UniteM Merged Id')
        fout.write('\tUniteM Id 1\tUniteM Id 2\tMarker Domain')
        fout.write('\tBin 1 GC\tBin 2 GC')
        fout.write('\tBin 1 Completeness\tBin 1 Contamination')
        fout.write('\tBin 2 Completeness\tBin 2 Contamination')
        fout.write('\tDelta Completeness\tDelta Contamination')
        fout.write('\tMerged Completeness\tMerged Contamination\n')
        
        # re-organize gene tables for easier iteration
        flat_gene_tables = []
        for bid, (bac_table, ar_table) in gene_tables.iteritems():
            flat_gene_tables.append((bid, bac_table, ar_table))
             
        # compare all pairs of bins
        merged_bids = set()
        merged_bins = {}
        final_bin_dir = os.path.join(output_dir, 'bins_final')
        for d1, d2 in itertools.combinations(flat_gene_tables, 2):
            bid1, bac1, ar1 = d1
            bid2, bac2, ar2 = d2
            
            # put into sorted order by bin ID
            if alphanumeric_sort([bid1, bid2])[0] == bid2:
                bid1, bid2 = bid2, bid1
                bac1, bac2 = bac2, bac1
                ar1, ar2 = ar2, ar1
            
            if bid1 in merged_bids or bid2 in merged_bids:
                continue
            
            domain1, comp1, cont1, = markers.evaluate(bac1, ar1)
            domain2, comp2, cont2, = markers.evaluate(bac2, ar2)
            if domain1 != domain2:
                continue

            # evaluate merge bins
            bac_merge = {}
            for mid in set(bac1.keys()).union(bac2.keys()):
                bac_merge[mid] = bac1.get(mid, [])+ bac2.get(mid, [])
            
            ar_merge = {}
            for mid in set(ar1.keys()).union(ar2.keys()):
                ar_merge[mid] = ar1.get(mid, [])+ ar2.get(mid, [])
                
            merge_domain, merge_comp, merge_cont = markers.evaluate(bac_merge, ar_merge)
            if merge_domain != domain1:
                continue

            if merge_comp < self.min_merged_comp or merge_cont > self.max_merged_cont:
                continue
                
            delta_comp = merge_comp - max(comp1, comp2)
            delta_cont = merge_cont - cont1 - cont2

            if delta_comp < self.min_delta_comp or delta_cont > self.max_delta_cont:
                continue

            pass_gc_test, gc1, gc2 = self.check_gc(bins[bid1], bins[bid2])
            if not pass_gc_test:
                continue
                
            # merge bin
            merged_bin = bins[bid1].copy()
            merged_bin.update(bins[bid2])
            domain, comp, cont = markers.bin_quality(merged_bin)
            quality = comp - report_weight*cont
            if quality < report_min_quality:
                continue
                
            bin_num1 = int(bid1.split('_')[-1])
            bin_num2 = int(bid2.split('_')[-1])
            merged_bin_name = '%s_merged_%d_%d' % (self.bin_prefix, bin_num1, bin_num2)
                
            merged_bids.add(bid1)
            merged_bids.add(bid2)
            merged_bins[merged_bin_name] = merged_bin
            
            self.logger.info("Merged %s and %s with quality = %.1f (comp. = %.1f%%, cont. = %.1f%%)." % (
                                    bid1,
                                    bid2,
                                    quality,
                                    comp,
                                    cont))

            # report merger
            fout.write('%s' % merged_bin_name)
            fout.write('\t%s\t%s\t%s' % (bid1, bid2, merge_domain))
            fout.write('\t%.2f\t%.2f' % (gc1*100, gc2*100))
            fout.write('\t%.2f\t%.2f' % (comp1, cont1))
            fout.write('\t%.2f\t%.2f' % (comp2, cont2))
            fout.write('\t%.2f\t%.2f' % (delta_comp, delta_cont))
            fout.write('\t%.2f\t%.2f' % (merge_comp, merge_cont))
            fout.write('\n')
            
            # write out bin
            merged_bin_file = os.path.join(final_bin_dir, merged_bin_name + '.fna')
            seq_io.write_fasta(merged_bin, merged_bin_file)
            
        fout.close()
        
        self.logger.info('Merged %d bins.' % len(merged_bins))
        
        return merged_bids, merged_bins
