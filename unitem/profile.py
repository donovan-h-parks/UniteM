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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import logging
from collections import defaultdict

from biolib.external.execute import check_on_path

from unitem.defaults import *


class Profile():
    """Profile genomes across different binning methods."""

    def __init__(self, cpus):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        check_on_path('checkm')
        
        self.cpus = cpus
        
    def _genome_quality(self, bac_quality_table, ar_quality_table):
        """Get CheckM estimates for each genome."""

        bac = {}
        with open(bac_quality_table) as f:
            f.readline()
            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[0]
                comp = float(line_split[5])
                cont = float(line_split[6])
                bac[gid] = (comp, cont)
                
        ar = {}
        with open(ar_quality_table) as f:
            f.readline()
            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[0]
                comp = float(line_split[5])
                cont = float(line_split[6])
                ar[gid] = (comp, cont)
        
        gq = {}
        for gid in set(bac.keys()).union(ar):
            bac_comp, bac_cont = bac[gid]
            ar_comp, ar_cont = ar[gid]
            if bac_comp + bac_cont > ar_comp + ar_cont:
                gq[gid] = ('Bacteria', bac_comp, bac_cont)
            else:
                gq[gid] = ('Archaea', ar_comp, ar_cont)
                
        return gq
        
    def _report_genome_quality(self, genome_quality, output_dir):
        """Summarize quality of genomes."""
        
        # create table for each binning method
        for bm in genome_quality:
            table = os.path.join(output_dir, bm + '_quality.tsv')
            fout = open(table, 'w')
            fout.write('Genome ID\tMarker Set Domain\tCompleteness (%)\tContamination (%)\tQuality\n')
            for gid in genome_quality[bm]:
                domain, comp, cont = genome_quality[bm][gid]
                fout.write('%s\t%s\t%.2f\t%.2f\t%.2f\n' % (gid, domain, comp, cont, comp-5*cont))
            fout.close()
            
        # report global results file
        summary = {}
        for bm in genome_quality:
            comp_cont_count = defaultdict(lambda: defaultdict(int))
            quality_count = defaultdict(int)
            total_comp = 0
            total_cont = 0
            total_q = 0
            for bid, (domain, comp, cont) in genome_quality[bm].items():
                quality = comp - 5*cont
                
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
                
                summary[bm] = (comp_cont_count, quality_count, total_comp, total_cont, total_q)
        
        fout = open(os.path.join(output_dir, 'bin_quality_summary.tsv'), 'w')
        fout.write('\t' + '\t'.join(sorted(summary)) + '\n')
        
        for comp in [90, 80, 70]:
            for cont in [5, 10]:
                fout.write('Comp. >= %d, cont. <= %d' % (comp, cont))
                for bm in sorted(summary):
                    fout.write('\t%d' % summary[bm][0][comp][cont])
                fout.write('\n')

        for q in [90, 70, 50]:
            fout.write('Quality >= %d' % q)
            for bm in sorted(summary):
                fout.write('\t%d' % summary[bm][1][q])
            fout.write('\n')
            
        fout.write('Total compl.')
        for bm in sorted(summary):
            fout.write('\t%.1f' % summary[bm][2])
        fout.write('\n')
            
        fout.write('Total cont.')
        for bm in sorted(summary):
            fout.write('\t%.1f' % summary[bm][3])
        fout.write('\n')
        
        fout.write('Total quality')
        for bm in sorted(summary):
            fout.write('\t%.1f' % summary[bm][4])
        fout.write('\n')

        fout.close()

    def run(self, bin_dirs, output_dir):
        """Profile genomes in each bin directory.

        Parameters
        ----------
        bin_dirs : list of str
            Directories containing bins from different binning methods.
        output_dir : str
            Output directory.
        """

        self.logger.info('Profiling genomes in %d directories.' % len(bin_dirs))
        
        num_processed = 0
        genome_quality = defaultdict(lambda: dict)
        for method_id, (bin_dir, bin_ext) in bin_dirs.items():
            num_processed += 1
            self.logger.info('Profiling %s (%d of %d).' % (method_id, num_processed, len(bin_dirs)))
            
            for d, ms_file in [(CHECKM_BAC_DIR, CHECKM_BAC_MS), (CHECKM_AR_DIR, CHECKM_AR_MS)]:
                cur_output_dir = os.path.join(output_dir, BINNING_METHOD_DIR, method_id, d)
                cmd = 'checkm analyze -t %d -x %s %s %s %s' % (self.cpus, 
                                                                bin_ext, 
                                                                ms_file, 
                                                                bin_dir, 
                                                                cur_output_dir)
                os.system(cmd)
                
                marker_gene_table = os.path.join(cur_output_dir, MARKER_GENE_TABLE)
                cmd = 'checkm qa -t %d -o 5 --tab_table -f %s %s %s' % (self.cpus,
                                                                        marker_gene_table, 
                                                                        ms_file, 
                                                                        cur_output_dir)
                os.system(cmd)
                
                genome_quality_table = os.path.join(cur_output_dir, GENOME_QUALITY_TABLE)
                cmd = 'checkm qa -t %d -o 2 --tab_table -f %s %s %s' % (self.cpus,
                                                                        genome_quality_table, 
                                                                        ms_file, 
                                                                        cur_output_dir)
                os.system(cmd)
                
            bac_quality_table = os.path.join(output_dir, 
                                                BINNING_METHOD_DIR, 
                                                method_id, 
                                                CHECKM_BAC_DIR, 
                                                GENOME_QUALITY_TABLE)
            ar_quality_table = os.path.join(output_dir, 
                                                BINNING_METHOD_DIR, 
                                                method_id, 
                                                CHECKM_AR_DIR, 
                                                GENOME_QUALITY_TABLE)
            if not os.path.exists(bac_quality_table) or not os.path.exists(ar_quality_table):
                self.logger.error('Missing quality table for %s.' % method_id)
                self.logger.error('Please verify there were bins in the bin directory specified for this method.')
                sys.exit()

            genome_quality[method_id] = self._genome_quality(bac_quality_table,
                                                                ar_quality_table)

        self._report_genome_quality(genome_quality, output_dir)