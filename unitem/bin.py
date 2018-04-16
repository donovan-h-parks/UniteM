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
import ntpath
import logging
from collections import defaultdict

from biolib.external.execute import check_on_path, run
from biolib.common import make_sure_path_exists


class Bin():
    """Apply binning methods to an assembly."""

    def __init__(self, assembly_file, output_dir, min_contig_len, cpus):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.assembly_file = assembly_file
        self.output_dir = output_dir
        self.min_contig_len = min_contig_len
        self.cpus = cpus
        
        self.failed_methods = []
        
    def coverage(self, bam_files, cov_file):
        """Calculate coverage file for use by different binning methods."""
        
        self.bam_files = bam_files
        
        if cov_file:
            self.logger.info('Using coverage information in %s.' % cov_file)
            # check coverage file has correct format
            header = open(cov_file).readline().split('\t')
            if header[0] != 'contigName' or header[1] != 'contigLen' or header[2] != 'totalAvgDepth':
                self.logger.error('Provided coverage file does not have the correct headers.')
                self.logger.error("Coverage file must have the format produced by 'jgi_summarize_bam_contig_depths'.")
                sys.exit(1)
            self.cov_file = cov_file
            
            if bam_files:
                self.logger.warning('BAM files are being ignored.')
        else:
            found = check_on_path('jgi_summarize_bam_contig_depths', exit_on_fail=False)
            if not found:
                self.logger.error('jgi_summarize_bam_contig_depths is not on the system path.')
                self.logger.error('This script is provide with MetaBAT v2.')
                sys.exit(1)

            self.logger.info('Calculating coverage for %d BAM files.' % len(bam_files))
            self.logger.info("Running 'jgi_summarize_bam_contig_depths'.")
            
            self.cov_file = os.path.join(self.output_dir, 'coverage.tsv')
            cmd = 'jgi_summarize_bam_contig_depths --minContigLength %d --minContigDepth 1 --outputDepth %s %s' % (self.min_contig_len,
                                                                                                                    self.cov_file, 
                                                                                                                    ' '.join(bam_files))

            success, exception = run(cmd)
            if not success:
                self.logger.error('Failed to execute: %s' % ' '.join(exception.cmd))
                sys.exit(1)
        
    def check_on_path(self, options):
        """Check that all binning methods are on the system path."""
        
        if options.mb2:
            self.logger.info('Checking MetaBAT v2 dependencies.')
            check_on_path('metabat2')
        if options.gm2:
            self.logger.info('Checking GroopM v2 dependencies.')
            check_on_path('groopm2')
        if options.max40 or options.max107:
            self.logger.info('Checking MaxBin dependencies.')
            check_on_path('run_MaxBin.pl')
            check_on_path('FragGeneScan')
        if options.bs:
            self.logger.info('Checking BinSanity dependencies.')
            check_on_path('Binsanity-wf')
            check_on_path('transform-coverage-profile')
            check_on_path('checkm')
        if (options.mb_verysensitive 
                or options.mb_sensitive
                or options.mb_specific
                or options.mb_veryspecific
                or options.mb_superspecific):
            self.logger.info('Checking MetaBAT dependencies.')
            check_on_path('metabat1')

    def run(self, options):
        """Run binning methods."""
        
        bin_file = os.path.join(self.output_dir, 'bin_dirs.tsv')
        bin_file_out = open(bin_file, 'w')
        
        if options.mb2:
            self.metabat2(bin_file_out)
        if options.gm2:
            self.groopm2(bin_file_out)
        if options.max40:
            self.maxbin(bin_file_out, 40)
        if options.max107:
            self.maxbin(bin_file_out, 107)
        if options.bs:
            self.binsanity(bin_file_out)
        if options.mb_verysensitive:
            self.metabat(bin_file_out, 'verysensitive')
        if options.mb_sensitive:
            self.metabat(bin_file_out, 'sensitive')
        if options.mb_specific:
            self.metabat(bin_file_out, 'specific')
        if options.mb_veryspecific:
            self.metabat(bin_file_out, 'veryspecific')
        if options.mb_superspecific:
            self.metabat(bin_file_out, 'superspecific')
            
        bin_file_out.close()

        self.logger.info('File with location of bin directories written to %s.' % bin_file)
        
        if self.failed_methods:
            self.logger.warning('The following methods failed to run: %s' % ' '.join(self.failed_methods))
            
    def _run_method(self, cmd, bin_dir, bin_file_out, binning_method):
        """Run binning method."""
        
        success, output = run(cmd)
        if not success:
            self.logger.warning('Failed to execute: %s' % ' '.join(output.cmd))
            self.failed_methods.append(binning_method)
            return
            
        fout = open(os.path.join(bin_dir, 'output.log'), 'w')
        fout.write(output)
        fout.close()

        bin_file_out.write('%s\t%s\n' % (binning_method, os.path.abspath(bin_dir)))
    
    def metabat2(self, bin_file_out):
        """Run MetaBAT v2."""
        
        self.logger.info("Running 'metabat2'.")
        bin_dir = os.path.join(self.output_dir, 'metabat2')
        bin_prefix = os.path.join(bin_dir, 'mb2')
        cmd = 'metabat2 -t %d -m %d -i %s -a %s -o %s' % (self.cpus,
                                                            self.min_contig_len,
                                                            self.assembly_file, 
                                                            self.cov_file, 
                                                            bin_prefix)
        
        self._run_method(cmd, bin_dir, bin_file_out, 'metabat2')

    def metabat(self, bin_file_out, preset):
        """Run MetaBAT."""
        
        self.logger.info("Running 'metabat' with the %s preset." % preset)
        bin_dir = os.path.join(self.output_dir, 'metabat_%s' % preset)
        bin_prefix = os.path.join(bin_dir, 'mb_%s' % preset)
        cmd = 'metabat1 -t %d -m %d -i %s -a %s -o %s --%s' % (self.cpus,
                                                            self.min_contig_len,
                                                            self.assembly_file, 
                                                            self.cov_file, 
                                                            bin_prefix,
                                                            preset)
                                                            
        self._run_method(cmd, bin_dir, bin_file_out, 'metabat_%s' % preset)
        
    def groopm2(self, bin_file_out):
        """Run GroopM v2."""

        self.logger.info("Running 'groopm2 parse'.")
        bin_dir = os.path.join(self.output_dir, 'groopm2')
        make_sure_path_exists(bin_dir)
        output_db = os.path.join(bin_dir, 'groopm.db')
        cmd = 'groopm2 parse -f -t %d -c %d --cov_file %s %s %s' % (self.cpus,
                                                                    self.min_contig_len,
                                                                    self.cov_file,
                                                                    output_db, 
                                                                    self.assembly_file)
        success, exception = run(cmd)
        if not success:
            self.logger.warning('Failed to execute: %s' % ' '.join(exception.cmd))
            self.failed_methods.append('GroopM v2')
            return

        self.logger.info("Running 'groopm2 core'.")
        cmd = 'groopm2 core -f %s -c %d --save_dists' % (output_db,
                                                        self.min_contig_len)
        success, exception = run(cmd)
        if not success:
            self.logger.warning('Failed to execute: %s' % ' '.join(exception.cmd))
            self.failed_methods.append('GroopM v2')
            return

        self.logger.info("Running 'groopm2 extract'.")
        bin_prefix = os.path.join(bin_dir, 'gm2')
        cmd = 'groopm2 extract -p %s %s %s' % (bin_prefix, 
                                                output_db, 
                                                self.assembly_file)
                                                
        self._run_method(cmd, bin_dir, bin_file_out, 'groopm2')
        
    def _create_maxbin_coverage_files(self, cov_file, output_dir):
        """Parse coverage information files required by MaxBin."""

        abund_list_file = os.path.join(output_dir, 'abund_files.lst')
        fout = open(abund_list_file, 'w')
        with open(cov_file) as f:
            headers = f.readline().rstrip().split('\t')
            bam_headers = headers[3::2]

            fhs = []
            for bh in bam_headers:
                abund_file = os.path.abspath(os.path.join(output_dir, bh + '.abund.tsv'))
                fhs.append(open(abund_file, 'w'))
                fout.write(abund_file + '\n')

            for line in f:
                line_split = line.rstrip().split('\t')
                contig_id = line_split[0]
                for fh_index, col_index in enumerate(range(3,len(line_split),2)):
                    fhs[fh_index].write(contig_id + '\t' + line_split[col_index] + '\n')

            for fh in fhs:
                fh.close()
        fout.close()

        return abund_list_file
        
    def maxbin(self, bin_file_out, num_markers):
        """Run MaxBin."""
        
        bin_dir = os.path.join(self.output_dir, 'maxbin_ms%d' % num_markers)
        make_sure_path_exists(bin_dir)
        cov_file_dir = os.path.join(bin_dir, 'coverage_files')
        make_sure_path_exists(cov_file_dir)
        abund_list_file = self._create_maxbin_coverage_files(self.cov_file, cov_file_dir)
        
        self.logger.info("Running 'run_MaxBin.pl' with %d markers." % num_markers)
        bin_prefix = os.path.join(bin_dir, 'max%d' % num_markers)
        cmd = 'run_MaxBin.pl -min_contig_length %d' % self.min_contig_len
        cmd += ' -thread %d -markerset %d -contig %s -out %s -abund_list %s' % (self.cpus,
                                                                                    num_markers,
                                                                                    self.assembly_file,
                                                                                    bin_prefix,
                                                                                    abund_list_file)
                                                                                    
        self._run_method(cmd, bin_dir, bin_file_out, 'maxbin_ms%d' % num_markers)
        
    def _create_binsanity_cov_file(self, cov_file, output_dir):
        """Create coverage file required by BinSanity."""
        
        bs_cov_file = os.path.join(output_dir, 'binsanity.cov')
        fout = open(bs_cov_file, 'w')
        with open(cov_file) as f:
            f.readline()
            
            for line in f:
                line_split = line.split('\t')
                fout.write(line_split[0])
                fout.write('\t' + '\t'.join(line_split[3::2]) + '\n')
        fout.close()
        
        bs_cov_scale_file = os.path.join(output_dir, 'binsanity.cov.x100.lognorm')
        cmd = 'transform-coverage-profile -t scale -c %s' % bs_cov_file
        os.system(cmd)
        
        return bs_cov_scale_file

    def binsanity(self, bin_file_out):
        """Run BinSanity."""
        
        self.logger.info("Running 'Binsanity transform-coverage-profile'.")
        bin_dir = os.path.join(self.output_dir, 'binsanity')
        make_sure_path_exists(bin_dir)
        bs_cov_file = self._create_binsanity_cov_file(self.cov_file, bin_dir)
        
        self.logger.info("Running 'Binsanity-wf'.")
        assembly_dir, assembly_filename = ntpath.split(self.assembly_file)
        cmd = 'Binsanity-wf --threads %d -x %d -l %s -c %s -o %s --Prefix %s' % (self.cpus,
                                                                        self.min_contig_len,
                                                                        assembly_filename,
                                                                        bs_cov_file,
                                                                        bin_dir,
                                                                        os.path.join(bin_dir, 'out'))
        if assembly_dir:
            cmd += ' -f %s' % assembly_dir
        else:
            cmd += ' -f %s' % '.'
        
        self._run_method(cmd, bin_dir, bin_file_out, 'binsanity')
        
        # ***HACK: BinSanity insists on putting bins in BinSanity-Final-bins
        os.system('mv %s %s' % (os.path.join(bin_dir, 'BinSanity-Final-bins', '*'), bin_dir))
