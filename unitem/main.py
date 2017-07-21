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
from collections import defaultdict, Counter

from unitem.bin_tools import BinTools
from unitem.profile import Profile
from unitem.consensus import Consensus

from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists,
                           query_yes_no,
                           remove_extension)


class OptionsParser():
    def __init__(self):
        """Initialization"""
        
        self.logger = logging.getLogger('timestamp')
        self.reporter = logging.getLogger('no_timestamp')
        
    def _bin_extension(self, bin_dir):
        """Determine extension of bins."""
        
        exts = []
        for f in os.listdir(bin_dir):
            f_split = f.split('.')
            
            if len(f_split) > 1:
                exts.append(f_split[-1])
                
        ext, count = Counter(exts).most_common(1)[0]
        
        return ext, count
        
    def _bin_dirs(self, options):
        """Get directories with bins from different binners."""
        
        bin_dirs = {}
        if hasattr(options, 'bin_dirs') and options.bin_dirs:
            for d in options.bin_dirs:
                check_dir_exists(d)
                method_id = os.path.basename(os.path.normpath(d))
                bin_ext, count = self._bin_extension(d)
                bin_dirs[method_id] = (d, bin_ext)
                
                self.logger.info("Processing %d genomes from %s with extension '%s'." % (count, 
                                                                                            method_id, 
                                                                                            bin_ext))
                
        if hasattr(options, 'bin_file') and options.bin_file:
            for line in open(options.bin_file):
                if line.strip():
                    line_split = map(str.strip, line.split('\t'))
                    if len(line_split) != 2:
                        self.logger.warning("Skipping invalid line: %s" % line.strip())
                        continue
                        
                    method_id = line_split[0]
                    d = line_split[1]
                    check_dir_exists(d)
                    bin_ext, count = self._bin_extension(d)
                    bin_dirs[method_id] = (d, bin_ext)
                    
                    self.logger.info("Processing %d genomes from %s with extension '%s'." % (count, 
                                                                                                method_id, 
                                                                                                bin_ext))
                    
        if len(bin_dirs) < 2:
            self.logger.error("Two or more bin directories must be specified.")
            sys.exit(-1)
                    
        return bin_dirs
                    
    def profile(self, options):
        """Profile command"""
        
        make_sure_path_exists(options.output_dir)
        
        bin_dirs = self._bin_dirs(options)
        
        profile = Profile(options.cpus)
        profile.run(bin_dirs, 
                    options.output_dir)

        self.logger.info("UniteM 'profile' results written to: %s" % options.output_dir)

            
    def consensus(self, options):
        """Consensus command"""
        
        check_dir_exists(options.profile_dir)
        make_sure_path_exists(options.output_dir)
        
        bin_dirs = self._bin_dirs(options)
        
        consensus = Consensus(options.bin_prefix)
        consensus.run(options.profile_dir,
                        bin_dirs,
                        options.sel_weight,
                        options.sel_min_quality,
                        options.remove_perc,
                        options.add_perc,
                        options.add_matches,
                        False,  # perform consensus bin selection
                        False,  # perform unanimous bin selection
                        False,  # do not perform merging
                        options.report_weight,
                        options.report_min_quality,
                        options.output_dir)
        
        self.logger.info("UniteM 'consensus' results written to: %s" % options.output_dir)
        
    def greedy(self, options):
        """Greedy command"""
        
        check_dir_exists(options.profile_dir)
        make_sure_path_exists(options.output_dir)
        
        bin_dirs = self._bin_dirs(options)
        
        consensus = Consensus(options.bin_prefix)
        consensus.run(options.profile_dir,
                        bin_dirs,
                        options.sel_weight,
                        options.sel_min_quality,
                        None,
                        None,
                        None,
                        True,   # perform greedy bin selection
                        False,  # perform unanimous bin selection
                        False,  # do not perform merging
                        options.report_weight,
                        options.report_min_quality,
                        options.output_dir)
        
        self.logger.info("UniteM 'greedy' results written to: %s" % options.output_dir)
        
    def unanimous(self, options):
        """Unanimous command"""
        
        check_dir_exists(options.profile_dir)
        make_sure_path_exists(options.output_dir)
        
        bin_dirs = self._bin_dirs(options)
        
        consensus = Consensus(options.bin_prefix)
        consensus.run(options.profile_dir,
                        bin_dirs,
                        options.sel_weight,
                        options.sel_min_quality,
                        None,
                        None,
                        None,
                        False,  # perform greedy bin selection
                        True,   # perform unanimous bin selection
                        False,  # do not perform merging
                        options.report_weight,
                        options.report_min_quality,
                        options.output_dir)
        
        self.logger.info("UniteM 'greedy' results written to: %s" % options.output_dir)
    
    def compare(self, options):
        """Compare command"""

        check_dir_exists(options.bin_dir1)
        check_dir_exists(options.bin_dir2)

        bt = BinTools()
        bin_files1 = bt.bin_files(options.bin_dir1, options.extension1)
        bin_files2 = bt.bin_files(options.bin_dir2, options.extension2)

        bt.compare(bin_files1, bin_files2, options.assembly_file, options.output_file)

        self.logger.info("UniteM 'compare' results written to: %s" % options.output_file)
        
    def unique(self, options):
        """Unique command"""

        check_dir_exists(options.bin_dir)
        
        bt = BinTools()
        bin_files = bt.bin_files(options.bin_dir, options.extension)
        bt.unique(bin_files)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if(options.subparser_name == 'profile'):
            self.profile(options)
        elif(options.subparser_name == 'consensus'):
            self.consensus(options)
        elif(options.subparser_name == 'greedy'):
            self.greedy(options)    
        elif(options.subparser_name == 'unanimous'):
            self.unanimous(options)   
        elif(options.subparser_name == 'compare'):
            self.compare(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        else:
            self.logger.error('Unknown UniteM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
