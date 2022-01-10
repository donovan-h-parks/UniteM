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
from collections import Counter

from unitem.utils import (make_sure_path_exists,
                          check_dir_exists,
                          check_file_exists)
from unitem.bin_tools import BinTools
from unitem.bin import Bin
from unitem.profile import Profile
from unitem.ensemble import Ensemble


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def _bin_extension(self, bin_dir):
        """Determine extension of bins."""

        exts = []
        for f in os.listdir(bin_dir):
            f_split = f.split('.')

            if len(f_split) > 1:
                ext = f_split[-1]
                if ext in ['fa', 'fna', 'fasta']:
                    exts.append(ext)

        if len(exts) == 0:
            return None, None

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
                if not bin_ext:
                    self.logger.warning(
                        f'No bins identified for {method_id} in {d}.')
                else:
                    bin_dirs[method_id] = (d, bin_ext)
                    self.logger.info(
                        f"Processing {count} genomes from {method_id} with extension '{bin_ext}'.")

        if hasattr(options, 'bin_file') and options.bin_file:
            check_file_exists(options.bin_file)
            for line in open(options.bin_file):
                if line.strip():
                    tokens = [token.strip() for token in line.split('\t')]
                    if len(tokens) != 2:
                        self.logger.warning(
                            f"Skipping invalid line: {line.strip()}")
                        continue

                    method_id = tokens[0]
                    d = tokens[1]
                    check_dir_exists(d)
                    bin_ext, count = self._bin_extension(d)
                    if not bin_ext:
                        self.logger.warning(
                            f'No bins identified for {method_id} in {d}.')
                    else:
                        bin_dirs[method_id] = (d, bin_ext)
                        self.logger.info(
                            f"Processing {count} genomes from {method_id} with extension '{bin_ext}'.")

        return bin_dirs

    def bin(self, options):
        """Bin command"""

        check_file_exists(options.assembly_file)
        make_sure_path_exists(options.output_dir)

        bin = Bin(options.assembly_file,
                  options.output_dir,
                  options.min_contig_len,
                  options.cpus)
        bin.check_on_path(options)
        bin.coverage(options.bam_files, options.cov_file)
        bin.run(options)

        self.logger.info(
            f"UniteM 'bin' results written to: {options.output_dir}")

    def profile(self, options):
        """Profile command"""

        make_sure_path_exists(options.output_dir)

        bin_dirs = self._bin_dirs(options)

        profile = Profile(options.cpus)
        profile.run(bin_dirs,
                    options.marker_dir,
                    options.keep_intermediate,
                    options.output_dir)

        self.logger.info(
            f"UniteM 'profile' results written to: {options.output_dir}")

    def consensus(self, options):
        """Consensus command"""

        check_dir_exists(options.profile_dir)
        make_sure_path_exists(options.output_dir)

        bin_dirs = self._bin_dirs(options)

        e = Ensemble(options.bin_prefix)
        e.run(options.profile_dir,
              bin_dirs,
              options.marker_dir,
              options.weight,
              options.sel_min_quality,
              options.sel_min_comp,
              options.sel_max_cont,
              options.remove_perc,
              options.add_perc,
              options.add_matches,
              False,  # perform consensus bin selection
              False,  # perform unanimous bin selection
              options.report_min_quality,
              options.simple_headers,
              options.output_dir)

        self.logger.info(
            f"UniteM 'consensus' results written to: {options.output_dir}")

    def greedy(self, options):
        """Greedy command"""

        check_dir_exists(options.profile_dir)
        make_sure_path_exists(options.output_dir)

        bin_dirs = self._bin_dirs(options)

        e = Ensemble(options.bin_prefix)
        e.run(options.profile_dir,
              bin_dirs,
              options.marker_dir,
              options.weight,
              options.sel_min_quality,
              options.sel_min_comp,
              options.sel_max_cont,
              None,
              None,
              None,
              True,   # perform greedy bin selection
              False,  # perform unanimous bin selection
              options.report_min_quality,
              options.simple_headers,
              options.output_dir)

        self.logger.info(
            f"UniteM 'greedy' results written to: {options.output_dir}")

    def unanimous(self, options):
        """Unanimous command"""

        check_dir_exists(options.profile_dir)
        make_sure_path_exists(options.output_dir)

        bin_dirs = self._bin_dirs(options)

        e = Ensemble(options.bin_prefix)
        e.run(options.profile_dir,
              bin_dirs,
              options.marker_dir,
              options.weight,
              options.sel_min_quality,
              options.sel_min_comp,
              options.sel_max_cont,
              None,
              None,
              None,
              False,  # perform greedy bin selection
              True,   # perform unanimous bin selection
              options.report_min_quality,
              options.simple_headers,
              options.output_dir)

        self.logger.info(
            f"UniteM 'unanimous' results written to: {options.output_dir}")

    def compare(self, options):
        """Compare command"""

        check_dir_exists(options.bin_dir1)
        check_dir_exists(options.bin_dir2)

        bt = BinTools()
        bin_files1 = bt.bin_files(options.bin_dir1, options.extension1)
        bin_files2 = bt.bin_files(options.bin_dir2, options.extension2)

        bt.compare(bin_files1, bin_files2,
                   options.assembly_file, options.output_file)

        self.logger.info("UniteM 'compare' results written to: %s" %
                         options.output_file)

    def unique(self, options):
        """Unique command"""

        check_dir_exists(options.bin_dir)

        bt = BinTools()
        bin_files = bt.bin_files(options.bin_dir, options.extension)
        bt.unique(bin_files)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if(options.subparser_name == 'bin'):
            self.bin(options)
        elif(options.subparser_name == 'profile'):
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
            self.logger.error(
                f'Unknown UniteM command: {options.subparser_name}\n')
            sys.exit(1)

        return 0
