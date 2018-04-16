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

file_dir = os.path.dirname(os.path.realpath(__file__))

CHECKM_BAC_MS = os.path.join(file_dir, 'checkm_ms', 'bacteria.ms')
CHECKM_AR_MS = os.path.join(file_dir, 'checkm_ms', 'archaea.ms')

CHECKM_BAC_DIR = 'checkm_bac'
CHECKM_AR_DIR = 'checkm_ar'
BINNING_METHOD_DIR = 'binning_methods'

MARKER_GENE_TABLE = 'marker_gene_table.tsv'
GENOME_QUALITY_TABLE = 'genome_quality.tsv'