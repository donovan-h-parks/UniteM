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
import logging
from collections import defaultdict

from unitem.defaults import (TIGRFAM_TOP_HIT_SUFFIX,
                             PFAM_TOP_HIT_SUFFIX,
                             BINNING_METHOD_DIR)
from unitem.external.tophits import TopHitFile
from unitem.external.tigrfam_search import TigrfamSearch
from unitem.external.pfam_search import PfamSearch


class Markers():
    """Read marker information."""

    def __init__(self, marker_dir):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        pfam_bac_ms, pfam_ar_ms = PfamSearch(marker_dir, 1).get_marker_genes()
        tigr_bac_ms, tigr_ar_ms = TigrfamSearch(
            marker_dir, 1).get_marker_genes()

        self.bac_ms = pfam_bac_ms.union(tigr_bac_ms)
        self.ar_ms = pfam_ar_ms.union(tigr_ar_ms)

    @staticmethod
    def estimate_comp_cont(hmm_hits, markers):
        """Estimate genome completeness and contamination.

        completeness = (# identified genes)/(total marker genes)
        contamination = (multi-copy count)/(total marker genes)
        """

        present = 0
        multi_copy_count = 0
        for hmm_id, hits in hmm_hits.items():
            if hmm_id in markers:
                present += 1
                multi_copy_count += len(hits) - 1

        if len(markers) > 0:
            perc_comp = 100 * float(present) / len(markers)
            perc_cont = 100 * float(multi_copy_count) / len(markers)
        else:
            perc_comp = perc_cont = 0

        return perc_comp, perc_cont

    @staticmethod
    def read_top_hit_tables(top_hit_dir):
        """Read top hit tables for all bins in directory."""

        markers = defaultdict(lambda: defaultdict(list))

        for f in os.listdir(top_hit_dir):
            if not (f.endswith(PFAM_TOP_HIT_SUFFIX) or f.endswith(TIGRFAM_TOP_HIT_SUFFIX)):
                continue

            top_hit_file = os.path.join(top_hit_dir, f)
            top_hits = TopHitFile(top_hit_file)
            top_hits.read()

            bin_id = f.replace(PFAM_TOP_HIT_SUFFIX, '')
            bin_id = bin_id.replace(TIGRFAM_TOP_HIT_SUFFIX, '')

            for gene_id, cur_hit in top_hits.iter_hits():
                scaffold_id = gene_id.rpartition('_')[0]
                markers[bin_id][cur_hit.hmm_id].append(scaffold_id)

        return markers

    def evaluate(self, gene_table):
        """Evaluate completeness and contamination of genome using best domain-level marker sets."""

        bac_comp, bac_cont = Markers.estimate_comp_cont(
            gene_table,
            self.bac_ms
        )

        ar_comp, ar_cont = Markers.estimate_comp_cont(
            gene_table,
            self.ar_ms
        )

        if bac_comp + bac_cont > ar_comp + ar_cont:
            return 'Bacteria', bac_comp, bac_cont

        return 'Archaea', ar_comp, ar_cont

    def create_gene_table(self, contig_ids):
        """Create gene tables for a set of contigs."""

        gene_table = defaultdict(list)
        for cid in contig_ids:
            for mid in self.markers_on_contigs[cid]:
                gene_table[mid].append(cid)

        return gene_table

    def bin_quality(self, bin):
        """Estimate quality of bin."""

        gene_table = self.create_gene_table(bin)

        return self.evaluate(gene_table)

    def marker_gene_tables(self, profile_dir, binning_methods=None):
        """Get marker genes for bins across all binning methods."""

        binning_methods_dir = os.path.join(profile_dir, BINNING_METHOD_DIR)
        gene_tables = {}
        for bm in os.listdir(binning_methods_dir):
            if binning_methods and bm not in binning_methods:
                continue  # not in set of desired binning methods

            top_hit_dir = os.path.join(binning_methods_dir, bm)
            gene_tables[bm] = Markers.read_top_hit_tables(top_hit_dir)

        self._markers_on_contigs(gene_tables)

        return gene_tables

    def _markers_on_contigs(self, gene_tables):
        """Get markers on each contig."""

        self.markers_on_contigs = defaultdict(list)
        processed_contigs = set()
        for gene_tables in gene_tables.values():
            # process all bins identified for a given binning method
            scaffolds_in_binning_method = set()
            for bin_id in gene_tables:
                for marker_id, scaffold_ids in gene_tables[bin_id].items():
                    for scaffold_id in scaffold_ids:
                        if scaffold_id in processed_contigs:
                            continue

                        self.markers_on_contigs[scaffold_id].append(marker_id)

                    scaffolds_in_binning_method.update(scaffold_ids)

            processed_contigs.update(scaffolds_in_binning_method)
