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

from unitem.defaults import *


class CheckM_MarkerSet():
    """A collection of marker genes organized into co-located sets."""

    def __init__(self, uid, lineage_str, num_genomes, marker_set):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.uid = uid                  # unique ID of marker set
        self.lineage_str = lineage_str  # taxonomic string associated with marker set
        self.num_genomes = num_genomes  # number of genomes used to calculate marker set
        self.marker_set = marker_set    # marker genes organized into co-located sets

    def size(self):
        """Number of marker genes and marker gene sets."""

        num_marker_genes = 0
        for m in self.marker_set:
            num_marker_genes += len(m)

        return num_marker_genes, len(self.marker_set)

    def num_markers(self):
        """Number of marker genes."""

        return self.size()[0]

    def num_sets(self):
        """Number of marker sets."""

        return len(self.marker_set)

    def get_marker_genes(self):
        """Get marker genes within marker set."""

        marker_genes = set()
        for m in self.marker_set:
            for marker in m:
                marker_genes.add(marker)

        return marker_genes

    def genome_check(self, hits, individual_markers):
        """Calculate genome completeness and contamination."""

        if individual_markers:
            present = 0
            multi_copy_count = 0
            for marker in self.get_marker_genes():
                if marker in hits:
                    present += 1
                    multi_copy_count += (len(hits[marker]) - 1)

            perc_comp = 100 * float(present) / self.num_markers()
            perc_cont = 100 * float(multi_copy_count) / self.num_markers()
        else:
            comp = 0.0
            cont = 0.0
            for ms in self.marker_set:
                present = 0
                multi_copy = 0
                for marker in ms:
                    count = len(hits.get(marker, []))
                    if count == 1:
                        present += 1
                    elif count > 1:
                        present += 1
                        multi_copy += (count - 1)

                comp += float(present) / len(ms)
                cont += float(multi_copy) / len(ms)

            perc_comp = 100 * comp / len(self.marker_set)
            perc_cont = 100 * cont / len(self.marker_set)

        return perc_comp, perc_cont


class CheckM_BinMarkerSets():
    """A collection of one or more marker sets associated with a bin."""

    # type of marker set
    TAXONOMIC_MARKER_SET = 1
    TREE_MARKER_SET = 2
    HMM_MODELS_SET = 3

    def __init__(self, bin_id, ms_type):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')
        self.marker_sets = []
        self.bin_id = bin_id
        self.ms_type = ms_type

    def most_specific_ms(self):
        """Get the most specific marker set."""

        return self.marker_sets[0]

    def read(self, line):
        """Construct bin marker set data from line."""

        tokens = line.split('\t')
        num_ms = int(tokens[1])
        for i in range(0, num_ms):
            uid = tokens[i * 4 + 2]
            lineage_str = tokens[i * 4 + 3]
            num_genomes = int(tokens[i * 4 + 4])
            marker_set = eval(tokens[i * 4 + 5])
            self.marker_sets.append(
                CheckM_MarkerSet(uid, lineage_str, num_genomes, marker_set))


class CheckM_MarkerSetParser():
    """Parse CheckM marker sets."""

    def __init__(self):
        """Initialization."""

        pass

    def parse_taxonomic_ms(self, marker_set_file):
        """Parse marker set from a taxonomic-specific marker set file."""

        with open(marker_set_file) as f:
            f.readline()  # skip header

            bin_line = f.readline()
            taxon_id = bin_line.split('\t')[0]
            bin_ms = CheckM_BinMarkerSets(
                taxon_id,
                CheckM_BinMarkerSets.TAXONOMIC_MARKER_SET)
            bin_ms.read(bin_line)

        return bin_ms


class Markers():
    """Read marker information."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        parser = CheckM_MarkerSetParser()
        bin_marker_sets = parser.parse_taxonomic_ms(CHECKM_BAC_MS)
        self.bac_ms = bin_marker_sets.most_specific_ms()

        bin_marker_sets = parser.parse_taxonomic_ms(CHECKM_AR_MS)
        self.ar_ms = bin_marker_sets.most_specific_ms()

        self.bac_markers_on_contigs = None
        self.ar_markers_on_contigs = None

    def read_table(self, marker_gene_table):
        """Read table indicating marker genes for each bin."""

        markers = defaultdict(lambda: defaultdict(list))
        if not os.path.exists(marker_gene_table):
            # did not identify any markers
            return markers

        with open(marker_gene_table) as f:
            f.readline()

            for line in f:
                line_split = line.strip().split('\t')

                bin_id = line_split[0]
                marker_id = line_split[1]
                gene_id = line_split[2]
                if '&&' in gene_id:
                    # indicates a marker gene is identified
                    # in adjacent genes and is likely an
                    # assembly or gene calling error
                    gene_ids = gene_id.split('&&')
                    gene_id = gene_ids[0]
                scaffold_id = gene_id[0:gene_id.rfind('_')]

                markers[bin_id][marker_id].append(scaffold_id)

        return markers

    def evaluate_bac(self, gene_table, individual_markers=False):
        """Evaluate completeness and contamination of genome using bacterial marker sets."""

        comp, cont = self.bac_ms.genome_check(gene_table, individual_markers)
        return comp, cont

    def evaluate_ar(self, gene_table, individual_markers=False):
        """Evaluate completeness and contamination of genome using archaeal marker sets."""

        comp, cont = self.ar_ms.genome_check(gene_table, individual_markers)
        return comp, cont

    def evaluate(self, bac_gene_table, ar_gene_table):
        """Evaluate completeness and contamination of genome using best domain-level marker sets."""

        bac_comp, bac_cont = self.evaluate_bac(bac_gene_table, True)
        ar_comp, ar_cont = self.evaluate_ar(ar_gene_table, True)

        if bac_comp + bac_cont > ar_comp + ar_cont:
            # select domain set with the larget number of identified markers
            # including those present multiple times
            bac_comp, bac_cont = self.evaluate_bac(bac_gene_table, False)
            return 'Bacteria', bac_comp, bac_cont

        ar_comp, ar_cont = self.evaluate_ar(ar_gene_table, False)
        return 'Archaea', ar_comp, ar_cont

    def bin_quality(self, bin):
        """Estimate quality of bin."""

        # create gene tables for bin
        bac_gene_table = defaultdict(list)
        ar_gene_table = defaultdict(list)
        for cid in bin:
            for marker_id in self.bac_markers_on_contigs[cid]:
                bac_gene_table[marker_id].append(cid)

            for marker_id in self.ar_markers_on_contigs[cid]:
                ar_gene_table[marker_id].append(cid)

        domain, comp, cont = self.evaluate(bac_gene_table,
                                           ar_gene_table)

        return domain, comp, cont

    def marker_gene_tables(self, profile_dir, binning_methods=None):
        """Get marker genes for bins across all binning methods."""

        markers = Markers()
        binning_methods_dir = os.path.join(profile_dir, BINNING_METHOD_DIR)
        gene_tables = {}
        for bm in os.listdir(binning_methods_dir):
            if binning_methods and bm not in binning_methods:
                continue  # not is set of desired binning methods

            bac_mg_table = os.path.join(
                binning_methods_dir, bm, CHECKM_BAC_DIR, MARKER_GENE_TABLE)
            bac_gene_tables = markers.read_table(bac_mg_table)

            ar_mg_table = os.path.join(
                binning_methods_dir, bm, CHECKM_AR_DIR, MARKER_GENE_TABLE)
            ar_gene_tables = markers.read_table(ar_mg_table)

            # fill in any bins with no marker genes
            bin_dir = os.path.join(binning_methods_dir,
                                   bm, CHECKM_BAC_DIR, 'bins')
            for bin_id in os.listdir(bin_dir):
                if not os.path.isdir(os.path.join(bin_dir, bin_id)):
                    continue

                if bin_id not in bac_gene_tables:
                    bac_gene_tables[bin_id] = {}

                if bin_id not in ar_gene_tables:
                    ar_gene_tables[bin_id] = {}

            gene_tables[bm] = (bac_gene_tables, ar_gene_tables)

        self._markers_on_contigs(gene_tables)

        return gene_tables

    def _markers_on_contigs(self, gene_tables):
        """Get markers on each contig."""

        self.bac_markers_on_contigs = defaultdict(list)
        self.ar_markers_on_contigs = defaultdict(list)
        processed_contigs = set()
        for binning_method, (bac_gene_tables, ar_gene_tables) in gene_tables.items():
            scaffolds_in_binning_method = set()
            for bin_id in bac_gene_tables:
                for marker_id, scaffold_ids in bac_gene_tables[bin_id].items():
                    for scaffold_id in scaffold_ids:
                        if scaffold_id in processed_contigs:
                            continue
                        self.bac_markers_on_contigs[scaffold_id].append(
                            marker_id)

                    scaffolds_in_binning_method.update(scaffold_ids)

            for bin_id in ar_gene_tables:
                for marker_id, scaffold_ids in ar_gene_tables[bin_id].items():
                    for scaffold_id in scaffold_ids:
                        if scaffold_id in processed_contigs:
                            continue
                        self.ar_markers_on_contigs[scaffold_id].append(
                            marker_id)

                    scaffolds_in_binning_method.update(scaffold_ids)

            processed_contigs.update(scaffolds_in_binning_method)

    def create_gene_table(self, contig_ids):
        """Create gene tables for a set of contigs."""

        bac_table = defaultdict(list)
        ar_table = defaultdict(list)
        for cid in contig_ids:
            for mid in self.bac_markers_on_contigs[cid]:
                bac_table[mid].append(cid)

            for mid in self.ar_markers_on_contigs[cid]:
                ar_table[mid].append(cid)

        return bac_table, ar_table
