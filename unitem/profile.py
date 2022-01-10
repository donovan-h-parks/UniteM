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
import shutil
import logging
from collections import defaultdict

from unitem.utils import check_dependencies
from unitem.external.prodigal import Prodigal
from unitem.external.tigrfam_search import TigrfamSearch
from unitem.external.pfam_search import PfamSearch
from unitem.external.tophits import TopHitFile
from unitem.markers import Markers
from unitem.defaults import (BINNING_METHOD_DIR,
                             TIGRFAM_TOP_HIT_SUFFIX,
                             PFAM_TOP_HIT_SUFFIX)


class Profile():
    """Profile genomes across different binning methods."""

    def __init__(self, cpus):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        check_dependencies(['prodigal', 'hmmsearch'])

        self.logger.info(f'Using Prodigal {Prodigal.version()}.')
        self.logger.info(f'Using hmmsearch {TigrfamSearch.version()}.')

        self.cpus = cpus

    def _genome_quality(self, gids, result_dir, tigr_search, pfam_search):
        """Determine genome quality based on presence/absence of marker genes."""

        pfam_bac_ms, pfam_ar_ms = pfam_search.get_marker_genes()
        tigr_bac_ms, tigr_ar_ms = tigr_search.get_marker_genes()

        bac_ms = pfam_bac_ms.union(tigr_bac_ms)
        ar_ms = pfam_ar_ms.union(tigr_ar_ms)

        fout_bac = open(os.path.join(result_dir, 'bac_hits.tsv'), 'w')
        fout_bac.write('Genome ID')
        for ms in bac_ms:
            fout_bac.write(f'\t{ms}')
        fout_bac.write('\n')

        fout_ar = open(os.path.join(result_dir, 'ar_hits.tsv'), 'w')
        fout_ar.write('Genome ID')
        for ms in ar_ms:
            fout_ar.write(f'\t{ms}')
        fout_ar.write('\n')

        genome_qual = {}
        for gid in gids:
            hmm_hits = defaultdict(list)

            # get top hits to PFAM HMMs
            pfam_top_hit_file = os.path.join(
                result_dir, f'{gid}{PFAM_TOP_HIT_SUFFIX}')
            pfam_top_hits = TopHitFile(pfam_top_hit_file)
            pfam_top_hits.read()

            for gene_id, cur_hit in pfam_top_hits.iter_hits():
                scaffold_id = gene_id.rpartition('_')[0]
                hmm_hits[cur_hit.hmm_id].append(scaffold_id)

            # get top hits to TIGRFAM HMMs
            tigr_top_hit_file = os.path.join(
                result_dir, f'{gid}{TIGRFAM_TOP_HIT_SUFFIX}')
            tigr_top_hits = TopHitFile(tigr_top_hit_file)
            tigr_top_hits.read()

            for gene_id, cur_hit in tigr_top_hits.iter_hits():
                scaffold_id = gene_id.rpartition('_')[0]
                hmm_hits[cur_hit.hmm_id].append(scaffold_id)

            bac_comp, bac_cont = Markers.estimate_comp_cont(
                hmm_hits,
                bac_ms
            )

            ar_comp, ar_cont = Markers.estimate_comp_cont(
                hmm_hits,
                ar_ms
            )

            if bac_comp + bac_cont > ar_comp + ar_cont:
                genome_qual[gid] = ('Bacteria', bac_comp, bac_cont)

                fout_bac.write(f'{gid}')
                for ms in bac_ms:
                    fout_bac.write(f'\t{len(hmm_hits[ms])}')
                fout_bac.write('\n')
            else:
                genome_qual[gid] = ('Archaea', bac_comp, bac_cont)

                fout_ar.write(f'{gid}')
                for ms in ar_ms:
                    fout_ar.write(f'\t{len(hmm_hits[ms])}')
                fout_ar.write('\n')

        fout_bac.close()
        fout_ar.close()

        return genome_qual

    def _report_genome_quality(self, genome_quality, output_dir):
        """Summarize quality of genomes."""

        # create table for each binning method
        for bm in genome_quality:
            table = os.path.join(output_dir, bm + '_quality.tsv')
            fout = open(table, 'w')
            fout.write(
                'Genome ID\tMarker Set Domain\tCompleteness (%)\tContamination (%)\tQuality\n')
            for gid in genome_quality[bm]:
                domain, comp, cont = genome_quality[bm][gid]
                fout.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\n'.format(
                           gid,
                           domain,
                           comp,
                           cont,
                           comp-5*cont))
            fout.close()

        # report global results file
        summary = {}
        for bm in genome_quality:
            comp_cont_count = defaultdict(lambda: defaultdict(int))
            quality_count = defaultdict(int)
            total_comp = 0
            total_cont = 0
            total_q = 0
            for _bid, (domain, comp, cont) in genome_quality[bm].items():
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

                summary[bm] = (comp_cont_count, quality_count,
                               total_comp, total_cont, total_q)

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

    def run(self, bin_dirs, marker_dir, keep_intermediate, output_dir):
        """Profile genomes in each bin directory.

        Parameters
        ----------
        bin_dirs : list of str
            Directories containing bins from different binning methods.
        marker_dir : str
            Directory containing Pfam and TIGRfam marker gene data.
        output_dir : str
            Output directory.
        """

        self.logger.info(f'Profiling genomes in {len(bin_dirs)} directories.')

        num_processed = 0
        genome_quality = defaultdict(lambda: dict)
        prodigal = Prodigal(self.cpus)
        for method_id, (bin_dir, bin_ext) in bin_dirs.items():
            num_processed += 1
            self.logger.info('Profiling {} ({} of {}):'.format(
                method_id,
                num_processed,
                len(bin_dirs)))

            # path to genomic FASTA files
            genomes = {}
            for gf in os.listdir(bin_dir):
                if not gf.endswith(bin_ext):
                    continue

                gid = gf.replace(f'.{bin_ext}', '')
                genomes[gid] = os.path.join(bin_dir, gf)

            # identify genes using Prodigal
            self.logger.info(
                f' - identifying genes in {len(genomes):,} genomes using Prodigal')

            method_dir = os.path.join(output_dir,
                                      BINNING_METHOD_DIR,
                                      method_id)
            tmp_marker_gene_dir = os.path.join(method_dir, 'intermediate')

            prodigal_files = prodigal.run(genomes, tmp_marker_gene_dir)

            # annotated genes against TIGRFAM database
            self.logger.info(' - identifying TIGRFAM markers using HMMER')
            aa_gene_files = [prodigal_files[gid]['aa_gene_path']
                             for gid in prodigal_files.keys()]
            tigr_search = TigrfamSearch(marker_dir, self.cpus)
            tigr_search.run(aa_gene_files, tmp_marker_gene_dir, method_dir)

            # annotate genes against Pfam database
            self.logger.info(' - identifying Pfam markers using HMMER')
            pfam_search = PfamSearch(marker_dir, self.cpus)
            pfam_search.run(aa_gene_files, tmp_marker_gene_dir,  method_dir)

            genome_quality[method_id] = self._genome_quality(
                genomes,
                method_dir,
                tigr_search,
                pfam_search)

            if not keep_intermediate:
                shutil.rmtree(tmp_marker_gene_dir)

        self._report_genome_quality(genome_quality, output_dir)
