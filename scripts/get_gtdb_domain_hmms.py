#! /usr/bin/env python3

__prog_name__ = "gtdb_gtdb_domain_hmms.py"
__prog_desc__ = "Get uniquitious, single-copy genes across GTDB species"

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2021'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'


import os
import sys
import logging
import argparse
from collections import defaultdict, namedtuple

from unitem.utils import CustomHelpFormatter, logger_setup

from biolib.common import canonical_gid
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies

from numpy import mean as np_mean, std as np_std


class GetHMMs():
    """Get uniquitious, single-copy genes across GTDB species."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        check_dependencies(['hmmpress'])

        # hard-coded paths to HMM files on ACE servers
        # (change these if using on your own system)
        self.TIGRFAM_HMMS = '/srv/db/tigrfam/15.0/TIGRFAMs_15.0_HMM/tigrfam.hmm'
        self.PFAM_HMMS = '/srv/db/pfam/33.1/Pfam-A.hmm'

    def read_hmms(self, hmm_file):
        """Read HMM models."""

        hmms = {}
        with open(hmm_file) as f:
            model = []
            for line in f:
                model.append(line)
                if line.startswith('//'):
                    # reached end of current model
                    hmms[model_acc] = ''.join(model)

                    model = []
                    model_acc = None
                elif line.startswith('ACC'):
                    model_acc = line.strip().split()[1]

        return hmms

    def read_pfam_dat(self, hmm_file):
        """Read PFAM data file."""

        dat_file = f'{hmm_file}.dat'

        hmm_dat = {}
        with open(dat_file) as f:
            dat = []
            for line in f:
                dat.append(line)
                if line.startswith('//'):
                    # reached end of current model
                    dat.append(line)
                    hmm_dat[model_acc] = ''.join(dat)

                    dat = []
                    model_acc = None
                elif line.startswith('#=GF AC'):
                    model_acc = line.strip().split()[2]

        return hmm_dat

    def get_hmms(self,
                 bac_tigr_markers,
                 bac_pfam_markers,
                 ar_tigr_markers,
                 ar_pfam_markers,
                 output_dir):
        """Pull out selected HMMs."""

        # write out selected markers
        with open(os.path.join(output_dir, f'tigr_bac.lst'), 'w') as f:
            f.write('\n'.join(bac_tigr_markers))
            f.write('\n')

        with open(os.path.join(output_dir, f'tigr_ar.lst'), 'w') as f:
            f.write('\n'.join(ar_tigr_markers))
            f.write('\n')

        with open(os.path.join(output_dir, f'pfam_bac.lst'), 'w') as f:
            f.write('\n'.join(bac_pfam_markers))
            f.write('\n')

        with open(os.path.join(output_dir, f'pfam_ar.lst'), 'w') as f:
            f.write('\n'.join(ar_pfam_markers))
            f.write('\n')

        # get TIGRFAM HMMs
        tigr_hmms = self.read_hmms(self.TIGRFAM_HMMS)
        tigr_markers = bac_tigr_markers.union(ar_tigr_markers)

        with open(os.path.join(output_dir, f'tigrfam.hmm'), 'w') as f:
            for marker in tigr_markers:
                f.write(tigr_hmms[marker])

        # get PFAM HMMs
        pfam_hmms = self.read_hmms(self.PFAM_HMMS)
        pfam_dat = self.read_pfam_dat(self.PFAM_HMMS)
        pfram_markers = bac_pfam_markers.union(ar_pfam_markers)

        with open(os.path.join(output_dir, f'pfam.hmm'), 'w') as f:
            for marker in pfram_markers:
                f.write(pfam_hmms[marker])

        with open(os.path.join(output_dir, f'pfam.hmm.dat'), 'w') as f:
            for marker in pfram_markers:
                f.write(pfam_dat[marker])

    def parse_metadata(self, metadata_file, qc_passed):
        """Parse GTDB genome metadata."""

        GenomeQuality = namedtuple('GenomeQuality', 'comp cont')

        genome_qual = {}
        gtdb_taxonomy = {}
        with open(metadata_file) as f:
            header = f.readline().strip().split('\t')

            comp_idx = header.index('checkm_completeness')
            cont_idx = header.index('checkm_contamination')
            gtdb_rep_idx = header.index('gtdb_representative')
            gtdb_taxonomy_idx = header.index('gtdb_taxonomy')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])

                comp = float(tokens[comp_idx])
                cont = float(tokens[cont_idx])
                gtdb_rep = tokens[gtdb_rep_idx].lower().startswith('t')

                if gtdb_rep and gid in qc_passed:
                    comp, cont = qc_passed[gid]
                    genome_qual[gid] = GenomeQuality(comp, cont)

                    gtdb_taxa = [t.strip()
                                 for t in tokens[gtdb_taxonomy_idx].split(';')]
                    gtdb_taxonomy[gid] = gtdb_taxa

        return genome_qual, gtdb_taxonomy

    def determine_marker_genes(self,
                               gtdb_bac_hits_file,
                               gids_to_consider,
                               min_single_copy):
        """Determining ubiquitous, single-copy marker genes."""

        # determine how frequent each marker gene is identified as
        # being single copy across the set of trusted genomes
        sc_count = defaultdict(int)
        genome_mgs = defaultdict(set)
        with open(gtdb_bac_hits_file) as f:
            gene_ids = []
            for gene_id in f.readline().strip().split('\t')[1:]:
                gene_ids.append(gene_id
                                .replace('PFAM_', '')
                                .replace('TIGR_', ''))

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                if gid not in gids_to_consider:
                    continue

                for idx, state in enumerate(tokens[1:]):
                    if state == 'Single':
                        sc_count[gene_ids[idx]] += 1

                        genome_mgs[gid].add(gene_ids[idx])

        # determine marker genes meeting desired single-copy rate
        marker_genes = {}
        for gene_id in gene_ids:
            sc_rate = 100.0 * sc_count[gene_id] / len(gids_to_consider)
            if sc_rate >= min_single_copy:
                marker_genes[gene_id] = sc_rate

        return marker_genes, genome_mgs

    def single_copy_phylum_table(self, marker_genes, gtdb_taxonomy, genome_mgs, min_single_copy, out_table):
        """Write out rate of marker genes for each phylum."""

        # get phyla in domain
        phylum_gids = defaultdict(list)
        for gid, taxa in gtdb_taxonomy.items():
            phylum = taxa[Taxonomy.PHYLUM_INDEX]
            phylum_gids[phylum].append(gid)

        # get percent of genomes with marker gene
        fout = open(out_table, 'w')
        sorted_mgs = sorted(marker_genes)
        sorted_mgs_str = '\t'.join(sorted_mgs)
        fout.write('Phylum\tNo. genomes\tAvg. single copy')
        fout.write(f'\tNo. genes >50% single copy')
        fout.write(f'\tNo. genes >={min_single_copy}% single copy')
        fout.write(f'\t{sorted_mgs_str}\n')

        for phylum in sorted(phylum_gids):
            cur_gids = phylum_gids[phylum]
            fout.write(f'{phylum}\t{len(cur_gids)}')

            mg_rates = []
            for mg in sorted_mgs:
                gid_count = 0
                for gid in cur_gids:
                    if mg in genome_mgs[gid]:
                        gid_count += 1

                mg_rates.append(100.0*gid_count/len(cur_gids))

            fout.write(f'\t{np_mean(mg_rates):.2f}')

            num_mv_rate = sum([1 for r in mg_rates if r > 50])
            fout.write(f'\t{num_mv_rate:,}')

            num_sc_rate = sum([1 for r in mg_rates if r >= min_single_copy])
            fout.write(f'\t{num_sc_rate:,}')

            for mg_rate in mg_rates:
                fout.write(f'\t{mg_rate:.2f}')
            fout.write('\n')

    def filter_mg_across_phyla(self, min_phyla_rate, taxonomy, marker_genes, genome_mgs):
        """Filter marker genes that are not predominately single copy across the majority of phyla."""

        phylum_gids = defaultdict(list)
        for gid, taxa in taxonomy.items():
            phylum = taxa[Taxonomy.PHYLUM_INDEX]
            phylum_gids[phylum].append(gid)

        filtered_mgs = set()
        for mg in marker_genes:
            phyla_sc_count = 0
            for phylum, gids in phylum_gids.items():
                gid_count = 0
                for gid in gids:
                    if mg in genome_mgs[gid]:
                        gid_count += 1

                if 100.0*gid_count/len(gids) >= min_phyla_rate:
                    phyla_sc_count += 1

            if 100.0*phyla_sc_count/len(phylum_gids) < min_phyla_rate:
                filtered_mgs.add(mg)

        return filtered_mgs

    def run(self,
            gtdb_bac_hits_file,
            gtdb_ar_hits_file,
            gtdb_bac_metadata_file,
            gtdb_ar_metadata_file,
            checkm_v2_rep_file,
            min_comp,
            max_cont,
            min_single_copy,
            min_phyla_rate,
            output_dir):
        """Get uniquitious, single-copy genes across GTDB species."""

        # get CheckM v2 quality estimates
        qc_passed = {}
        with open(checkm_v2_rep_file) as f:
            header = f.readline().rstrip().split('\t')

            comp_idx = header.index('Completeness')
            cont_idx = header.index('Contamination')

            for line in f:
                tokens = line.rstrip().split('\t')

                comp = float(tokens[comp_idx])
                cont = float(tokens[cont_idx])

                if comp >= min_comp and cont <= max_cont:
                    gid = canonical_gid(tokens[0])
                    qc_passed[gid] = (comp, cont)

        # get representative genomes meeting quality criteria
        bac_qual = None
        ar_qual = None
        bac_taxonomy = None
        ar_taxonomy = None
        for domain, domain_metadata_file in [
            ('bacterial', gtdb_bac_metadata_file),
                ('archaeal', gtdb_ar_metadata_file)]:

            self.logger.info(
                f'Identifying {domain} genomes meeting quality criteria:')

            cur_qual, cur_taxonomy = self.parse_metadata(domain_metadata_file, qc_passed)
            self.logger.info(
                f' - identified {len(cur_qual):,} genomes passing QC')

            comp = [q.comp for q in cur_qual.values()]
            self.logger.info(
                f' - completeness: {np_mean(comp):.1f} +/- {np_std(comp):.1f}')

            cont = [q.cont for q in cur_qual.values()]
            self.logger.info(
                f' - contamination: {np_mean(cont):.1f} +/- {np_std(cont):.1f}')

            if domain == 'bacterial':
                bac_qual = cur_qual
                bac_taxonomy = cur_taxonomy
            else:
                ar_qual = cur_qual
                ar_taxonomy = cur_taxonomy

        # get ubiqutious, single-copy bacterial marker genes across
        # all GTDB species representatives
        bac_mg = None
        ar_mg = None
        bac_genome_mgs = None
        ar_genome_mgs = None
        for domain, domain_qual, domain_hit_file in [
            ('bacterial', bac_qual, gtdb_bac_hits_file),
                ('archaeal', ar_qual, gtdb_ar_hits_file)]:

            self.logger.info(
                f'Identifying ubiqutious, single-copy {domain} marker genes:')
            cur_mg, cur_genome_mgs = self.determine_marker_genes(
                domain_hit_file,
                domain_qual,
                min_single_copy)
            self.logger.info(
                f' - identified {len(cur_mg):,} marker genes')

            sc_rates = [sc for sc in cur_mg.values()]
            self.logger.info(
                f' - single-copy rate: {np_mean(sc_rates):.1f} +/- {np_std(sc_rates):.1f}')

            if domain == 'bacterial':
                bac_mg = cur_mg
                bac_genome_mgs = cur_genome_mgs
            else:
                ar_mg = cur_mg
                ar_genome_mgs = cur_genome_mgs

        # remove marker genes that are not predominately single copy across all phyla
        self.logger.info('Identifying marker genes not predominately single copy across the majority of phyla:')

        filtered_bac_mgs = self.filter_mg_across_phyla(min_phyla_rate, bac_taxonomy, bac_mg, bac_genome_mgs)
        for mg in filtered_bac_mgs:
            del bac_mg[mg]

        self.logger.info(f' - removed {len(filtered_bac_mgs):,} bacterial marker genes')
        self.logger.info(f' - retained {len(bac_mg):,} bacterial marker genes')

        filtered_ar_mgs = self.filter_mg_across_phyla(min_phyla_rate, ar_taxonomy, ar_mg, ar_genome_mgs)
        for mg in filtered_ar_mgs:
            del ar_mg[mg]

        self.logger.info(f' - removed {len(filtered_ar_mgs):,} archaeal marker genes')
        self.logger.info(f' - retained {len(ar_mg):,} archaeal marker genes')

        # create table indicating single-copy rate of genes for each phylum
        bac_table = os.path.join(output_dir, 'phylum_mg_table_bac.tsv')
        self.single_copy_phylum_table(
            bac_mg,
            bac_taxonomy,
            bac_genome_mgs,
            min_single_copy,
            bac_table)

        ar_table = os.path.join(output_dir, 'phylum_mg_table_ar.tsv')
        self.single_copy_phylum_table(
            ar_mg,
            ar_taxonomy,
            ar_genome_mgs,
            min_single_copy,
            ar_table)

        # pull out bacterial TIGRFAM and Pfam HMMs
        self.logger.info(f'Pulling out TIGRfam and Pfam HMMs to: {output_dir}')
        bac_tigr_markers = set([mg for mg in bac_mg if mg.startswith('TIGR')])
        bac_pfam_markers = set([mg for mg in bac_mg if mg.startswith('PF')])
        ar_tigr_markers = set([mg for mg in ar_mg if mg.startswith('TIGR')])
        ar_pfam_markers = set([mg for mg in ar_mg if mg.startswith('PF')])
        self.get_hmms(bac_tigr_markers,
                      bac_pfam_markers,
                      ar_tigr_markers,
                      ar_pfam_markers,
                      output_dir)

        # create binary index for PFAM HMMs
        cmd = f"hmmpress {os.path.join(output_dir, 'pfam.hmm')}"
        os.system(cmd)


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('gtdb_bac_hits_file',
                        help="GTDB multi-hit marker file (e.g., gtdb_r202_bac120_multi_hits.tsv)")
    parser.add_argument('gtdb_ar_hits_file',
                        help="GTDB multi-hit marker file (e.g., gtdb_r202_ar122_multi_hits.tsv)")
    parser.add_argument('gtdb_bac_metadata_file',
                        help="GTDB metadata for bacterial genomes (e.g., bac120_metadata_r202.tsv)")
    parser.add_argument('gtdb_ar_metadata_file',
                        help="GTDB metadata for bacterial genomes (e.g., ar122_metadata_r202.tsv)")
    parser.add_argument('checkm_v2_rep_file',
                        help="file with CheckM v2 estimates for GTDB representative genomes")
    parser.add_argument('output_dir',
                        help="output director")

    parser.add_argument('--min_comp',
                        help="minimum completeness to consider genome",
                        type=float, default=95.0)
    parser.add_argument('--max_cont',
                        help="maximum contamination to consider genome",
                        type=float, default=5.0)

    parser.add_argument('--min_single_copy',
                        help="minimum single-copy rate of marker gene",
                        type=float, default=95.0)

    parser.add_argument('--min_phyla_rate',
                        help="minimum single-copy rate of marker gene within and across phyla",
                        type=float, default=75.0)

    parser.add_argument('--silent',
                        help="suppress output of logger", action='store_true')

    args = parser.parse_args()

    logger_setup(args.output_dir,
                 __prog_name__.replace('.py', '.log'),
                 __prog_name__,
                 __version__,
                 args.silent)

    try:
        p = GetHMMs()
        p.run(
            args.gtdb_bac_hits_file,
            args.gtdb_ar_hits_file,
            args.gtdb_bac_metadata_file,
            args.gtdb_ar_metadata_file,
            args.checkm_v2_rep_file,
            args.min_comp,
            args.max_cont,
            args.min_single_copy,
            args.min_phyla_rate,
            args.output_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
        raise
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
