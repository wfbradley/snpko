#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
import os
import utils_snpko as utils
from version_snpko import __version__

logger = utils.logger


# Note: may return NaN.
def correlation(x, y):
    r = pearsonr(x, y)[0]

    if np.isnan(r):
        logger.warn('Pearson R failure; raising exception')
        raise Exception
    return(r)


def prune(args):
    '''
    Adjacent SNPs on a genome tend to be highly correlated.  Therefore,
    we can't really distinguish nearby SNPs, so we should really consider
    "genetic loci" rather than SNPs.  (This problem is sometimes called
    "linkage disequilibrium".)

    To address this problem, we look for correlated clusters of SNPs, and only
    take one representative from each cluster.  Following Sesia et al. and 
    Candes et al., we consider two SNPs correlated if the correlation 
    coefficient > 0.5.

    In principle, we have several choices for what data we use to estimate correlation.
    (1) Sesia et al reserve a fraction of their experimental data for this purpose;
        they then use a t-test to select the best representative from a correlated set
        of SNPs.  They can't reuse that fraction of data freely, but can do so in
        a crippled fashion.
    (2) We can use the underlying corpus of ENSEMBL data to estimate the correlation.
        There is typically much more data here so.
    (3) We can use all the experimental data to estimate correlation, choose the
        representative randomly, and then reuse it for classification.  Candes et al
        discuss a variant of this method, but feel that (1) is stronger.

    Because we wish to handle cases with severely restricted data, for which we cannot
    afford to reserve any data for t-tests, we will use method (3).

    Also, SNPs that are constant are of limited value, so we also remove those.
    '''

    logger.info("####################################")
    logger.info('Pruning SNP list to remove correlations within loci.')
    df_SNP = pd.read_csv(os.path.join(args.working_dir, 'SNP_facts.csv'))
    df_wild = pd.read_csv(os.path.join(args.working_dir,'wild_types.csv'))
    SNP_to_wild_type = dict(
        zip(df_wild['SNP'].values, df_wild['wild_type'].values))

    df = pd.read_csv(os.path.join(args.working_dir, 'cleaned_input.csv'))
    data_labels = [x for x in df.columns if x.startswith(args.data_prefix)]
    df_ensembl = pd.read_csv(os.path.join(
        args.working_dir, 'genotypes_ensembl.csv'))
    ensembl_SNP_list = [x for x in df_ensembl.columns if x.startswith('rs')]

    # Remove any SNPs not in ENSEMBL list, or with constant number of genotypes (e.g.,
    # if all were wild_type/wild_type, or all were wild_type/non-wild_type)
    good_SNP_vector = np.ones(len(df_SNP)).astype(bool)
    good_SNP_list = []
    bad_SNP_list = []
    for i, SNP in enumerate(df_SNP.SNP.values):
        if SNP not in ensembl_SNP_list:
            good_SNP_vector[i] = False
            continue
        reference_count = utils.genotype_to_nonwild_type_count(
            df[SNP].values, SNP_to_wild_type[SNP])
        if np.all(reference_count == reference_count[0]):
            good_SNP_vector[i] = False
            bad_SNP_list.append(SNP)
        good_SNP_list.append(SNP)

    logger.info('Dropping %d SNPs with constant genotypes:'%(len(bad_SNP_list)))
    logger.info(bad_SNP_list)

    df_SNP = df_SNP.iloc[good_SNP_vector].reset_index()
    df = df[good_SNP_list + data_labels]
    df_ensembl = df[good_SNP_list]

    grouped_by_chromosome = df_SNP.groupby('chromosome')
    logger.info('Considering %d chromosomes' %
                (len(grouped_by_chromosome.groups)))

    distinct_loci = []
    locus_count = 0
    locus_SNP_count = 0
    for chromosome in grouped_by_chromosome.groups:
        indices = grouped_by_chromosome.groups[chromosome]
        df_SNP_chromo = df_SNP.iloc[indices].sort_values('chromosome_position')
        SNPs_on_chromosome = df_SNP_chromo['SNP'].values
        SNP_run = []
        if len(SNPs_on_chromosome) == 1:
            # Only one SNP on the chromosome, so no possible correlation...
            distinct_loci.append(SNPs_on_chromosome[0])
            logger.info('Only one SNP on chromosome %d' % chromosome)
            continue
        for i in xrange(len(SNPs_on_chromosome)):
            SNP_run.append(SNPs_on_chromosome[i])
            if i < len(SNPs_on_chromosome) - 1:
                reference_count_1 = utils.genotype_to_nonwild_type_count(
                    df[SNPs_on_chromosome[i]].values,
                    SNP_to_wild_type[SNPs_on_chromosome[i]])
                reference_count_2 = utils.genotype_to_nonwild_type_count(
                    df[SNPs_on_chromosome[i + 1]].values,
                    SNP_to_wild_type[SNPs_on_chromosome[i + 1]])

                corr = correlation(reference_count_1, reference_count_2)
            else:
                corr = 0.0
            if corr < args.locus_threshold:
                SNP = SNP_run[np.random.randint(len(SNP_run))]
                distinct_loci.append(SNP)
                if len(SNP_run) > 1:
                    locus_count += 1
                    locus_SNP_count += len(SNP_run)
                SNP_run = []

    logger.info("Created %d loci from %d underlying SNPs" %
                (locus_count, locus_SNP_count))

    df[distinct_loci + data_labels].to_csv(os.path.join(args.working_dir, 'pruned_experiment.csv'),
                             index=False)

    df_ensembl[distinct_loci].to_csv(os.path.join(
        args.working_dir, 'pruned_ensembl.csv'), index=False)

    index = np.zeros(len(df_SNP)).astype(bool)
    for i in xrange(len(df_SNP)):
        if df_SNP['SNP'].values[i] in distinct_loci:
            index[i] = True
    df_SNP.iloc[index].reset_index().to_csv(os.path.join(args.working_dir,
                                                         'pruned_SNP_facts.csv'),
                                            index=False)


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.initialize_logger(args)
    prune(args)
