#!/usr/bin/env python


import pandas as pd
import os
import numpy as np
import SNPknock.fastphase as fp
from SNPknock import knockoffHMM
from joblib import Parallel, delayed
import utils_snpko as utils
from version_snpko import __version__

logger = utils.logger


def make_knockoff(chromosome=None, grouped_by_chromosome=None, df_SNP=None,
                  df_geno_experiment=None, df_geno_ensembl=None,
                  SNP_to_wild_type=None, cache_dir=None, path_to_fp=None,
                  em_iterations=25, random_seed=123):
    #assert chromosome!=None and grouped_by_chromosome!=None and df_SNP!=None
    assert chromosome is not None
    assert grouped_by_chromosome is not None
    assert df_SNP is not None

    logger.debug("################")
    logger.debug("Chromosome %2d #" % chromosome)
    logger.debug("################")

    num_experiment_people = len(df_geno_experiment)
    num_ensembl_people = len(df_geno_ensembl)

    indices = grouped_by_chromosome.groups[chromosome]
    df_SNP_chromo = df_SNP.iloc[indices].sort_values('chromosome_position')
    SNPs_on_chromosome = df_SNP_chromo['SNP'].values

    X_experiment = np.empty((num_experiment_people, len(SNPs_on_chromosome)))
    X_ensembl = np.empty((num_ensembl_people, len(SNPs_on_chromosome)))
    for X, df in [
            (X_experiment, df_geno_ensembl),
            (X_ensembl, df_geno_ensembl)]:

        for j, SNP in enumerate(SNPs_on_chromosome):
            X[:, j] = utils.genotype_to_nonwild_type_count(
                df[SNP].values, SNP_to_wild_type[SNP] )

    out_path = '%s/chrom_%d' % (cache_dir, chromosome)

    # If all relevant files are found in cache, skip EM recomputation; otherwise,
    # redo the whole thing.
    target_file_suffix_list = [
        'alphahat.txt', 'finallikelihoods', 'origchars', 'rhat.txt', 'thetahat.txt']
    already_in_cache = True
    for suffix in target_file_suffix_list:
        target_path = os.path.join(
            cache_dir, 'chrom_%d_%s' % (chromosome, suffix))
        if not os.path.exists(target_path):
            already_in_cache = False
            break
    if already_in_cache:
        logger.debug("Found chrom %d HMM in cache" % chromosome)
    else:
        # Write array to file
        Xfp_file = '%s/X_%d.inp' % (cache_dir, chromosome)
        fp.writeX(X_ensembl, Xfp_file)

        # Run fastPhase on data (which runs EM)
        fp.runFastPhase(path_to_fp, Xfp_file, out_path,
                        K=12, numit=em_iterations)

    # Read in fastPhase results (i.e., HMM parameters) from file:
    r_file = out_path + "_rhat.txt"
    alpha_file = out_path + "_alphahat.txt"
    theta_file = out_path + "_thetahat.txt"
    # Why is X_ensembl[0, :] in the function arguments below?
    hmm = fp.loadFit(r_file, theta_file, alpha_file, X_ensembl[0, :])

    # Actually produce the knockoffs
    knockoffs = knockoffHMM(hmm["pInit"], hmm["Q"], hmm[
                            "pEmit"], seed=random_seed)
    X_knockoffs = knockoffs.sample(X_experiment)

    return(X_knockoffs, X_experiment, SNPs_on_chromosome)


def make_all_knockoffs(args):
    '''
    For each chromosome, independently:
       Sort SNPs according to position on genome.
       Train HMM parameters with EM on ENSEMBL data.
       Generate knockoffs of experimentals SNP data.

    For now, we ignore sex of persons, although that is
    available in ENSEMBL
    '''

    logger.info("####################################")
    logger.info("Fitting HMM and generating knockoffs")

    path_to_fp = os.path.join(args.fastPHASE_path, 'fastPHASE')
    if not(os.path.exists(path_to_fp)):
        logger.info("Cannot find fastPHASE at %s" % path_to_fp)
        raise Exception

    cache_dir = os.path.join(args.working_dir, 'fastphase_cache')
    utils.safe_mkdir(cache_dir)

    df_geno_ensembl = pd.read_csv(os.path.join(
        (args.working_dir), 'pruned_ensembl.csv'))

    # SNP,wild_type,chromosome,chromosome_position
    df_SNP = pd.read_csv(os.path.join((args.working_dir), 'pruned_SNP_facts.csv'))
    df_wild = pd.read_csv(os.path.join(args.working_dir,'wild_types.csv'))
    SNP_to_wild_type = dict(
        zip(df_wild['SNP'].values, df_wild['wild_type'].values))

    chromosome_list = np.sort(np.unique(df_SNP['chromosome']))
    for chromosome in chromosome_list:
        assert chromosome in np.arange(1,24)

    df_geno_experiment = pd.read_csv(os.path.join(
        (args.working_dir), 'pruned_experiment.csv'))

    # Make sure we have the same SNPs everywhere.
    assert set(df_geno_ensembl.columns) == set(df_geno_ensembl.columns)
    for SNP in df_SNP.SNP.values:
        assert SNP in df_geno_ensembl.columns

    grouped_by_chromosome = df_SNP.groupby('chromosome')
    num_experiment_people = len(df_geno_experiment)

    knockoff_SNP_list = []

    utils.safe_mkdir(os.path.join(args.working_dir,'knockoffs'))

    em_iterations = 500
    logger.info('Number of EM iterations: %d'%em_iterations)

    for knockoff_trial_count in xrange(args.num_knockoff_trials):
        random_seed = knockoff_trial_count + args.random_seed
        if ((args.num_knockoff_trials <= 20) or
                knockoff_trial_count % ((args.num_knockoff_trials) // 20) == 0):
            logger.info("Knockoff sampling %d of %d" % (
                knockoff_trial_count, args.num_knockoff_trials))

        if False:
            # Serial version; code preserved for debugging purposes
            for chromosome in chromosome_list:
                knockoff_SNP_list.append(
                    make_knockoff(
                        chromosome=chromosome,
                        grouped_by_chromosome=grouped_by_chromosome, df_SNP=df_SNP,
                        df_geno_experiment=df_geno_experiment, df_geno_ensembl=df_geno_ensembl,
                        SNP_to_wild_type=SNP_to_wild_type, cache_dir=cache_dir,
                        path_to_fp=path_to_fp, em_iterations=em_iterations, random_seed=random_seed))
        else:
            knockoff_SNP_list = Parallel(n_jobs=args.num_workers)(
                delayed(make_knockoff)(
                    chromosome=i,
                    grouped_by_chromosome=grouped_by_chromosome, df_SNP=df_SNP,
                    df_geno_experiment=df_geno_experiment, df_geno_ensembl=df_geno_ensembl,
                    SNP_to_wild_type=SNP_to_wild_type, cache_dir=cache_dir, path_to_fp=path_to_fp,
                    em_iterations=em_iterations, random_seed=random_seed)
                for i in chromosome_list)

        # Stitch results for each chromosome back together into a single dataframe
        # Knockoff results
        SNP_columns = [
            x for x in df_geno_ensembl.columns if x.startswith('rs')]
        df_knockoffs = pd.DataFrame(
            columns=SNP_columns, index=np.arange(num_experiment_people))

        # Matched experimental observations + knockoffs in one dataframe
        matched_columns = []
        data_labels = []
        for field in df_geno_experiment.columns:
            if field.startswith('rs'):
                matched_columns.append(field)
                matched_columns.append(field + '_knockoff')
            elif field.startswith(args.data_prefix):
                data_labels.append(field)
            else:
                continue
        df_matched = pd.DataFrame(columns=matched_columns+data_labels,
                                  index=np.arange(num_experiment_people))

        for (X_knockoffs, X_experiment, SNPs_on_chromosome) in knockoff_SNP_list:
            for i in xrange(num_experiment_people):
                for j, SNP in enumerate(SNPs_on_chromosome):
                    df_knockoffs[SNP].values[i] = X_knockoffs[i, j]
                    df_matched[SNP].values[i] = int(X_experiment[i, j])
                    df_matched[
                        SNP + '_knockoff'].values[i] = int(X_knockoffs[i, j])
        for data_label in data_labels:
            df_matched[data_label]=df_geno_experiment[data_label]

        # Sanity check that all fields are filled in.
        for field in df_knockoffs:
            for i in xrange(num_experiment_people):
                assert pd.notnull(df_knockoffs[field].values[i])

        df_matched.to_csv(os.path.join((args.working_dir), 'knockoffs',
                                       'knockoffs_%03d.csv' % knockoff_trial_count),
                          index=False)

    logger.info("Done making knockoffs!!!")


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.initialize_logger(args)
    make_all_knockoffs(args)
