#!/usr/bin/env python

# Run complete SNPKO trial on fake data, and confirm that it works.

import numpy as np
import pandas as pd
import os
import master_snpko
import utils_snpko as utils


logger = utils.logger

# Some parameters that we'll use throughout:
SNP_list = [
    u'rs3766606_G', u'rs6426833_A', u'rs12568930_T', u'rs11209026_G',
    u'rs2651244_G', u'rs17391694_C', u'rs6679677_C', u'rs2641348_A',
    u'rs4845604_G', u'rs670523_A', u'rs4656958_G', u'rs1801274_A']
Label_list = [u'Label_bucolic_plague', u'Label_hypochrondria']
N = 500

# We spike the data with some correlations:
correlation_indices = {0: 0, 1: 0, 2: 1, 4: 1}
SNP_likelihood_ratio_indices = {0: 3, 1: 3.7, 2: 3, 4: 3.7}
correlated = {}
SNP_LR = {}
for i in correlation_indices:
    SNP = SNP_list[i]
    label = Label_list[correlation_indices[i]]
    correlated[SNP] = label
for i in SNP_likelihood_ratio_indices:
    SNP = SNP_list[i]
    lr = SNP_likelihood_ratio_indices[i]
    SNP_LR[SNP] = lr


input_filename_base = 'test_csv'


def make_test_data(args):

    logger.info("####################################")
    logger.info("Making synthetic test data.")

    df = pd.DataFrame(columns=SNP_list + Label_list, index=range(N))

    # Set SNPs (features) to correlate or be independent, as we decide:
    # Note: Each SNP is correlated with at most one label, so we can
    # keep the probabilities independent.
    assert len(SNP_list) >= np.max(np.array(correlation_indices.keys()) + 1)
    assert len(Label_list) >= np.max(
        np.array(correlation_indices.values()) + 1)

    # Set SNPs randomly
    for SNP in SNP_list:
        wild_type = SNP.split('_')[1]  # E.g., 'A' or 'G'
        p = np.random.rand(N)  # uniform [0,1] sample
        df[SNP].values[p < 0.5] = '%s|%s' % (wild_type, wild_type)
        df[SNP].values[p >= 0.5] = 'X|X'

    # Set all labels randomly (at first), more weight to 0
    for label in Label_list:
        df[label].values[:] = (np.random.rand(N) < .25).astype(int)

    # Now overwrite some labels to instill correlation
    for label in Label_list:
        likelihood = np.ones(N)
        for SNP in correlated:
            if correlated[SNP] != label:
                continue
            likelihood[df[SNP].values == 'X|X'] *= SNP_LR[SNP]
            likelihood[df[SNP].values != 'X|X'] /= SNP_LR[SNP]
        prob = likelihood / (likelihood + 1.0)
        p = np.random.rand(N)
        df[label].values[p < prob] = 1

    logger.info("Writing test data file %s" % (args.input_file))
    df.to_csv(args.input_file, index=False)


def validate_results(args):
    logger.info("####################################")
    logger.info("Validating results on test data")

    # Check that output file discovers the expected correlations
    df_sig = pd.read_csv(os.path.join(
        args.working_dir, 'results', 'sig_results.csv'))

    true_correlations = {}
    for SNP_ in correlated:
        SNP, wild_type = SNP_.split('_')
        true_correlations[SNP] = correlated[SNP_]

    discovered_correlations = {}
    for i in xrange(len(df_sig)):
        SNP = df_sig.SNP.values[i]
        label = df_sig.label.values[i]
        discovered_correlations[SNP] = label

    # Check that all discoveries are valid
    false_positive_count = 0
    for SNP in discovered_correlations:
        if SNP not in true_correlations:
            logger.info("Invalid correlation: %s => %s" % (
                SNP, discovered_correlations[SNP]))
            false_positive_count += 1

    # Check that all true correlations are discovered
    false_negative_count = 0
    for SNP in true_correlations:
        if SNP not in discovered_correlations:
            logger.info("Missing correlation: %s => %s" % (
                SNP, true_correlations[SNP]))
            false_negative_count += 1

    precision = 1.0 - 1.0 * false_positive_count / \
        (len(discovered_correlations))
    recall = 1.0 - 1.0 * false_negative_count / (len(true_correlations))
    logger.info("Precision: %.3f" % precision)
    logger.info("Recall:    %.3f" % recall)

    if precision < 0.6 or recall < 0.6:
        logger.info("Performance too bad!")
        raise Exception
    logger.info("Test passed successfully.")


if __name__ == '__main__':
    args = utils.parse_arguments()

    # Set some of the command-line options to our liking.

    # If default working directory, redirect somewhere safer:
    if args.working_dir == 'data':
        args.working_dir = '/tmp/test_snpko'
    try:
        utils.safe_mkdir(args.working_dir)
    except Exception:
        print()
        print("Could not create working directory %s" % (args.working_dir))
        print()
        raise
    args.num_knockoff_trials = 19
    args.data_prefix = 'Label_'
    if args.input_file is None:
        # If default working directory, redirect to /tmp :
        input_file = os.path.join(args.working_dir, input_filename_base)
        args.input_file = input_file
    args.random_seed = 234

    utils.initialize_logger(args)

    # Make synthetic data with known correlations
    make_test_data(args)

    # Try to find the correlations
    master_snpko.master(args)

    # See if we got it right
    validate_results(args)
