#!/usr/bin/env python


import numpy as np
import pandas as pd
import os
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import SGDClassifier
import utils_snpko as utils
import operator
import itertools
from joblib import Parallel, delayed
import datetime
import halt_machine
import traceback


logger = utils.logger


def single_FDR(child_num, SGD_max_iterations, args, one_label_field, knockoff_trial):
    '''
    Computes both modified FDR (mFDR) and classical fdr (cFDR) for
    a single feature, trained on a single knockoff trial.
    '''
    df = pd.read_csv(os.path.join(
        args.working_dir, 'knockoffs',
        'knockoffs_%03d.csv' % knockoff_trial))

    feature_fields = [field for field in df.columns if field.startswith('rs')]

    features = np.array(df[feature_fields]).astype(float)
    labels = np.array(df[one_label_field]).astype(float)

    # Set the parameters by cross-validation
    tuned_parameters = [{'alpha': np.power(
        10.0, np.linspace(-4, 0, num=args.alpha_count)),
        'l1_ratio': np.linspace(0.05, 1.00, num=args.l1_count)}]
    # Avoid reuse of random seed:
    assert child_num < 10000
    seed = args.random_seed + child_num

    # loss='log' is logistic regression
    clf = GridSearchCV(SGDClassifier(loss='log', penalty='elasticnet',
                                     max_iter=SGD_max_iterations,
                                     random_state=seed, tol=args.tol,
                                     n_iter_no_change=args.n_iter_no_change),
                       tuned_parameters, cv=args.cv, n_jobs=1, iid=True)

    clf.fit(features, labels)

    logger.debug("Best parameters set found on development set:")
    logger.debug("GRID_PARAMS: %s" % (clf.best_params_))
    logger.debug("Best score: %.4f" % (clf.best_score_))
    # Disabling this for now; it's just too verbose.
    if False and args.verbose:
        logger.debug("Grid scores on development set:")
        means = clf.cv_results_['mean_test_score']
        stds = clf.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, clf.cv_results_['params']):
            logger.debug("%0.3f (+/-%0.03f) for %r"
                         % (mean, std * 2, params))

    # Number of features must be even, because every true SNP is paired with
    # a knockoff:
    assert len(feature_fields) % 2 == 0

    stat_list = []
    for i in xrange(0, len(feature_fields), 2):
        W = np.abs(clf.best_estimator_.coef_[0][
            i]) - np.abs(clf.best_estimator_.coef_[0][i + 1])
        stat_list.append(W)

    stat_list = np.array(stat_list)

    # Implement Equation 3.11 of Candes et al 2017
    tau_mFDR = np.inf  # Modified FDR
    tau_cFDR = np.inf  # Classical FDR
    for W in stat_list:
        if W <= 0:
            continue
        if W >= tau_mFDR and W >= tau_cFDR:
            continue

        ratio_mFDR = 1.0 * (np.sum(stat_list <= -W)) / \
            (np.sum(stat_list >= W))
        ratio_cFDR = 1.0 * (1 + np.sum(stat_list <= -W)) / \
            (np.sum(stat_list >= W))

        if ratio_mFDR < args.fdr:
            tau_mFDR = W
        if ratio_cFDR < args.fdr:
            tau_cFDR = W
    logger.debug("tau_mFDR = %.3f, tau_cFDR = %.3f" %
                 (tau_mFDR, tau_cFDR))

    SNP_list_mFDR = []
    SNP_list_cFDR = []
    for i, stat in enumerate(stat_list):
        SNP = feature_fields[2 * i]

        if stat >= tau_mFDR:
            SNP_list_mFDR.append(SNP)
        if stat >= tau_cFDR:
            SNP_list_cFDR.append(SNP)
    return(one_label_field, SNP_list_mFDR, SNP_list_cFDR)


def significant_SNPs(args):
    '''
    Determine which SNPs are actually significant predictors of features.
    '''

    logger.info("####################################")
    logger.info("Classifier for significance.")

    logger.info("SGD iterations: %d" % args.SGD_max_iterations)

    logger.info("Target FDR: %.2f" % args.fdr)

    # Extract list of data labels (i.e., the dependent variables we're trying to
    # predict)
    df_for_field_names = pd.read_csv(
        os.path.join(args.working_dir, 'pruned_experiment.csv'))
    label_fields = [
        field for field in df_for_field_names.columns if field.startswith(args.data_prefix)]
    feature_fields = [
        field for field in df_for_field_names.columns if field.startswith('rs')]
    del df_for_field_names

    logger.info("Num features=%d, num labels=%d" %
                (len(feature_fields), len(label_fields)))

    logger.info('Label fields:')
    logger.info(label_fields)

    utils.safe_mkdir(os.path.join(args.working_dir, 'results'))

    # Do the work (in parallel)
    results = (Parallel(n_jobs=args.num_workers)
               (delayed(single_FDR)(child_num, args.SGD_max_iterations, args, *x)
                for child_num, x in enumerate(itertools.product(label_fields,
                                                                xrange(args.num_knockoff_trials)))))

    summarized = {}
    for (one_label_field, SNP_list_mFDR, SNP_list_cFDR) in results:
        if one_label_field not in summarized:
            summarized[one_label_field] = {'mFDR': {}, 'cFDR': {}}
        for l, fdr in [(SNP_list_mFDR, 'mFDR'), (SNP_list_cFDR, 'cFDR')]:
            for SNP in l:
                if SNP not in summarized[one_label_field][fdr]:
                    summarized[one_label_field][fdr][SNP] = 0
                summarized[one_label_field][fdr][SNP] += 1

    out_fp = open(os.path.join(args.results_dir, 'knockoff_trials.txt'), 'w')
    out_fp.write(
        'Using the HMM knockoff framework, and applying the method %d times\n'
        'with independent knockoff samples, determine which SNPs are significant\n'
        'predictors of which data labels (i.e., dependent variables).\n\n'
        'We examine both a classical FDR (cFDR) and a modified FDR (mFDR),\n'
        'per Candes 2017, Equations 3.10 and 3.11.\n\n' % (
            args.num_knockoff_trials))
    out_fp.write('Target FDR: %.1f%%\n\n' % (100.0 * (args.fdr)))
    out_fp.write(str(datetime.datetime.now()))
    out_fp.write('\n')
    for one_label_field in summarized:
        out_fp.write('Label: %s\n' % one_label_field)
        for fdr in ['mFDR', 'cFDR']:
            out_fp.write('Type of FDR: %s\n' % fdr)
            if len(summarized[one_label_field][fdr]) == 0:
                out_fp.write('  No significant SNPs.\n')
            else:
                sorted_SNPs = sorted(
                    summarized[one_label_field][fdr].items(),
                    key=operator.itemgetter(1), reverse=True)
                for (SNP, count) in sorted_SNPs:
                    percentage = 100.0 * count / args.num_knockoff_trials
                    out_fp.write("   %s : %d%%\n" %
                                 (SNP, np.round(percentage)))
    out_fp.close()

    logger.info('Done with classifier!')


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.initialize_logger(args)

    # Because this step takes so long, we allow halting here if
    # called from command line.
    try:
        significant_SNPs(args)
        halt_machine.possibly_halt(args)
    except Exception:
        logger.warn(traceback.format_exc())

        halt_machine.possibly_halt(args)
        raise
