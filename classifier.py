#!/usr/bin/env python


import numpy as np
import pandas as pd
import os
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import SGDClassifier
import snpko_utils as utils
import operator

logger = utils.logger


def signficant_SNPs(args):
    '''
    Determine which SNPs are actually significant predictors of features.
    '''

    logger.info("####################################")
    logger.info("Classifier for significance.")

    max_SGD_iterations = 500
    logger.info("SGD iterations: %d" % max_SGD_iterations)

    q = 0.20
    logger.info("Target FDR: %.2f" % q)

    # Extract list of data labels (i.e., the dependent variables we're trying to
    # predict)
    df_for_field_names = pd.read_csv(
        os.path.join(args.working_dir, 'pruned_experiment.csv'))
    label_fields = [
        field for field in df_for_field_names.columns if field.startswith(args.data_prefix)]
    feature_fields = [
        field for field in df_for_field_names.columns if field.startswith('rs')]
    del df_for_field_names

    # features = np.array(df[df.columns[1:]]).astype(float)
    logger.info("Num features=%d, num labels=%d" %
                (len(feature_fields), len(label_fields)))
    # labels = np.random.randint(2,size=len(df))

    logger.info('Label fields:')
    logger.info(label_fields)

    utils.safe_mkdir(os.path.join(args.working_dir,'results'))
    for one_label_field in label_fields:
        logger.info('')
        logger.info("## %s ##" % (one_label_field))

        results_mFDR = {}
        results_FDR = {}

        for knockoff_trial in xrange(args.num_knockoff_trials):
            df = pd.read_csv(os.path.join(
                args.working_dir, 'knockoffs',
                'knockoffs_%03d.csv' % knockoff_trial))

            features = np.array(df[feature_fields]).astype(float)
            labels = np.array(df[one_label_field]).astype(float)

            # Set the parameters by cross-validation
            tuned_parameters = [{'alpha': np.power(
                10.0, np.linspace(-4, 0, num=9)),
                'l1_ratio': np.linspace(0.05, 1.00, num=20)}]
            clf = GridSearchCV(SGDClassifier(loss='log', penalty='elasticnet',
                                             max_iter=max_SGD_iterations),
                               tuned_parameters, cv=9, n_jobs=-1)
            clf.fit(features, labels)

            logger.debug("Best parameters set found on development set:")
            logger.debug(clf.best_params_)
            logger.debug("Best score: %.4f" % (clf.best_score_))
            if args.verbose:
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
            tau_mFDR = np.inf
            tau_FDR = np.inf
            for W in stat_list:
                if W <= 0:
                    continue
                if W >= tau_mFDR and W >= tau_FDR:
                    continue

                ratio_mFDR = 1.0 * (np.sum(stat_list <= -W)) / \
                    (np.sum(stat_list >= W))
                ratio_FDR = 1.0 * (1 + np.sum(stat_list <= -W)) / \
                    (np.sum(stat_list >= W))

                if ratio_mFDR < q:
                    tau_mFDR = W
                if ratio_FDR < q:
                    tau_FDR = W
            logger.debug("tau_mFDR = %.3f, tau_FDR = %.3f" %
                         (tau_mFDR, tau_FDR))

            for i, stat in enumerate(stat_list):
                SNP = feature_fields[2 * i]

                if stat >= tau_mFDR:
                    if SNP not in results_mFDR:
                        results_mFDR[SNP] = 0
                    results_mFDR[SNP] += 1
                if stat >= tau_FDR:
                    if SNP not in results_FDR:
                        results_FDR[SNP] = 0
                    results_FDR[SNP] += 1

        out_fp = open(os.path.join(args.working_dir,'results','knockoff_trials.txt'),'a')
        out_fp.write('#################\n')
        out_fp.write('Label: %s\n'%one_label_field)
        for results, fdr_type in [(results_mFDR, 'mFDR'), (results_FDR, 'FDR')]:
            out_fp.write("FDR type: %s\n" % fdr_type)
            if len(results) == 0:
                out_fp.write("No significant SNPs found.\n")
                continue
            sorted_SNPs = sorted(
                results.items(), key=operator.itemgetter(1), reverse=True)
            for (SNP, count) in sorted_SNPs:
                percentage = 100.0 * count / args.num_knockoff_trials
                out_fp.write("   %s : %d%%\n" % (SNP, percentage))
        out_fp.close()
    logger.info('Done with classifier!')

if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.initialize_logger(args)
    signficant_SNPs(args)
