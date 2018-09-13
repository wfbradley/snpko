#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import utils_snpko as utils

from scipy.stats import fisher_exact


logger = utils.logger


def stats(args):
    '''
    Compute some simple statistics on the data:
    * Univariate (uncorrected) p-value
    * (Uncorrected) likelihood ratio
    * Bonferroni corrected p-value
    '''

    df = pd.read_csv(os.path.join(args.working_dir, 'cleaned_input.csv'))
    utils.safe_mkdir(os.path.join(args.working_dir, 'results'))

    df_wild = pd.read_csv(os.path.join(args.working_dir, 'wild_types.csv'))
    SNP_to_wild_type = dict(
        zip(df_wild['SNP'].values, df_wild['wild_type'].values))

    # "features" are SNPs
    feature_list = [field for field in df.columns if field.startswith('rs')]
    # "labels" are the dependent variable (e.g., MRI observations)
    label_list = [
        field for field in df.columns if field.startswith(args.data_prefix)]
    N = len(df)

    feature_array = np.zeros((N, len(feature_list)))
    for i, feature in enumerate(feature_list):
        feature_array[:, i] = utils.genotype_to_nonwild_type_count(
            df[feature].values, SNP_to_wild_type[feature])

    label_array = np.zeros((N, len(label_list)))
    for i, label in enumerate(label_list):
        label_array[:, i] = df[label].values[:]

    # The above counts the number of non-wild-type haplotypes, so the values are
    # 0 (wild type diploid), 1, or 2.  To analyze with 2x2 contingency table, we
    # will combine 1 and 2 into a single state, so we either have "diploid wild type"
    # or not.
    feature_array[feature_array == 2] = 1

    # Uncorrected p-value
    with open(os.path.join(args.working_dir, 'results', 'uncorrected.csv'), 'w') as f:
        f.write('SNP,label,uncorrected_p_value,uncorrected_odds_ratio,'
                'bonferroni_corrected_p_value,empirical_ratio_with_imaging_feature,'
                'empirical_ratio_without_imaging_feature\n')

        contingency_table = np.zeros((2, 2))
        p_raw_array = np.zeros((len(label_list), len(feature_list)))

        for label_index, label in enumerate(label_list):
            for feature_index, feature in enumerate(feature_list):
                for label_state in [0, 1]:
                    for feature_state in [0, 1]:
                        contingency_table[feature_state, label_state] = (
                            np.sum(np.logical_and(
                                feature_array[
                                    :, feature_index] == feature_state,
                                df[label].values == label_state)))
                oddsratio, pvalue = fisher_exact(contingency_table)
                p_raw_array[label_index, feature_index] = pvalue
                bonferroni = pvalue * len(feature_list) * len(label_list)
                if bonferroni > 1.0:
                    bonferroni = 1.0

                # Unfortunately, an "imaging feature" is what we call a "label" in the
                # contingency table, not a "feature".
                empirical_ratio_with_feature = '%d/%d' % (contingency_table[1, 1],
                                                          contingency_table[1, 1] + contingency_table[0, 1])
                empirical_ratio_without_feature = '%d/%d' % (contingency_table[1, 0],
                                                             contingency_table[1, 0] + contingency_table[0, 0])
                f.write('%s,%s,%f,%f,%f,%s,%s\n' %
                        (feature, label, pvalue, oddsratio, bonferroni,
                            empirical_ratio_with_feature, empirical_ratio_without_feature))


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    stats(args)
