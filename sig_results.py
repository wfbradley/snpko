#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import utils_snpko as utils

logger = utils.logger


def parse_knockoff_results(args, df_uncorrected=None):
    if df_uncorrected is None:
        df_uncorrected = pd.read_csv(os.path.join(
            args.results_dir, 'uncorrected.csv'))

    grouped_uncorrected = df_uncorrected.groupby(['SNP', 'label'])

    label = None
    fdr_type = None
    SNP = None
    result_table = []
    with open(os.path.join(args.results_dir, 'knockoff_trials.txt')) as fp:
        for line in fp:
            # chomp:
            line = line[:-1]
            if line.startswith('Target FDR: '):
                # The observed FDR takes precedence over args.fdr
                fdr = float(line[12:-1]) / 100.0
            elif line.startswith('Label: '):
                label = line[7:]
            elif line.startswith('Type of FDR: '):
                fdr_type = line[13:]
            elif line.startswith('   rs'):
                (SNP, obs_freq) = line[3:].split(' : ')
                obs_freq = float(obs_freq[:-1]) / 100.0

                index = grouped_uncorrected.groups[(SNP, label)][0]
                uncorrected_p_value = df_uncorrected[
                    'uncorrected_p_value'].values[index]
                uncorrected_odds_ratio = df_uncorrected[
                    'uncorrected_odds_ratio'].values[index]

                results = (label, fdr_type, SNP, obs_freq,
                           uncorrected_p_value, uncorrected_odds_ratio, fdr)
                result_table.append(results)

    # Produce simpler summary output
    (label_list, fdr_type_list, SNP_list, obs_freq_list,
        uncorrected_p_value_list, uncorrected_odds_ratio_list, fdr_list) = zip(*result_table)
    df_results = pd.DataFrame(
        {'label': label_list, 'fdr_type': fdr_type_list,
         'SNP': SNP_list, 'obs_freq': obs_freq_list,
         'uncorrected_p_value': uncorrected_p_value_list,
         'uncorrected_odds_ratio': uncorrected_odds_ratio_list,
         'fdr': fdr_list
         })

    df_results.to_csv(os.path.join(args.results_dir, 'all_results.csv'),
                      index=False)
    return(df_results, fdr)


def summarize(args):
    '''
    Summarize results about significantly predictive SNPs.
    '''

    logger.info("####################################")
    logger.info("Summarizing final results")

    df_uncorrected = pd.read_csv(os.path.join(
        args.results_dir, 'uncorrected.csv'))

    # Filter uncorrected down to uncorrected p-value <0.05
    df_uncorrected.iloc[df_uncorrected['uncorrected_p_value'].values < 0.05].to_csv(
        os.path.join(args.results_dir, 'exploratory.csv'), index=False)

    (df_results, fdr) = parse_knockoff_results(args, df_uncorrected=df_uncorrected)

    # Restrict to SNPs that occur more frequently than a target threshold
    df_sig_threshold = df_results.iloc[
        df_results.obs_freq.values > args.obs_freq]
    df_sig_threshold.to_csv(os.path.join(args.results_dir, 'sig_results.csv'),
                            index=False)

    # Alternately, extract the single most-frequently occuring SNP of each type
    grouped = df_results.groupby(['fdr_type', 'label'])
    max_index = []
    for _, df_fdr_label in grouped:
        local_index_of_biggest = np.argmax(df_fdr_label.obs_freq.values)
        max_index.append(df_fdr_label.index[local_index_of_biggest])
    df_sig_max = df_results.iloc[max_index]
    df_sig_max = df_sig_max.sort_values(by='obs_freq', ascending=False)
    df_sig_max.to_csv(os.path.join(args.results_dir, 'sig_max.csv'),
                      index=False)

    # Add the probabilities across all trials (gives something like the expected number
    # of times that a particular SNP shows up in *any* trial)
    df_expected = df_results[['SNP', 'fdr_type', 'label', 'obs_freq']]
    df_expected = df_expected.groupby(['SNP', 'fdr_type']).sum().reset_index()
    df_expected['fdr'] = fdr
    df_expected.rename(columns={'obs_freq': 'expected_obs_freq'}, inplace=True)
    df_expected = df_expected.sort_values(by='expected_obs_freq', ascending=False)
    df_expected.to_csv(os.path.join(args.results_dir, 'expected_appearance.csv'),
                       index=False)


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    summarize(args)
