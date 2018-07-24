#!/usr/bin/env python

import pandas as pd
import os
import utils_snpko as utils
from version_snpko import __version__


logger = utils.logger


def summarize(args):
    '''
    Summarize results about significantly predictive SNPs.
    '''

    logger.info("####################################")
    logger.info("Summarizing final results")

    df_uncorrected = pd.read_csv(os.path.join(
        args.working_dir, 'results', 'uncorrected.csv'))
    grouped_uncorrected = df_uncorrected.groupby(['SNP', 'label'])

    # Filter uncorrected down to uncorrected p-value <0.05
    df_uncorrected.iloc[df_uncorrected['uncorrected_p_value'].values < 0.05].to_csv(
        os.path.join(args.working_dir, 'results', 'exploratory.csv'), index=False)

    # Parse Knockoff results
    label = None
    fdr_type = None
    SNP = None
    result_table = []
    with open(os.path.join(args.working_dir, 'results', 'knockoff_trials.txt')) as fp:
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

                if obs_freq > args.obs_freq:
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
    df_results.to_csv(os.path.join(args.working_dir,
                                   'results', 'sig_results.csv'), index=False)


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    summarize(args)
