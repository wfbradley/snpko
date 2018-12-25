#!/usr/bin/env python

import os
import shutil
import pandas as pd
import numpy as np
import utils_snpko as utils

logger = utils.logger


# Running "master_snpko.py" will run a series of knockoff trials aimed to
# generate lists of SNPs with a target FDR.  For a particular SNP and FDR
# threshold, then, we might discover that it appears in, say, 70% of all
# knock-off runs.

# When master_snpko.py is run, it runs multiple knockout trials.  We are
# then in the position of knowing that, e.g., SNP rs12345 appears in 85% of
# the FDR lists for label "Imaging: Fistula", with a 10% mFDR for each list.
# So, should we consider 85% to be significant?

# Let's refer to this probability (the "85%" above) as Q.  So, should we
# consider the fact that Q=85% to be significant?

# To address this question, we do the following.
#     We examine the experimental data to determine
#         The SNPs of interest,
#         The number of people in the trials (say, X subjects), and
#         The labels for those people
#     We pull background data from ENSEMBL, with SNPs from Y>>X people.
#     for trial=1, ..., T=100:
#         Randomly partition ENSEMBLE into X subjects vs Y-X subjects.
#         Replace the experimental SNPs with the SNPs from the X
#             ENSEMBL subjects
#         Restrict ENSEMBL data to the Y-X subjects.
#         Proceed as usual with SNPKO.
#         For each label (e.g., "Imaging: Fistula"), each SNP shows up
#             in some fraction of the FDR lists.  Record the maximum,
#             across all SNPS, of Q.  I.e., for a given label,
#                 max_{SNP} Q_{SNP, label}

# Once all T=100 trials are done, we now have a distribution for the maximum
# Q value for the null hypothesis (that there is no correlation between the
# SNPs and the binary output).  So, for example, if we're interested in a
# p=0.05 significance level, we consider all <SNP,label> pairs that are
# greater than or equal to the 95th largest Q.

def preserve_original_files(args):
    '''
    Save original files in sub-directory; we'll be rewriting them
    in "prepare_files()".  (This function should be called once.)
    '''
    orig_dir = os.path.join(args.working_dir, 'original')
    utils.safe_mkdir(orig_dir)
    file_prefix_list = ['ensembl', 'experiment']
    for file_prefix in file_prefix_list:
        filebase = 'pruned_%s.csv' % file_prefix
        filename_old = os.path.join(args.working_dir, filebase)
        filename_new = os.path.join(orig_dir, filebase)
        os.rename(filename_old, filename_new)
    os.rename(os.path.join(args.results_dir, 'uncorrected.csv'),
              os.path.join(orig_dir, 'uncorrected.csv'))


def prepare_files(args, p_trial_num):
    '''
    Slice and dice the files saved in "preserve_original_files()" and extract
    a random subset for generating one p-value sample.  (This function
    will be called many times.)
    '''

    logger.info('#############################')
    logger.info('#############################')
    logger.info('Preparing next set of files')

    # Remove any old cruft
    file_prefix_list = ['ensembl', 'experiment']
    for file_prefix in file_prefix_list:
        filebase = 'pruned_%s.csv' % file_prefix
        filename_old = os.path.join(args.working_dir, filebase)
        try:
            os.unlink(filename_old)
        except Exception:
            pass
    if False:
        dir_list = ['fastphase_cache', 'knockoffs']
        for d in dir_list:
            try:
                shutil.rmtree(os.path.join(args.working_dir, d),
                              ignore_errors=True)
            except Exception:
                pass

    # Change random seed
    args.random_seed += 10000
    logger.info('Random seed is now %d' % args.random_seed)

    # Create and point to bespoke output directory
    args.results_dir = os.path.join(
        args.working_dir, 'results_%03d' % (p_trial_num))
    utils.safe_mkdir(args.results_dir)

    # Read in real (original) experiment and ENSEMBL data
    orig_dir = os.path.join(args.working_dir, 'original')
    df_experiment = pd.read_csv(os.path.join(
        orig_dir, 'pruned_experiment.csv'))
    df_ensembl = pd.read_csv(os.path.join(orig_dir, 'pruned_ensembl.csv'))

    # Partition off random subset of data to replace experiment data
    num_subjects = len(df_experiment)
    num_ensembl_total = len(df_ensembl)
    assert num_ensembl_total >= num_subjects
    fake_subject_index = np.random.permutation(
        np.arange(num_ensembl_total).astype(int))[:num_subjects]
    for i, j in enumerate(fake_subject_index):
        for col in df_ensembl:
            df_experiment[col].values[i] = df_ensembl[col].values[j]
    remaining_ensembl_index = np.ones(num_ensembl_total).astype(bool)
    remaining_ensembl_index[fake_subject_index] = False
    df_ensembl = df_ensembl.iloc[remaining_ensembl_index].reset_index()
    assert len(df_ensembl) == num_ensembl_total - num_subjects

    # Write new, doctored version of data.
    shutil.copyfile(os.path.join(orig_dir, 'uncorrected.csv'),
                    os.path.join(args.results_dir, 'uncorrected.csv'))
    df_experiment.to_csv(os.path.join(
        args.working_dir, 'pruned_experiment.csv'), index=False)
    df_ensembl.to_csv(os.path.join(
        args.working_dir, 'pruned_ensembl.csv'), index=False)


def upload_p_value_files(args, p_trial_num):
    source_filename = os.path.join(args.results_dir, 'all_results.csv')
    destination_name = os.path.join('all_results_%d_%d_%03d.csv' % (
        args.original_random_seed, args.machine_num, p_trial_num))
    if args.upload_gcloud:
        utils.upload_file_to_gcloud(bucket_name=args.bucket_name,
                                    source_filename=source_filename,
                                    destination_name=destination_name)
    if args.upload_aws:
        utils.upload_file_to_aws(bucket_name=args.bucket_name,
                                 source_filename=source_filename,
                                 destination_name=destination_name)


def extract_null_distribution(args):
    p_dir = os.path.join(args.working_dir, 'p_values')
    utils.safe_mkdir(p_dir)

    if args.download_gcloud:
        logger.info('Gcloud storage download not implemented.')
        pass
    elif args.download_aws:
        logger.info('AWS S3 download not implemented.')
        pass
    else:
        # Just use local files
        for p_trial_num in xrange(args.p_samples):
            try:
                src_file = os.path.join(args.working_dir,
                                        'results_%03d' % (p_trial_num),
                                        'all_results.csv')
                dst_file = os.path.join(p_dir,
                                        'all_results_%d_%d_%03d.csv' % (
                                            args.original_random_seed,
                                            args.machine_num,
                                            p_trial_num))
                shutil.copyfile(src_file, dst_file)
            except Exception:
                logger.info("Failed to copy %s" % src_file)

    # Figure out how many null-hypothesis samples we *actually* have
    null_hypo_files = [f for f in os.listdir(p_dir)
                       if os.path.isfile(os.path.join(p_dir, f)) and f.startswith('all_results')]
    logger.info('Found %d files for the null hypothesis' %
                (len(null_hypo_files)))
    if len(null_hypo_files) == 0:
        return

    # Peek at an arbitrary file to determine the label types (e.g., "Imaging:
    # Fistula")
    df = pd.read_csv(os.path.join(p_dir, null_hypo_files[0]))
    label_types = np.unique(df.label.values)
    max_obs_freq = []

    # Extract max obs_freq from each file for each label type
    for f in null_hypo_files:
        df = pd.read_csv(os.path.join(p_dir, f))
        for label in label_types:
            M = np.max(df.iloc[df.label.values == label].obs_freq.values)
            max_obs_freq.append((label, M))

    (label_list, max_list) = zip(*max_obs_freq)
    df_null_hypo = pd.DataFrame(
        {'label': label_list, 'max_obs_freq': max_list})
    df_null_hypo.to_csv(os.path.join(
        args.original_results_dir, 'null_hypothesis.csv'), index=False)

    # For fun, report the p=0.05 values.
    logger.info("p=%.2f threshold for each label" % (args.p_thresh))
    for label in sorted(label_list):
        samples = np.sort(df_null_hypo.iloc[
                          df_null_hypo.label.values == label].max_obs_freq.values)
        v = samples[int(np.ceil((1.0 - args.p_thresh) * len(samples)))]
        logger.info("   %s : %.1f%%" % (label, v))


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    preserve_original_files(args)
