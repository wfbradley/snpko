# Create directory only if needed.  (Makes nested directories if needed,
#   e.g. "new1/new2/new3/".)
import numpy as np
import pandas as pd
import os
import sys
import logging
import argparse
import sklearn
import scipy
import google.cloud.storage
import google.auth
import boto3
from version_snpko import __version__


def safe_mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


logger = logging.getLogger('SNP Knockoff logger')

logger_initialized = False


def check_permissions(args):
    if os.getuid() == 0:
        logger.info('Running as root.')
        return
    # If here, we are not root.
    if args.halt:
        logger.info('Invoking "--halt" flag without "sudo".')
        logger.info('  Try: sudo ./master.py --input my_input_file.csv --halt')
        raise Exception
    return


def initialize_logger(args):
    global logger_initialized

    if logger_initialized:
        return
    logger_initialized = True

    safe_mkdir(args.working_dir)
    log_file = os.path.join(args.working_dir, 'run.log')

    if not (args.keep_old_logs):
        try:
            os.remove(log_file)
        except Exception:
            pass
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s: %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("##############  Beginning logging  ##############")
    logger.info("SNPKO version %s" % __version__)
    logger.info("Command line: %s" % (' '.join(sys.argv)))
    logger.info("Runtime parameters:")
    for f in sorted(args.__dict__.keys()):
        if f.startswith('__'):
            continue
        logger.info("   %s  :  %s" % (f, args.__dict__[f]))
    logger.info("Library versions:")
    for f in [sklearn, np, pd, scipy]:
        logger.info("   %s  :  %s" % (f.__name__, f.__version__))


def upload_file_to_gcloud(bucket_name=None, source_name=None, destination_name=None):
    credentials, project = google.auth.default()
    storage_client = google.cloud.storage.Client(credentials=credentials)
    bucket = storage_client.get_bucket(bucket_name)
    blob = bucket.blob(destination_name)
    try:
        blob.upload_from_filename(source_name)
    except Exception:
        logger.info(
            "Gcloud upload failure; check that instance has write access to 'Google storage'.")
        raise


def download_file_from_gcloud(bucket_name=None, source_name=None, destination_name=None):
    credentials, project = google.auth.default()
    storage_client = google.cloud.storage.Client(credentials=credentials)
    bucket = storage_client.get_bucket(bucket_name)
    blob = bucket.blob(source_name)
    blob.download_to_filename(destination_name)


def list_files_in_gcloud(bucket_name=None, prefix=None, delimiter=None):
    storage_client = google.cloud.storage.Client()
    bucket = storage_client.get_bucket(bucket_name)

    blobs = bucket.list_blobs(prefix=prefix, delimiter=delimiter)
    list_of_names = [blob.name for blob in blobs]
    return(list_of_names)


def download_prefix_from_gcloud(bucket_name=None, prefix=None, destination_dir=None):
    safe_mkdir(destination_dir)
    gcloud_file_list = list_files_in_gcloud(
        bucket_name=bucket_name,
        prefix=prefix)
    for f in gcloud_file_list:
        destination_name = os.path.join(destination_dir, os.path.basename(f))
        download_file_from_gcloud(bucket_name=bucket_name,
                                  source_name=f,
                                  destination_name=destination_name)


def upload_file_to_aws(bucket_name=None, source_name=None, destination_name=None,
                       SSE=True):
    s3 = boto3.client('s3')
    if SSE:
        s3.upload_file(source_name, bucket_name, destination_name,
                       ExtraArgs={'ServerSideEncryption': "AES256"})
    else:
        s3.upload_file(source_name, bucket_name, destination_name)


def download_file_from_aws(bucket_name=None, source_name=None, destination_name=None,
                           SSE=True):
    s3 = boto3.client('s3')
    s3.download_file(bucket_name, source_name, destination_name)


def upload_results_dir(args):
    # Possibly upload to cloud
    result_files = [f for f in os.listdir(args.results_dir) if os.path.isfile(
        os.path.join(args.results_dir, f))]

    for f in result_files:
        source_name = os.path.join(args.results_dir, f)
        destination_name = os.path.join('causal_%d/%s' % (
            args.original_random_seed, f))
        if args.upload_gcloud:
            upload_file_to_gcloud(bucket_name=args.bucket_name,
                                  source_name=source_name,
                                  destination_name=destination_name)
        if args.upload_aws:
            upload_file_to_aws(bucket_name=args.bucket_name,
                               source_name=source_name,
                               destination_name=destination_name)


def parse_arguments():

    parser = argparse.ArgumentParser(description='SNP Knockoffs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_file', type=str,
                        help='CSV file with SNPs and dependent variables.')
    parser.add_argument('--working_dir', type=str, default='data',
                        help='Directory for all working files and final output.')
    parser.add_argument('--results_dir', type=str, default=None,
                        help='Directory for final output (default = "results/" as a '
                        'subdirectory of the working_dir.')
    parser.add_argument('--skip_rows', type=str, default=None,
                        help='Skip rows of data file, 0-up indexing.')
    parser.add_argument('--na_threshold', type=float, default=0.5,
                        help='If fraction of N/A entries in a row or column is > this threshold, we '
                        'discard it.  "--skip_rows" is applied first.  First drop columns then rows.')
    parser.add_argument('--never_na', action='store_true', default=False,
                        help='By default, convert remaining N/As to double mutants; '
                        'if set, N/As raise exception.')
    parser.add_argument('--keep_old_logs', action='store_true', default=False,
                        help='By default, old log file is deleted; this preserves it.')
    parser.add_argument('--data_prefix', type=str, default='Imaging',
                        help='All dependent variables should be labelled with a common prefix (not '
                        'beginning with "rs" or "r"), so we can automatically detect them.')
    parser.add_argument('--num_workers', type=int, default=-1,
                        help='Number of parallel threads to use.  (Default = number of cores.)')
    parser.add_argument('--snp_weight', type=float, default=2.0,
                        help='Weight for Pareto-optimal tradeoff between population and SNP count.')
    parser.add_argument('--fastPHASE_path', type=str, default='.',
                        help='Path to "fastPHASE executable"')
    parser.add_argument('--random_seed', type=int, default=123,
                        help='Random seed for (reproducible) PRNGs')
    parser.add_argument('--num_knockoff_trials', type=int, default=100,
                        help='Because the knockoff process draws random samples, it can be'
                        ' helpful to repeat it multiple times.')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='Enable verbose logging (debug level)')
    parser.add_argument('--locus_threshold', type=float, default=0.5,
                        help='Correlation threshold for declaring two SNPs to be in the same locus.')
    parser.add_argument('--fdr', type=float, default=0.1,
                        help='Target false discover rate (FDR).')
    parser.add_argument('--obs_freq', type=float, default=0.5,
                        help='Only trust SNPs that show up in >obs_freq of the knockoff trials.')
    parser.add_argument('--p_samples', type=int, default=100,
                        help='Number of null-hypothesis samples to generate for estimating p-values.')
    parser.add_argument('--machine_num', type=int, default=0,
                        help='This is the index for this machine (in case using multiple machines.')
    parser.add_argument('--cv', type=int, default=9,
                        help='cv-fold cross-validation for the classifier.')
    parser.add_argument('--alpha_count', type=int, default=6,
                        help='Number of alpha values to try in classifier grid search.')
    parser.add_argument('--l1_count', type=int, default=20,
                        help='Number of l1-ratios to try in the classifier grid search.')
    parser.add_argument('--tol', type=str, default='1e-6',
                        help='Tolerance threshold for classifier.  Disable with "--tol=-inf".')
    parser.add_argument('--n_iter_no_change', type=int, default=25,
                        help='If classifier does not improve for this many consecutive steps, then stop.')
    parser.add_argument('--SGD_max_iterations', type=int, default=500,
                        help='Maximum iterations until SGD classifier terminates.')
    parser.add_argument('--p_values', action='store_true', default=False,
                        help='Compute p-values.  (Will make computation *much* slower.)')
    parser.add_argument('--p_thresh', type=float, default=0.05,
                        help='For convenience, summarize threshold for this p-value.')
    parser.add_argument('--bucket_name', type=str, default='snpko',
                        help='Cloud bucket name')
    parser.add_argument('--upload_gcloud', action='store_true', default=False,
                        help='Upload p-value files to Google cloud storage.')
    parser.add_argument('--upload_aws', action='store_true', default=False,
                        help='Upload p-value files to S3 on AWS.')
    parser.add_argument('--download_gcloud', action='store_true', default=False,
                        help='Download p-value files from Google cloud storage.')
    parser.add_argument('--download_aws', action='store_true', default=False,
                        help='Download p-value files from S3 on AWS.')
    parser.add_argument('--skip_p_value_accumulation', action='store_true', default=False,
                        help='Do not collect final p_values (because of machine distribution)')
    parser.add_argument('--halt', action='store_true', default=False,
                        help='Halt machine at completion.')

    args = parser.parse_args()

    # Fully qualify "~"s from path names
    if args.input_file is not None:
        args.input_file = os.path.expanduser(args.input_file)
    args.working_dir = os.path.expanduser(args.working_dir)
    args.fastPHASE_path = os.path.expanduser(args.fastPHASE_path)
    if args.results_dir is None:
        args.results_dir = os.path.join(args.working_dir, 'results')
    else:
        args.results_dir = os.path.expanduser(args.results_dir)
    args.original_results_dir = args.results_dir

    args.original_random_seed = args.random_seed
    args.random_seed += 10000000 * args.machine_num

    args.tol = float(args.tol)

    return(args)


def genotype_to_nonwild_type_count(x, wild_type, on_null=2):
    '''
    Given array with entries like 'G|T', return array
    with number of nonwild_type haplotypes per genotype.
    E.g., if wild_type='G', then
    'G|G' => 0
    'G|A' => 1
    'G|X' => 1
    'A|A' => 2
    'A|T' => 2

    If entry is null, return value specified by "on_null"
    '''
    out = np.zeros(len(x)).astype(int)
    for i in xrange(len(x)):
        if pd.isnull(x[i]):
            out[i] = on_null
            continue
        out[i] += (x[i][0] != wild_type)
        out[i] += (x[i][2] != wild_type)

    return(out)
