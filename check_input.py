#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import utils_snpko as utils
import re

logger = utils.logger


def check_and_convert_input(args):
    '''
    Check that the input file exists and is well-formed.  On failure,
    raise exception.  On failure, returns data frame with SNP data.
    '''

    if args.input_file is None:
        logger.info('Must specify "--input_file" argument.')
        raise Exception
    assert not(args.data_prefix.startswith('rs'))
    assert not(args.data_prefix == 'r')
    assert (args.na_threshold) <= 1.0

    logger.info("####################################")
    logger.info('Validating original input file %s, sanitizing, and rewriting.' %
                args.input_file)

    sanitized_outfile = os.path.join(args.working_dir, 'cleaned_input.csv')

    if not os.path.exists(args.input_file):
        logger.error("Failure to find input file %s" % (args.input_file))
        raise OSError
    if args.skip_rows is not None:
        args.skip_rows = [int(x) for x in args.skip_rows.split(',')]
    if (args.input_file).lower().endswith('xlsx'):
        df = pd.read_excel(args.input_file, skiprows=args.skip_rows)
    else:
        df = pd.read_csv(args.input_file, skiprows=args.skip_rows)

    # SNPs are columns of the form "rs12345*" where the "*" can be any
    # other stuff (that doesn't begin with a digit, of course).
    #
    # In particular, it may be "rs12345_G", where "G" is the wild type; if
    # so, we record that information, then normalize the column name.
    #
    # If value looks like e.g. "GT", we normalize to "G|T" to be consistent
    # with ENSEMBL format.
    #
    # Input may have base "X" indicating "non-ancestral"; if so, we
    # pass that along.

    # Truncate column names if necessary to "rs12345".
    p = re.compile(r'^rs\d+')
    SNP_columns = {}
    for c in df.columns:
        x = p.match(c)
        if x != None:
            SNP_columns[c] = c[:x.end()]

    # If "rs12345_G" format, then extract wild-type info.
    p2 = re.compile(r'^rs\d+_[ACGT]$')
    all_rs_fields_match = True
    for c in df.columns:
        if p.match(c):
            if not(p2.match(c)):
                all_rs_fields_match = False
                break
    if all_rs_fields_match:
        SNP_list_1 = []
        wild_type_list_1 = []
        logger.info('All SNP fields of form "rs12345_G"; extracting wild-type.')
        for c in df.columns:
            if p2.match(c):
                (SNP, wild_type) = c.split('_')
                SNP_list_1.append(SNP)
                wild_type_list_1.append(wild_type)
        df_wild = pd.DataFrame(
            {'SNP': SNP_list_1, 'wild_type': wild_type_list_1})
        df_wild.to_csv(os.path.join(args.working_dir,
                                    'wild_types.csv'), index=False)

    df.rename(columns=SNP_columns, inplace=True)
    SNP_columns = SNP_columns.values()

    # Convert "GT" to "G|T" (to be consistent with ENSEMBL formatting)
    p3 = re.compile(r'^[ACGTX]{2}$')
    p4 = re.compile(r'^[ACGTX]|[ACGTX]$')
    for SNP in SNP_columns:
        for i in xrange(len(df)):
            if pd.notnull(df[SNP].values[i]):
                geno = df[SNP].values[i].upper()

                if (not p3.match(geno)) and (not p4.match(geno)):
                    logger.info(
                        '[col=%s, row=%d] Expect genotype like "GT" or "G|T" but got %s' % (
                            SNP, i, geno))
                    raise Exception
                if len(geno) == 2:
                    geno = geno[0] + '|' + geno[1]
                df[SNP].values[i] = geno

    data_columns = [f for f in df.columns if f.startswith(args.data_prefix)]

    if len(SNP_columns) == 0:
        logger.error(
            'SNP columns must be of form "rs12345"; no such columns found.')
        raise Exception
    if len(data_columns) == 0:
        logger.error(
            'Some data columns must begin with "%s" (set with --data_prefix); '
            'none found.' % args.data_prefix)
        raise Exception
    relevant_columns = SNP_columns + data_columns
    df = df[relevant_columns]

    # Drop columns with too many N/As
    threshold = int((args.na_threshold) * len(df))
    drop_col = []
    for f in relevant_columns:
        if df[f].isnull().sum() > threshold:
            drop_col.append(f)
    df.drop(columns=drop_col, inplace=True)
    logger.info('Dropping %d columns because of N/As:' % (len(drop_col)))
    logger.info(drop_col)

    # Drop rows with too many N/As
    threshold = int((args.na_threshold) * len(relevant_columns))
    drop_row = []
    for i in xrange(len(df)):
        if df.iloc[i].isnull().sum() > threshold:
            drop_row.append(i)
    df.drop(labels=drop_row, inplace=True)
    logger.info('Dropping %d rows because of N/As:' % (len(drop_row)))
    logger.info(drop_row)

    # Confirm that data columns have binary (0/1) data
    for f in data_columns:
        if np.sum(df[f].values == 0) + np.sum(df[f].values == 1) != len(df):
            logger.error('Field %s has non-binary data!' % f)
            raise Exception

    if len(df.columns) == 0:
        logger.error('After dropping N/A columns, no SNPs are left!.')
        raise Exception
    if len(df) == 0:
        logger.error('After dropping N/A rows, no patients are left!.')
        raise Exception

    num_na = df.isnull().sum(axis=1).sum()
    logger.info("Detected %d N/A entries in input data" % num_na)
    if num_na > 0 and args.never_na:
        logger.error('Failing because of %d N/A values!' % num_na)
        raise Exception

    df.to_csv(sanitized_outfile, index=False)
    logger.info('Detected %d SNPs and %d patients.' %
                (len(df.columns), len(df)))
    logger.info('Sanitized data rewritten as %s' % sanitized_outfile)
    return


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    check_and_convert_input(args)
