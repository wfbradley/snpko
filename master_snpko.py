#!/usr/bin/env python

import utils_snpko as utils
import check_input
import ensembl_miner
import simple_stats
import population_refiner
import find_loci
import make_knockoffs
import classifier
import sig_results
import halt_machine
import traceback
import p_values

logger = utils.logger


def master(args):
    '''
    Given a collection of SNP data and some observed features (like presence
    or absence of disease), detemine which SNPs are significant

    This is the master script that will programmatically run all of the
    constituent steps in determining which SNPs are predictive of
    which observations.  For general context for this problem, see
    README.md
    '''

    utils.check_permissions(args)

    try:
        check_input.check_and_convert_input(args)
        ensembl_miner.download_SNPs(args)
        simple_stats.stats(args)
        population_refiner.refine(args)
        find_loci.prune(args)
        if (args.machine_num == 0) or not(args.p_values):
            logger.info("####################################")
            logger.info("CAUSAL KNOCKOFFS")
            # If computing p-values and using multiple machines,
            # only compute causal knockoffs on one machine.
            # (Rest are for computing p-values.)
            make_knockoffs.make_all_knockoffs(args)
            classifier.significant_SNPs(args)
            sig_results.summarize(args)
            utils.upload_results_dir(args)
        if args.p_values:
            logger.info("####################################")
            logger.info("P-VALUE KNOCKOFFS")
            p_values.preserve_original_files(args)
            for p_trial_num in xrange(args.p_samples):
                p_values.prepare_files(args, p_trial_num)
                make_knockoffs.make_all_knockoffs(args)
                classifier.significant_SNPs(args)
                sig_results.parse_knockoff_results(args)
                p_values.upload_p_value_files(args, p_trial_num)
            p_values.extract_null_distribution(args)
        halt_machine.possibly_halt(args)
    except Exception:
        logger.warn(traceback.print_exc())

        halt_machine.possibly_halt(args)
        raise


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.initialize_logger(args)
    master(args)
