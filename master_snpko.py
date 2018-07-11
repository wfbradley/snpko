#!/usr/bin/env python

import utils_snpko as utils
import check_input
import ensembl_miner
import population_refiner
import find_loci
import make_knockoffs
import classifier
from version_snpko import __version__

def master(args):
    '''
    Given a collection of SNP data and some observed features (like presence
    or absence of disease), detemine which SNPs are significant

    This is the master script that will programmatically run all of the
    constituent steps in determining which SNPs are predictive of
    which observations.  For general context for this problem, see
    README.md
    '''
    try:
        check_input.check_and_convert_input(args)
        ensembl_miner.download_SNPs(args)
        simple_stats.stats(args)
        population_refiner.refine(args)
        find_loci.prune(args)
        make_knockoffs.make_all_knockoffs(args)
        classifier.signficant_SNPs(args)
    except:
        halt_machine.possibly_halt(args)
        raise

if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.initialize_logger(args)
    master(args)
