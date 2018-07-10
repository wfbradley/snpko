#!/usr/bin/env python

import cPickle as pickle
import os
import snpko_utils as utils
import numpy as np
from version_snpko import __version__

logger = utils.logger


def refine(args):
    '''
    Not all SNPs are sequenced for all people.  Therefore, we need to choose
    a subset of SNPs and a subset of people.  There is a trade-off-- if we use
    a larger set of SNPs, we reduce the population of people, and vice versa.

    This function tries to find a reasonably good trade-off between the two.
    '''

    logger.info("####################################")
    logger.info("Balancing SNP coverage with population size.")

    ensembl_filename = os.path.join(args.working_dir, 'ensembl.pkl')
    genotypes = pickle.load(open(ensembl_filename, "rb"))

    logger.info("Initial number of SNPs: %d"%(len(genotypes)))
    SNP_list = genotypes.keys()

    person_count = {}
    for SNP in SNP_list:
        for person in genotypes[SNP]:
            if person not in person_count:
                person_count[person] = 0
            person_count[person] += 1
    logger.info("Initial number of people: %d"%(len(person_count)))

    # We are
    #   (1) going to grab a swathe of population with many SNPs (i.e., restrict
    #       the number of genomes)
    #   (2) then figure out which of our target SNP list is represented across
    #       the whole population (i.e., restrict the number of SNPs)
    #
    # This will probably all boil down to extracting the people from the 1000 Genomes
    # project, but we should hold out the hope that we might opportunistically
    # find other genomes that will work.

    def pareto_point(c):
        proposed_population = {}
        for person in person_count:
            if person_count[person] >= c:
                proposed_population[person] = True

        proposed_SNP_list = []
        for SNP in SNP_list:
            matched_all = True
            for person in proposed_population:
                if person not in genotypes[SNP]:
                    matched_all = False
                    break
            if matched_all:
                proposed_SNP_list.append(SNP)

        return(proposed_population, proposed_SNP_list)

    # Find best compromise between population size and SNP size
    best_score = -np.inf
    best_c = None
    possible_c = set(person_count.values())
    for c in possible_c:
        (proposed_population, proposed_SNP_list) = pareto_point(c)
        pop_size = len(proposed_population)
        snp_size = len(proposed_SNP_list)
        score = (1.0 * pop_size / len(person_count)
                 + (args.snp_weight) * snp_size / len(genotypes))
        if score > best_score:
            best_score = score
            best_c = c

    (proposed_population, proposed_SNP_list) = pareto_point(best_c)
    logger.info("Pareto-optimal point: score=%.3f, c=%d, |pop|=%d, |SNPs|=%d" % (
        best_score, best_c, len(proposed_population), len(proposed_SNP_list)))

    logger.info("Missing SNPs:")
    for SNP in SNP_list:
        if SNP not in proposed_SNP_list:
            logger.info("   %s" % SNP)

    f = open(os.path.join(args.working_dir, 'genotypes_ensembl.csv'), 'w')
    f.write('id')
    for SNP in proposed_SNP_list:
        f.write(',%s' % SNP)
    f.write('\n')
    for person in proposed_population:
        f.write('%s' % person)
        for SNP in proposed_SNP_list:
            f.write(',%s' % (genotypes[SNP][person]))
        f.write('\n')
    f.close()

if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    refine(args)
