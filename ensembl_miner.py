#!/usr/bin/env python

# ENSEMBL, the European data consortium, has an API for extracting SNP genotype
# data from various corpora, including 1000genomes.  We use it to
#   (*) Get sample genotypes from a large population of people
#   (*) Determine what the wild type base pair is.
#
# ENSEMBL's API is described here:
# https://rest.ensembl.org/documentation/info/variation_id

import requests
import cPickle as pickle
import pandas as pd
import os
import snpko_utils as utils
from joblib import Parallel, delayed

logger = utils.logger


def grab_individual_genotypes(SNP, cache_dir):
    '''
    Extract information about a single SNP using ENSEMBL's RESTful API.
    '''
    json_file = os.path.join(cache_dir, str(SNP))
    logger.debug("Downloading SNP %s" % SNP)

    if os.path.exists(json_file):
        # Try to grab ENSEMBL data from cache...
        decoded = pickle.load(open(json_file))
    else:
        # ...if not in cache, fetch it from server
        server = "https://rest.ensembl.org"
        ext = "/variation/human/%s?genotypes=1" % (SNP)

        r = requests.get(
            server + ext, headers={"Content-Type": "application/json"})

        if not r.ok:
            # Calling code can capture exceptions (e.g., 404) here.
            r.raise_for_status()

        decoded = r.json()
        pickle.dump(decoded, open(json_file, 'w'))

    # Determine chromosome and genomic location for SNP
    if 'mappings' not in decoded:
        raise LookupError
    for m in decoded['mappings']:
        required_fields = ['assembly_name', 'location', 'start', 'allele_string']
        missing_field=False
        for f in required_fields:
            if f not in m:
                missing_field=True
                break
        if missing_field:
            continue

        if m['assembly_name'] != 'GRCh38':
            continue
        # decoded['mappings']['location'] should look like '9:133256042-133256042'
        # If may also look like 'CHR_HG2030_PATCH:133256189-133256189';
        # we want to skip that kind.
        try:
            chromosome = int(m['location'].split(':')[0])
        except:
            continue

        loc_start = int(m['start'])
        wild_type = m['allele_string']
    try:
        # Make sure these are defined
        chromosome
        loc_start
    except:
        raise LookupError

    # Don't use ancestral type...
    # # Determine ancestral genotype
    # if 'ancestral_allele' not in decoded:
    #     raise LookupError
    # ancestral = decoded['ancestral_allele']



    # Determine genotypes of individual people
    if 'genotypes' not in decoded:
        raise LookupError
    out = {}
    for x in decoded['genotypes']:
        person = x['sample']
        genotype = x['genotype']
        if person in out:
            if out[person] != genotype:
                logger.debug("Duplicate differing samples: %s: %s: %s %s" % (
                    SNP, person, out[person], genotype))
        out[person] = genotype
    return(chromosome, loc_start, wild_type, out)


def download_SNPs(args):
    '''
    Extract list of SNPs and extract data for each.  Do it in parallel.
    '''

    logger.info("####################################")
    logger.info("Downloading SNP data from ENSEMBL")


    cache_dir = os.path.join(args.working_dir, 'ensembl_cache')
    utils.safe_mkdir(cache_dir)

    # Input is, e.g., "rs56116432"; output is a dictionary mapping
    # an individual (like '1000GENOMES:phase_3:HG00096') to
    # a genotype at the target SNP (like 'C|C')

    # SNP_plus_list=['rs12103_A','rs6667605_C']
    df_in = pd.read_csv(os.path.join(args.working_dir, 'cleaned_input.csv'))
    SNP_list = [x for x in df_in.columns if x.startswith('rs')]

    logger.info("Start download of genotype data for %d SNPs" %
                (len(SNP_list)))

    results = Parallel(n_jobs=args.num_workers)(
        delayed(grab_individual_genotypes)(SNP, cache_dir)
        for SNP in SNP_list)

    (chr_list, chr_loc_list, wild_type_list, geno_list) = zip(*results)
    df = pd.DataFrame({'SNP': SNP_list, 'chromosome': chr_list,
                       'chromosome_position': chr_loc_list, 'wild_type': wild_type_list})
    df.to_csv(os.path.join(args.working_dir, 'SNP_facts.csv'), index=False)

    # If original data did not specify wild_type, then extract reasonable
    # values from ENSEMBL.
    if not os.path.exists(os.path.join(args.working_dir,'wild_types.csv')):
        unique_wild_type=[]
        for i in xrange(len(wild_type_list)):
            unique_wild_type[i] = wild_type_list[i].split('/')[0]
        df_wild = pd.DataFrame({'SNP':SNP_list, 'wild_type': unique_wild_type})
        df_wild.to_csv(os.path.join(args.working_dir,'wild_types.csv'), index=False)

    genotypes = dict(zip(SNP_list, geno_list))
    outfile = os.path.join(args.working_dir, "ensembl.pkl")
    pickle.dump(genotypes, open(outfile, "wb"))
    logger.info("Download completed.")


if __name__ == '__main__':
    args = utils.parse_arguments()
    utils.safe_mkdir(args.working_dir)
    utils.initialize_logger(args)
    download_SNPs(args)
