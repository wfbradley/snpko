# A Pipeline for Statistical Knockouts with SNPs

##### Table of Contents  
[The Problem](#the-problem)  
[A Little History](#a-little-history)  
[Installation](#installation)  
[Quick Start Guide](#quick-start-guide)  
[More Detailed Guide](#more-detailed-guide)  
[Running Time and Space](#running-time-and-space)  
[Results](#results)  
[Author](#author)  
[License](#license)  
[Acknowledgements](#acknowledgements)  


## The Problem

The specific problem we addressed involved a population of 54 pediatric
patients suffering from Crohn's disease. For each patient, we measured 168
SNPs and observed 8 different radiological imaging features (such as "Lumen
narrowing" or "Small bowel disease").  We were interested to determine which
SNPs, if any, were significant predictors of particular imaging features.

This software provides tools for answering this type of question more
generally: given any set of SNPs and labels (i.e., dependent variables), these
tools will:
* Download relevant information about the SNPs from the ENSEMBL online repository, 
* Train an HMM on the data,
* Construct statistical knockouts,
* Train a classifier on the data, and
* Produce a set of SNPs with a controlled false discovery rate.

Some of these steps are computationally intensive, so much of the code is
parallelized.  The code will certainly run on a laptop, but a 32 or 64 core
server may be more appropriate.

## A Little History

Given a set of SNPs and some dependent variables, how should we determine which 
SNPs are significant as predictors?

At first blush, one might be tempted to use a univariate correlation analysis
on each SNP independently to determine if it is correlated with a predictor.
However, because there are multiple hypotheses in play (i.e., more than one SNP),
that doesn't work-- one risks "p-hacking" oneself.

An initial approach is to "correct" the p-value for the multiple hypotheses.
The simplest way to do this is to compute the univariate p-value, then
multiply the p-value by the number of hypotheses to "correct" it.  This is
called the "Bonferroni correction".  Unfortunately, there can easily be
thousands of hypotheses involved, so it is nearly impossible for a correlation
to appear significant (even if it really is).

An alternative, introduced in 1995, is the Benjamini-Hochberg procedure.  In
this approach, we construct a list of possibly significant SNPs, but control a
"false discovery rate" (FDR): we guarantee that the expected fraction of false
positives in our list is less than some target rate (like 10%).

The Benjamini-Hochberg procedure was a major advance in statistical handling
of situations with many hypotheses (as is the case with many SNPs).
However, there is ongoing interest in developing new techniques with more statistical
power; that is, two techniques may both guarantee a FDR less than 10%, but one
of them may produce more hypotheses than the other (i.e., is "more powerful").

A recent entrant into this field is the use of "statistical knockouts".  This
approach was introduced in [Panning for Gold: Model-X Knockoffs for High-
dimensional Controlled Variable Selection](https://arxiv.org/abs/1610.02351)
by Candes et al in 2016, and extended to hidden Markov models in [Gene Hunting
with Knockoffs for Hidden Markov Models](https://arxiv.org/abs/1706.04677) by
Sesia et al in 2017.  (Chromosomal crossover is well-modeled as a hidden Markov
model, so the latter set of techniques are appropriate.)

Sesia has open-sourced his code.  Our code largely acts to collect the relevant
pieces of data required by his code.

## Installation

For **Ubuntu 18.04 instance**.

```
sudo apt update
sudo apt -y install git
git clone https://github.com/wfbradley/snpko.git
sudo snpko/install.sh
```

Other **Linux** installs should be similar.

For **Mac/OS X**:  This has run this on a Mac under High Sierra.  In addition to 
the obvious install differences (e.g., replace `sudo apt get` with OS X equivalents), `pip install`ing SNPknock did not work as such.  Instead,
* git clone the SNPknock repo
* Find the `setup.py` file, and 
* Replace line 49 as follows:
    * Old: EXTRA_COMPILE_ARGS = ["-O3", "-std=c++11"]
    * New: EXTRA_COMPILE_ARGS = ["-O3", "-std=c++11", "-stdlib=libc++"]

## Quick Start Guide

You, the user, need to provide a file with SNPs and dependent variables for 
each of the people in your experiment.  An example is provided as `data/fake_SNP_data.csv`.

Given such a file, the entire pipeline can be run with 

```
python master_snpko.py --input my_SNP_data.csv
```

We provide a sample input file (with randomized data), so should be able to see the code in action by running
```
python master_snpko.py --input fake_SNP_data.csv
```

For the record, the command line for our original experimental data is:
```
python master_snpko.py --input_file 6.28.18.xlsx --skip_rows 1
```
(We cannot include the data file itself because of privacy concerns.)

## More Detailed Guide

The `master_snpko.py` module runs a series of individual that process the data in a series of stages.  The modules, in order, are:
*    **check_input**: Convert raw input data into a standardized form.
*    **ensembl_miner**: Query the ENSEMBL database for relevant genomic data about the SNPs, including genotypes of individuals.
*    **simple_stats**: Compute some naive univariate statistics with uncorrected p-values, along with Bonferroni-corrections.
*    **population_refiner**: Balance SNPs and population.  Not all individuals will have all SNPs sequenced, so we need to choose a subset of SNPs and a subset of the population so that both sets are relatively large.
*    **find_loci**: Remove correlated SNPs (i.e., deal with linkage disequilibrium.)
*    **make_knockoffs**:  Train hidden Markov Models for SNPs on each chromosome (with EM), and use each HMM to construct (multiple) knockoffs of the data.
*    **classifier**:  Run a classifier on the multiple knockoffs and determine which SNPs are significant predictors of which dependent variables given a target false discovery rate.

It is possible to run any one of these scripts individually; all take the same
command-line arguments (which are specified in `snpko_utils.py` and listed below in this document).  For example, the default false discovery rate is 0.1 (i.e., 10%), and you may want to 
experiment with several values.  The FDR only effects the classifier, so after running the code once, you can rerun just the classifier with a false discovery rate of 20% by running only the last step, e.g.:
```
python classifier.py --input my_SNP_data.csv --fdr 0.2
```

The code tries to cache partial results where possible; in particular, it will cache data from the ENSEMBL server (so it only needs to perform remote queries on the first pass, assuming the set of SNPs doesn't change) and it will cache the parameters from training the HMMs.

On the other hand, the entire set output (both intermediate files, cache files and results) exist in the same working directory (by default, `data/`).  If you want to preserve your current results in their entirety and start fresh, just move the directory:
```
mv data/ data_saved_1/
```

For easy reference, here is the full set of command-line arguments:

```
./master_snpko.py -h
usage: master_snpko.py [-h] [--input_file INPUT_FILE]
                       [--working_dir WORKING_DIR] [--skip_rows SKIP_ROWS]
                       [--na_threshold NA_THRESHOLD] [--never_na]
                       [--keep_old_logs] [--data_prefix DATA_PREFIX]
                       [--num_workers NUM_WORKERS] [--snp_weight SNP_WEIGHT]
                       [--fastPHASE_path FASTPHASE_PATH]
                       [--random_seed RANDOM_SEED]
                       [--num_knockoff_trials NUM_KNOCKOFF_TRIALS] [--verbose]
                       [--locus_threshold LOCUS_THRESHOLD] [--fdr FDR]
                       [--halt]

SNP Knockouts

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        CSV file with SNPs and dependent variables. (default:
                        None)
  --working_dir WORKING_DIR
                        Directory for all working files and final output.
                        (default: data)
  --skip_rows SKIP_ROWS
                        Skip rows of data file, 0-up indexing. (default: None)
  --na_threshold NA_THRESHOLD
                        If fraction of N/A entries in a row or column is >
                        this threshold, we discard it. "--skip_rows" is
                        applied first. First drop columns then rows. (default:
                        0.5)
  --never_na            By default, convert remaining N/As to double mutants;
                        if set, N/As raise exception. (default: False)
  --keep_old_logs       By default, old log file is deleted; this preserves
                        it. (default: False)
  --data_prefix DATA_PREFIX
                        All dependent variables should be labelled with a
                        common prefix (not beginning with "rs" or "r"), so we
                        can automatically detect them. (default: Imaging)
  --num_workers NUM_WORKERS
                        Number of parallel threads to use. (Default = number
                        of cores.) (default: -1)
  --snp_weight SNP_WEIGHT
                        Weight for Pareto-optimal tradeoff between population
                        and SNP count. (default: 2.0)
  --fastPHASE_path FASTPHASE_PATH
                        Path to "fastPHASE executable" (default: .)
  --random_seed RANDOM_SEED
                        Random seed for (reproducible) PRNGs (default: 123)
  --num_knockoff_trials NUM_KNOCKOFF_TRIALS
                        Because the knockoff process draws random samples, it
                        can be helpful to repeat it multiple times. (default:
                        100)
  --verbose             Enable verbose logging (debug level) (default: False)
  --locus_threshold LOCUS_THRESHOLD
                        Correlation threshold for declaring two SNPs to be in
                        the same locus. (default: 0.5)
  --fdr FDR             Target false discover rate (FDR). (default: 0.1)
  --halt                Enable verbose logging (debug level) (default: False)
```

## Running Time and Space

We benchmarked on a cloud instance running Ubuntu 18.04, with 96 cores, 10 GB of disk space and 480 GB of RAM.  (360 GB appeared insufficient.)  Our problem involved about 150 SNPs and 50 patients.  In that case, the `master_snpko` module took about 1 hour and 45 minutes to run.  Most of the time is spent in the last function, `classifier.py`.

10 GB was sufficient disk space for OS + temporary files. 

## Results

The script will produce a variety of output files in the `data/` directory and may take several hours to run.  The final output files will appear in `data/results/`.  In particular, output includes:
* `knockoff_trials.txt`: By default, we run 100 independent knockoffs for each experiment, and measure the percentage of knockoff trials in which a particular SNP shows up, for each label that we are predicting.  (For example, we might find that `rs12345` is a significant predictor for `symptom4`.)


## Author

This code was written by William (Bill) Bradley in 2017 and 2018.

## License

This project is licensed under the MIT License:

Copyright (c) 2018 William Bradley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Acknowledgements

This project was performed in collaboration with Dr. Michael S. Gee and Dr.
Cinthia C. Romero, both at Massachusetts General Hospital (MGH).

I would also like to thank [Matteo Sesia](http://web.stanford.edu/~msesia/)
for originally investigating the HMM case, for sharing his software tools with
the open source community, and for providing helpful advice.



