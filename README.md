# A Pipeline for Statistical Knockoffs with SNPs

##### Table of Contents  
-  [The Problem](#the-problem)  
    *  [The General Problem](#the-general-problem)  
    *  [The Specific Problem](#the-specific-problem)  
    *  [A Little History](#a-little-history)  
-  [Installation](#installation)  
-  [Quick Start Guide](#quick-start-guide)  
-  [More Detailed Guide](#more-detailed-guide)  
    *  [General Structure](#general-structure)
    *  [Logging](#logging)
    *  [Rerunning the Code](#rerunning-the-code)
    *  [Command Line Arguments](#command-line-arguments)
    *  [Halting on Completion](#halting-on-completion)
    *  [Running Time and Space](#running-time-and-space)  
    *  [Results](#results)  
    *  [P-values for the Selection Frequency](#p-values-for-the-selection-frequency)
    *  [Regression Tests](#regression-tests)
-  [Author](#author)  
-  [License](#license)  
-  [Acknowledgements](#acknowledgements)  


## The Problem

### The General Problem

Suppose we have selected a set of single nucleotide polymorphisms (SNPs) along
with a set of binary observations (such as the presence or absence of a
particular symptom). We measure these SNPs and observations across a
population of (human) subjects. Our goal is to determine which SNPs are
predictive of which observations.

As one can easily consider hundreds or thousands of SNPs, we are faced with a
multiple hypothesis problem.  In particular, forming a 2x2 contingency table
for each (SNP, observation) pair is statistically unsound-- we would consider
so many tables that we would "accidentally" see some correlation.

This software provides tools for answering this type of question using a new
statistical framework called a "knockoff", which appears to be a major advance
over classical tools for handling multiple hypotheses. Given any set of SNPs
and labels (i.e., dependent variables), these tools will:
* Download relevant information about the SNPs from the ENSEMBL online repository, 
* Train an HMM on the data,
* Construct statistical knockoffs,
* Train a classifier on the data, and
* Produce a set of SNPs with a controlled false discovery rate.

Some of these steps are computationally intensive, so much of the code is
parallelized.  The code will certainly run on a laptop, but a 32 or 64 core
server may be more appropriate.

### The Specific Problem

We developed this general pipeline to address a specific instance of the problem.

We examined a population of 54 pediatric patients suffering from Crohn's
disease. For each patient, we measured 168 SNPs and observed 8 different
radiological imaging features (such as "Lumen narrowing" or "Small bowel
disease").  We were interested to determine which SNPs, if any, were
predictors of particular imaging features.

This Python module includes all the statistics used in our paper:
[Novel Associations Between Genome‐Wide Single Nucleotide Polymorphisms and MR Enterography Features in Crohn's Disease Patients](https://onlinelibrary.wiley.com/doi/abs/10.1002/jmri.27250).


### A Little History

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
of situations with many hypotheses (as is the case with many SNPs). However,
there is ongoing interest in developing new techniques with more statistical
power; that is, two techniques may both guarantee a FDR less than 10%, but one
of them may produce more hypotheses than the other (i.e., is "more powerful").

A recent entrant into this field is the use of "statistical knockoffs".  This
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

If you need to install `git` and then clone this git repo, run:
```
sudo apt update
sudo apt upgrade
sudo apt -y install git
git clone https://github.com/wfbradley/snpko.git
```
Then to build and install the dependencies for this repo itself, run:

```
cd snpko
sudo ./install.sh
```
Optionally, to enable the module to halt the machine, see [Halting on
Completion](#halting-on-completion).

Other **Linux** installs should be similar.

For **Mac/OS X**:  We have run this software on a Mac under High Sierra.  We
only sketch details for installation.  In addition to  the obvious install
differences (e.g., replace `sudo apt get` with OS X equivalents), `pip
install`ing SNPknock did not work as such.  Instead,
* git clone the `SNPknock` repo, as above.
* Find the `setup.py` file in the `SNPknock` rep.
* Replace line 49 as follows:
    * Old: EXTRA_COMPILE_ARGS = ["-O3", "-std=c++11"]
    * New: EXTRA_COMPILE_ARGS = ["-O3", "-std=c++11", "-stdlib=libc++"]
* Run `setup.py` to install the `SNPknock` module.

## Quick Start Guide

You, the user, need to provide a file with SNPs and dependent variables for
each of the people in your experiment.  An example is provided as
`data/fake_SNP_data.csv`.

Given such a file, the entire pipeline can be run with 

```
./master_snpko.py --input my_SNP_data.csv
```

We provide a sample input file (with randomized data), so should be able to see the code in action by running
```
./master_snpko.py --input fake_SNP_data.csv
```

For the record, the command line for our original experimental data is:
```
./master_snpko.py --input_file 6.28.18.xlsx --skip_rows 1
```
(We cannot include the data file itself because of privacy concerns.)

## More Detailed Guide

### General Structure

The `master_snpko.py` module runs a series of individual that process the data in a series of stages.  The modules, in order, are:
*    **check_input**: Convert raw input data into a standardized form.
*    **ensembl_miner**: Query the ENSEMBL database for relevant genomic data about the SNPs, including genotypes of individuals.
*    **simple_stats**: Compute some naive univariate statistics with uncorrected p-values, along with Bonferroni-corrections.
*    **population_refiner**: Balance SNPs and population.  Not all individuals will have all SNPs sequenced, so we need to choose a subset of SNPs and a subset of the population so that both sets are relatively large.
*    **find_loci**: Remove correlated SNPs (i.e., deal with linkage disequilibrium.)
*    **make_knockoffs**:  Train hidden Markov Models for SNPs on each chromosome (with EM), and use each HMM to construct (multiple) knockoffs of the data.
*    **classifier**:  Run a classifier on the multiple knockoffs and determine which SNPs are significant predictors of which dependent variables given a target false discovery rate.
*    **sig_results**:  Filter the results to the relevant ones (i.e., the significant SNPs).

It is possible to run any one of these scripts individually; all take the same
command-line arguments (which are specified in `snpko_utils.py` and listed below in this document).

### Logging

Detailed logging information is written in the working directory to `run.log`.  By default, this file is located at `data/run.log`.

### Rerunning the Code

You may wish to rerun the code in part or in whole.  You'll need to run the `master.py` once to produce all the intermediate files, but from that point, you can rerun portions selectively.

For example, the default false discovery rate is 0.1 (i.e., 10%), but you may
want to experiment with several values.   Since FDR only effects the
classifier, so after running the code once, you can rerun just the classifier
with a false discovery rate of 20% by running only the last step, e.g.:
```
python classifier.py --input my_SNP_data.csv --fdr 0.2
```

To accelerate computation during subsequent reprocessing, the code caches
partial results where possible. In particular, it will cache data from the
ENSEMBL server (so it only needs to perform remote queries on the first pass,
assuming the set of SNPs doesn't change) and it will cache the HMM parameters
fit from the EM run.

Note that this process will overwrite old files in the working directory including
(unless you use the `--keep_old_logs` flag) the old `run.log` file.

If you wish to preserve all your data as is start afresh, just
move the working directory to a new location, e.g.,:
```
mv data/ data_saved_1/
```

### Command Line Arguments
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

SNP Knockoffs

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
  --obs_freq OBS_FREQ   Only trust SNPs that show up in >obs_freq of the
                        knockoff trials. (default: 0.5)  
  --halt                Enable verbose logging (debug level) (default: False)
```

### Halting on Completion

The following situation can occur: We run the module on some data, expecting it to take many hours to complete.  We use on a large, relatively expensive cloud instance to speed up the computation.  The code completes in the middle of the night.  We would really like the machine to shut itself down at this point so we do not need to keep paying for an idle instance.  In that situation, use the `--halt` flag, like so:
```
sudo ./master_snpko.py --input my_SNP_data.csv --halt
```
(Note the `sudo`.)  However, there is a complexity: halting an instance requires superuser access, and `sudo` privileges only lasts for 15 minutes (by default).  So by the time the module finishes running (many hours after launch) and tries to shut itself down, it will no longer be allowed to do so.

This problem can be addressed on Ubuntu 18.04 (and similar platforms) as follows:
```
# Install your favorite editor, e.g.,
sudo apt-get install emacs
# Set it as default
export VISUAL=emacs
# Tinker with sudoers file to make sudo last longer
sudo visudo
```
This will launch your editor on the /etc/sudoers file.  (It is preferable to use `visudo`, rather than editing directly, because `visudo` will guarantee that your edited file parses correctly.)  You will find a bunch of lines like:
```
Defaults        env_reset
Defaults        mail_badpass
...
```
Add a new line:
```
Defaults        timestamp_timeout=-1
```
and save the change.  This change causes sudo privilege never to expire.  There is an obvious security risk associated with this, although if an adversary can execute commands with `sudo` priveleges, the time limit is probably the least of your worries.

Finally, reboot your system:
```
sudo shutdown -r now
```
Once the instance restarts, log in again, and proceed to run something like:
```
sudo nohup ./master_snpko.py --input my_SNP_data.csv --halt &
```

Note that if a submodule raises an exception, this will be caught and the machine will still halt.  This is usually a feature (e.g., if your script fails after 3 hours, you still want to stop the instance), but it means that if you have an initial error (e.g., you mistype your input file name), your instance will immediately stop.  So it is probably wise to make sure your command runs for a few minutes *without* the `--halt` option first.  If it seems to be running successfully, then terminate the job (CTRL-C), and restart *with* the `--halt` flag.

### Running Time and Space

Much of of this code is parallelized, so you will benefit from running on a multi-core machine.

We benchmarked on a cloud instance running Ubuntu 18.04, with 96 cores, 10 GB of disk space and 480 GB of RAM.  (360 GB appeared insufficient.)  Our problem involved about 150 SNPs and 50 patients.  In that case, the `master_snpko` module took about 1 hour and 45 minutes to run.  Most of the time is spent in the last function, `classifier.py`; the combined running time for all other sections was only 3 minutes.

10 GB was sufficient disk space for OS + temporary files. 

### Results

When the module runs, it produces a variety of files in the working directory (default: `data/`) that may be of interest.  On completion, final output is written to a `results` subdirectory (default: `data/results/`).  Output files are:
* `knockoff_trials.txt`:  This file records the frequency with which particular SNPs were correlated with particular dependent variables.  The file should be fairly human-readable.  It has a section for each label (i.e., dependent variable); within a section, it has a section for the modified FDR (`mFDR`), followed by the classical FDR (`cFDR`).  Within a subsection, we list the percentage of knockoff trials for which the SNP appeared (suppressing the SNP if it never showed up).  For example:
```
Label: Imaging: Colon disease
Type of FDR: mFDR
   rs6088765 : 73%
   rs11742570 : 69%
   rs174537 : 43%
```
* `uncorrected.csv`: This is a convenience file listing the *uncorrected* p-values and odds ratios for all the <SNP, label> pairs, along with the (statistically weak) Bonferroni corrected p-values.  Here are the CSV fields with one sample row:
```
SNP,label,uncorrected_p_value,uncorrected_odds_ratio,bonferroni_corrected_p_value
rs10486483,Imaging: Activity,0.009536,9.285714,1.0
```
* `exploratory.csv`: This is the subset of <SNP, label>s from `uncorrected.csv` where the uncorrected p-values are <0.05.  Although these p-values are not reliable in themselves, they may be useful from the standpoint of exploratory statistics to help guide future studies.
* `sig_results.csv`: This is the the subset of <SNP, label> pairs that the knockoff procedure has deemed to be significant.  Specifically, we require the mFDR to be > `--fdr` (default: 0.1=10%); and we require this threshold to be crossed for > `--obs_freq` of the knockoff trials (default: 0.5=50%).  Here are the CSV fields with one sample row:
```
SNP,fdr,fdr_type,label,obs_freq,uncorrected_odds_ratio,uncorrected_p_value
rs6088765,0.1,mFDR,Imaging: Colon disease,0.73,6.111111,0.083833
```

By default, we run 100 independent knockoffs for each experiment, and measure the percentage of knockoff trials in which a particular SNP shows up, for each label that we are predicting.  (For example, we might find that `rs12345` is a significant predictor for `symptom4` in 37 of the .)

### P-values for the Selection Frequency
 
 Running "master_snpko.py" will run a series of knockoff trials aimed to
 generate lists of SNPs with a target FDR.  For a particular SNP and FDR
 threshold, then, we might discover that it appears in, say, 85% of all
 knock-off runs.  Let's call this fraction the *selection frequency*, and denote it by *Q*.

 So, should we consider a selection frequency of 85% to be significant?

 To address this question, we can do the following.
```
     We examine the experimental data to determine
         The SNPs of interest,
         The number of people in the trials (say, X subjects), and
         The labels for those people
     We pull background data from ENSEMBL, with SNPs from Y>>X people.
     for trial=1, ..., T=100:
         Randomly partition ENSEMBLE into X subjects vs Y-X subjects.
         Replace the experimental SNPs with the SNPs from the X
             ENSEMBL subjects
         Restrict ENSEMBL data to the Y-X subjects.
         Proceed as usual with SNPKO.
         For each label (e.g., "Imaging: Fistula"), each SNP shows up
             in some fraction of the FDR lists.  Record the maximum,
             across all SNPS, of Q.  I.e., for a given label,
                 max_{SNP} Q_{SNP, label}
```

 Once all T=100 trials are done, we now have a distribution for the maximum
 selection frequency Q for the null hypothesis (namely, that there is no correlation between the
 SNPs and the binary output).  So, for example, if we're interested in a
 p=0.05 significance level, we consider all <SNP,label> pairs that are
 greater than or equal to the 95th largest value of Q.

 The `snpko` module includes code for computing p-values for selection frequency.  Note that this
 is on the order of 100x slower than the regular SNPKO processing, so is best done by splitting
 the work across multiple large machines.  We also support sharing the data via Gcloud or AWS
 by command line arguments to keep the bookkeeping simpler.


### Regression Tests

To confirm that the code is functioning correctly, we provide a regression test.  The test generates some synthetic data with known correlations, runs `master_snpko.master()`, and confirms that the correct SNPs were recovered.  It can be run as follows:

```
cd snpko/  # change to main SNPKO directory, wherever that is
python test/test_snpko.py
```

On a 2017 laptop with 4 cores, the test took about 15 minutes.  Results are written to STDOUT and to `/tmp/test_snpko/run.log`.  Success ends with "Test passed successfully"; failure should throw an exception.

Sometimes it can be difficult to debug a problem because the data may be too sensitive to share, but to reproduce the problem we need to mimic the structure of the input file.  To address that problem, we also provide `tests/anonymize_data.py`.  This script produces an "anonymized" version of a target input file in which the entries of each row are scrambled.  Be aware that the original data is still present (e.g., if patients' names were present, they will still be present, just in a random order.)

## Author

This code was written by William (Bill) Bradley in 2017 and 2018.

## License

This project is licensed under the MIT License.  See `LICENSE` file for the complete license.

Copyright (c) 2018 William Bradley


## Acknowledgements

This project was performed in collaboration with Dr. Michael S. Gee and Dr.
Cinthia C. Romero, both at Massachusetts General Hospital (MGH).

I would also like to thank [Matteo Sesia](http://web.stanford.edu/~msesia/)
for originally investigating the HMM case, for sharing his software tools with
the open source community, and for providing helpful advice.



