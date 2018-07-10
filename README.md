# A Pipeline for Statistical Knockouts with SNPs

##### Table of Contents  
[The Problem](#the-problem)  
[A Little History](#a-little-history)  
[Installation](#installation)  
[Author](#author)  
[License](#license)  
[Acknowledgements](#acknowledgements)  


## The Problem

We observed a population of 54 pediatric patients suffering from Crohn's disease.
For each patient, we measured 168 SNPs and observed 8 different radiological
imaging features (such as "Lumen narrowing" or "Small bowel disease").  We
were interested to determine which SNPs, if any, were significant predictors
of particular imaging features.

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

These installation instructions worked on an Ubuntu 18.04 instance (10GB of disk space).

```
sudo apt update
sudo apt -y install git
git clone https://github.com/wfbradley/snpko.git
cd snpko
sudo ./install.sh
```

Install and compile fastPhase...

## Quick Start Guide

You, the user, need to provide a file with SNPs and dependent variables for 
each of the people in your experiment.  An example is provided as `data/fake_SNP_data.csv`.

Given such a file, the entire pipeline can be run with 

```
python master.py --input data/fake_SNP_data.csv
```

On an experiment with 150 SNPs and 50 patients, running on an instance with 96 cores, the script
takes about 15-20 minutes to run.  Most of the time is in the last function, `classifier.py`

The script will produce a variety of output files in the `data/` directory and may take several hours to run.  The final output files will appear in `data/results/`.  In particular, output includes:
* `knockoff_trials.txt`: By default, we run 100 independent knockoffs for each experiment, and measure the percentage of knockoff trials in which a particular SNP shows up, for each label that we are predicting.  (For example, we might find that `rs12345` is a significant predictor for `symptom4`.)

By examining `master.py`, you will see that there are a series of individual scripts run in
order.  It is possible to run any one of these scripts individually; all take the same
command-line arguments (which are specified in `snpko_utils.py`)

For the record, the command for our original experimental data (which we cannot include :
```
python master.py --input_file 6.28.18.xlsx --skip_rows 1
```

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



