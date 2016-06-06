# AnnoPred

Part of the code is modified from LDpred (https://bitbucket.org/bjarni_vilhjalmsson/ldpred). We thank Dr. Bjarni J. Vilhjalmsson for sharing his code.

## Introduction
This tool predict disease risk from genotype data using large GWAS summary statistics as training data and integrating functional annotations.

## Prerequisites
The software is developed and tested in Linux. You will need Python 2.7 and several pacakges to run it:
* h5py
* plinkio
* scipy
* numpy

Besides these, you also need to have [LDSC](https://github.com/bulik/ldsc) installed. 

## Input Data
GWAS Summary statistics with a fixed format, for example:
reference/validation genotype files, plink binary format, http://pngu.mgh.harvard.edu/~purcell/plink/

## Setup and Usage Example
Clone this repository

'''
git clone https://github.com/yiminghu/AnnoPred.git
'''
