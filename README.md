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
* GWAS Summary statistics with a fixed format, for example:
* Reference/Validation genotype files, plink binary format, http://pngu.mgh.harvard.edu/~purcell/plink/

## Setup and Usage Example
1) Clone this repository
```
git clone https://github.com/yiminghu/AnnoPred.git
```
2) Download reference data
```
cd AnnoPred
wget http://genocanyon.med.yale.edu/AnnoPredFiles/AnnoPred_ref.tar.gz
tar -zxvf AnnoPred_ref.tar.gz
```
This step will generated a folder named ref containing functional annotations for AnnoPred.
3) Setup LDSC: open LDSC.config and change 
