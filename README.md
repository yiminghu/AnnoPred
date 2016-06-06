## AnnoPred

Part of the code is modified from LDpred (https://bitbucket.org/bjarni_vilhjalmsson/ldpred). We thank Dr. Bjarni J. Vilhjalmsson for sharing his code.

# 1. Introduction
This tool predict disease risk from genotype data using large GWAS summary statistics as training data and integrating functional annotations.

# 2. Prerequisites
The software is developed and tested in Linux. You will need Python 2.7 and several pacakges to run it:
* h5py
* plinkio
* scipy
* numpy
Besides these, you also need to have [LDSC](https://github.com/bulik/ldsc) installed. 

3. Input Data
GWAS Summary statistics with a fixed format, for example:
reference/validation genotype files, plink binary format, http://pngu.mgh.harvard.edu/~purcell/plink/

4. Setup and Usage Example
Clone this repository
'''
git clone
'''
