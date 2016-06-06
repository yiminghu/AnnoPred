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

3) Setup LDSC: open LDSC.config and change /absolute/path/to/ldsc to the absolute path to LDSC in your local directory

4) Example:
```python
python AnnoPred.py\
  --sumstats=GWAS_sumstats.txt\
  --ref_gt=validation\
  --val_gt=validation\
  --coord_out=test_output/coord_out\
  --N_case=12171\
  --N_ctrl=56862\
  --P=0.1\
  --local_ld_prefix=tmp_test/local_ld\
  --out=test_output/test\
  --temp_dir=tmp_t2d\
```
