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
```
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
Keep in mind that this will generate intermediate data at tmp_test/ and test_output/. Contents in These folders are reused on different runs, not deleted: you might want to delete this folder before running AnnoPred on a new dataset, or specify a different folder on each run.

The example command parameters mean:
* --sumstats=GWAS_sumstats.txt: GWAS summary statistics, with seven fields: hg19chrc, snpid, a1, a2, bp, or and p
* --ref_gt=validation: path to the reference genotype data. We suggest also using validation data as reference data. Plink binary format (.bed, .bim, .fam), description can be * found in .
* --val_gt=validation: path to the validation genotype data. Plink binary format, the sixth column in fam file cannot be missing.
* --coord_out=test_output/coord_out: path for saving a h5py file, which contains validation genotypes, summary statistics and standardized effect sizes of SNPs in common.
* --N_case=12171: number of cases in GWAS
* --N_ctrl=56862: number of controls in GWAS
* --P=0.1: pre-spesified parameter, the proportion of causal variants
* --local_ld_prefix=tmp_test/local_ld: a path for saving a cPickle file, which contains LD matrix
* --out=test_output/test: a path for output files
* --temp_dir=tmp_t2d: a path for saving temporary files generated during the procedure. We suggest using different temp_dir for different dataset.

