# AnnoPred

## Introduction
This tool predicts disease risk from genotype data using large GWAS summary statistics as training data and integrating functional annotations. The online version of manuscript can be found at http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005589.

For a quick start, you may run the pipeline.sh!

## Prerequisites
The software is developed and tested in Linux. You will need Python 2.7 and several pacakges to run it:
* h5py
* plinkio
* scipy
* numpy

To install these packages, you can conveniently use **pip**. For example,
```
pip install h5py
```
or 
```
pip install --user h5py
```

Besides these, you also need to have [LDSC](https://github.com/bulik/ldsc) installed (LDSC itself also has a list of prerequisites, please make sure they are also installed). If you only want to use your own heritability estimation for each SNP, you can skip this.

## Input Data
* GWAS Summary statistics with a fixed format, for example: test_data/GWAS_sumstats.txt
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

3) Setup LDSC: open LDSC.config and change /absolute/path/to/ldsc to the absolute path to LDSC in your local directory or use
```
LDSC_path="LDSC_path /your/path/to/ldsc" ## change to your ldsc path
echo $LDSC_path > LDSC.config
```
Instruction on installing LDSC can be found at https://github.com/bulik/ldsc.

4) Example (when heritability estimation not provided):
```
mkdir test_output ## create dir for output files
mkdir tmp_test ## create dir for temporary files
python AnnoPred.py\
  --sumstats=test_data/GWAS_sumstats.txt\
  --ref_gt=test_data/test\
  --val_gt=test_data/test\
  --coord_out=test_output/coord_out\
  --N_case=12171\
  --N_ctrl=56862\
  --P=0.1\
  --local_ld_prefix=tmp_test/local_ld\
  --out=test_output/test\
  --temp_dir=tmp_test\
```
Keep in mind that this will generate intermediate data at tmp_test/ and test_output/. Contents in These folders are reused on different runs, not deleted: you might want to delete this folder before running AnnoPred on a new dataset, or specify a different folder on each run.

The example command parameters mean:
* --sumstats=test_data/GWAS_sumstats.txt: GWAS summary statistics, with seven fields: hg19chrc, snpid, a1, a2, bp, or and p. test_data/GWAS_sumstats.txt is a subset of DIAGRAM summary statistics. We thank DIAGRAM consortium for making the data publicly available. The oringinal download link is http://diagram-consortium.org/downloads.html
* --ref_gt=test_data/test: path to the reference genotype data. We suggest also using validation data as reference data. Plink binary format (.bed, .bim, .fam), description can be * found in .
* --val_gt=test_data/test: path to the validation genotype data. Plink binary format, the sixth column in fam file cannot be missing.
* --coord_out=test_output/coord_out: path for saving a h5py file, which contains validation genotypes, summary statistics and standardized effect sizes of SNPs in common.
* --N_case=12171: number of cases in GWAS
* --N_ctrl=56862: number of controls in GWAS
* --P=0.1: pre-spesified parameter, the proportion of causal variants
* --local_ld_prefix=tmp_test/local_ld: a path for saving a cPickle file, which contains LD matrix
* --out=test_output/test: a path for output files
* --temp_dir=tmp_test: a path for saving temporary files generated during the procedure. We suggest using different temp_dir for different dataset.

5) Example (when heritability estimation provided):
```
python AnnoPred.py\
  --sumstats=test_data/GWAS_sumstats.txt\
  --ref_gt=test_data/test\
  --val_gt=test_data/test\
  --coord_out=test_output/coord_out\
  --N_case=12171\
  --N_ctrl=56862\
  --P=0.1\
  --user_h2=test_data/user_h2_est.txt
  --local_ld_prefix=tmp_test/local_ld\
  --out=test_output/test\
  --temp_dir=tmp_test\
```
* --user_h2=test_data/user_h2_est.txt: user-provided heritability estimation for each SNP. A text file with three fields: chr, snpid and per-snp heritability estimation. See test_data/user_h2_est.txt for example.

## Output files
When --user_h2 is not provided, AnnoPred will output a set of files including two types of AnnoPred polygenic risk scores, phenotypes of testing data, prediction accuracy, posterior expectation estimation of the effect size of each snp. Take the output of code shown in 4) of previous section as an example:
* test_output/test_h2_non_inf_y_0.1.txt: phenotypes of testing data.
* test_output/test_h2_non_inf_prs_0.1.txt: AnnoPred PRS using the first type of priors (see manuscript for details).
* test_output/test_h2_non_inf_auc_0.1.txt: prediction accuracy of AnnoPred PRS using the first prior: AUC for binary traits and correlation between PRS and y for continuous traits.
* test_output/test_h2_non_inf_betas_0.1.txt: posterior expectation estimation of the effect size of each snps using the second type of prior.
* test_output/test_pT_non_inf_prs_0.1.txt: AnnoPred PRS using the second type of priors (see manuscript for details).
* test_output/test_pT_non_inf_auc_0.1.txt: prediction accuracy of AnnoPred PRS using the first prior: AUC for binary traits and correlation between PRS and y for continuous traits.
* test_output/test_pT_non_inf_betas_0.1.txt: posterior expectation estimation of the effect size of each snps using the second type of prior.

## A pipeline for setting up cross-validation and get the average prediction accuracy (using the test data as an example)
### Input files and paths
```
sumstats_path="test_data/GWAS_sumstats.txt" ## change to sumstats path
individual_gt_path="test_data/test" ## change to individual level data (assuming using the same genotype file for reference and validation)
Ncase=12171
Nctrl=56862
```
### Divide individual genotype data into two parts
```
## set your own paths and file names for these files
list_cv1="test_data/cv1.txt"
list_cv2="test_data/cv2.txt"
gt_cv1="test_data/cv1"
gt_cv2="test_data/cv2"

Rscript --vanilla split_cv.R $individual_gt_path".fam" $list_cv1 $list_cv2

plink --bfile $individual_gt_path --keep $list_cv1 --make-bed --out $gt_cv1
plink --bfile $individual_gt_path --keep $list_cv2 --make-bed --out $gt_cv2
```
### Set your own paths and file names for temporary and output files
```
tmp_path1="tmp_files" ## including ldsc results and prior files for cv1
tmp_path2="tmp_files2" ## including ldsc results and prior files for cv2
results_output="res_output" ## including prs, phenotypes, aucs and betas for cv1 and 2
coord1=$tmp_path1"/coord"
coord2=$tmp_path2"/coord"

mkdir $tmp_path1
mkdir $tmp_path2
mkdir $results_output
```
### Run AnnoPred on different cv with a sequence of tuning parameters
```
## could be paralellized
for p in 1.0 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 3e-05 1e-05
do
	python AnnoPred.py\
	  --sumstats=$sumstats_path\
	  --ref_gt=$gt_cv1\
	  --val_gt=$gt_cv1\
	  --coord_out=$coord1\
	  --N_case=$Ncase\
	  --N_ctrl=$Nctrl\
	  --P=$p\
	  --local_ld_prefix=$tmp_path1"/local_ld"\
	  --out=$results_output"/cv1"\
	  --temp_dir=$tmp_path1
done

for p in 1.0 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 3e-05 1e-05
do
	python AnnoPred.py\
	  --sumstats=$sumstats_path\
	  --ref_gt=$gt_cv2\
	  --val_gt=$gt_cv2\
	  --coord_out=$coord2\
	  --N_case=$Ncase\
	  --N_ctrl=$Nctrl\
	  --P=$p\
	  --local_ld_prefix=$tmp_path2"/local_ld"\
	  --out=$results_output"/cv2"\
	  --temp_dir=$tmp_path2
done
```
### Get average cv results
```
Rscript --vanilla results_cv.R $results_output"/cv1" $results_output"/cv2" "1.0 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 3e-05 1e-05"
```

## Acknowledgement
**Part of the code is modified from LDpred (https://bitbucket.org/bjarni_vilhjalmsson/ldpred). We thank Dr. Bjarni J. Vilhjalmsson for sharing his code.**


