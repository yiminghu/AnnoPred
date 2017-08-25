## Dependency ##
pip install --user h5py
pip install --user plinkio
pip install --user scipy
pip install --user numpy

## install ##
git clone https://github.com/yiminghu/AnnoPred.git
cd AnnoPred
wget http://genocanyon.med.yale.edu/AnnoPredFiles/AnnoPred_ref.tar.gz
tar -zxvf AnnoPred_ref.tar.gz

LDSC_path="LDSC_path /your/path/to/ldsc" ## change to your ldsc path

## Setup LDSC ##
echo $LDSC_path > LDSC.config

## Input files and paths ##
sumstats_path="test_data/GWAS_sumstats.txt" ## change to sumstats path
individual_gt_path="test_data/test" ## change to individual level data (assuming using the same genotype file for reference and validation)
Ncase=12171
Nctrl=56862

## divide individual genotype data into two parts ##
## generate two list of individuals (randomly divide original genotype file into half)
## set your own paths and file names for these files
list_cv1="test_data/cv1.txt"
list_cv2="test_data/cv2.txt"
gt_cv1="test_data/cv1"
gt_cv2="test_data/cv2"

Rscript --vanilla split_cv.R $individual_gt_path".fam" $list_cv1 $list_cv2

plink --bfile $individual_gt_path --keep $list_cv1 --make-bed --out $gt_cv1
plink --bfile $individual_gt_path --keep $list_cv2 --make-bed --out $gt_cv2

## set your own paths and file names for these files
tmp_path1="tmp_files" ## including ldsc results and prior files for cv1
tmp_path2="tmp_files2" ## including ldsc results and prior files for cv2
results_output="res_output" ## including prs, phenotypes, aucs and betas for cv1 and 2
coord1=$tmp_path1"/coord"
coord2=$tmp_path2"/coord"

mkdir $tmp_path1
mkdir $tmp_path2
mkdir $results_output

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

## get average cv results
Rscript --vanilla results_cv.R $results_output"/cv1" $results_output"/cv2" "1.0 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 3e-05 1e-05"



