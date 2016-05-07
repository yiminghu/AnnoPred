#!/usr/bin/env python

from argparse import ArgumentParser
from annopred import pred_main, LD_PyWrapper, prior_generating, coord_trimmed

# Create the master argparser and returns the argparser object
def get_argparser():
  parser = ArgumentParser(prog="AnnoPred", 
                          description="Genetic Risk Prediction Software")
  ## Input Files
  #################### 
  # GWAS sumstats
  parser.add_argument('--sumstats', required=True, help="GWAS summary stats")
  # Reference Genotype
  parser.add_argument('--ref_gt', required=True, 
                      help="Reference genotype, plink bed format")
  # Validation Genotype
  parser.add_argument('--val_gt', required=True, 
                      help="Validation genotype, plink bed format")

  # For LDSC
  parser.add_argument('--N_case', required=True, type=int,
                      help="Number of cases in GWAS training, for LDSC")
  parser.add_argument('--N_ctrl', required=True, type=int,
                      help="Number of ctrls in GWAS training, for LDSC")

  # For Pred
  parser.add_argument('--P', required=True, type=float,
                      help="Tuning parameter in (0,1)"
                           ", the proportion of causal snps")
  # Local LD file
  parser.add_argument('--local_ld_prefix', default="test_local_ld",
                      help="A local LD file name prefix"
                           ", will be created if not present")
  # Optional
  parser.add_argument('--ld_radius', type=int,
                      help="If not provided, will use the number of SNPs in" 
                           " common divided by 3000")
  parser.add_argument('--per_SNP_h2', 
                      help="Path to per-SNP heritability."
                           " If not provided, will use LDSC with 53 baseline"
                           " annotations, GenoCanyon and GenoSkyline.")

  ## Output Files
  #####################
  # Coord output H5 file 
  parser.add_argument('--coord_out', default="coord_out.h5", 
                      help="Output H5 File for coord_genotypes")
  # Output results
  parser.add_argument('--out', default="AnnoPred_out",
                      help="Output filename prefix for AnnoPred")

  return parser


if __name__ == '__main__':
  args = get_argparser().parse_args()
  print(args)
  print(pred_main.main)
  print(LD_PyWrapper.callLDSC)
  print(prior_generating.generate_h2_pT)
  print(prior_generating.generate_h2_from_user)
  print(coord_trimmed.main)
