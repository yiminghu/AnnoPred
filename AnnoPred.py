#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import isfile, isdir, join
from sys import exit

from annopred import prior_generating, coord_trimmed, pre_sumstats
from annopred import pred_main, LD_PyWrapper

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

  ## Parameters
  # For LDSC
  parser.add_argument('--N_sample', required=True, type=int,
                      help="Sample size of GWAS training, for LDSC")
  parser.add_argument('--annotation_flag', required=True,
                      help="Annotation flag: Tier0, Tier1, Tier2 and Tier3")

  # For Pred
  parser.add_argument('--P', required=True, type=float,
                      help="Tuning parameter in (0,1]"
                           ", the proportion of causal snps")
  # Local LD file
  parser.add_argument('--local_ld_prefix', required=True,
                      help="A local LD file name prefix"
                           ", will be created if not present")
  # Optional
  parser.add_argument('--ld_radius', type=int,
                      help="If not provided, will use the number of SNPs in" 
                           " common divided by 3000")
  parser.add_argument('--user_h2', 
                      help="Path to per-SNP heritability."
                           " If not provided, will use LDSC with 53 baseline"
                           " annotations, GenoCanyon and GenoSkyline.")
  ## Temporary file output directory
  parser.add_argument('--temp_dir', default=".",
                      help="Directory to output all temporary files."
                           " If not specified, will use the current directory.")
  parser.add_argument('--num_iter', type=int, default=60, 
                      help="Number of iterations for MCMC, default to 60.")
  ## Output Files
  #####################
  # Coord output H5 file 
  parser.add_argument('--coord_out', default="coord_out.h5", 
                      help="Output H5 File for coord_genotypes")
  # Output results
  parser.add_argument('--out', default="AnnoPred_out",
                      help="Output filename prefix for AnnoPred")

  return parser

# Check if all three files for PLINK exists
def check_plink_exist(prefix):
  suffices = ['bed', 'bim', 'fam']
  result = True
  for s in suffices:
    result = result and isfile(prefix + '.' + s)
  return result

# Validate Arguments in args and returns a dictionary
def process_args(args):
  pdict = {}
  
  # sumstats
  if (isfile(args.sumstats)):
    pdict['sumstats'] = args.sumstats
  else:
    exit("sumstats file does not exists!")

  # plink formats
  if (check_plink_exist(args.ref_gt)):
    pdict['ref_gt'] = args.ref_gt
  else:
    exit("Cannot find all reference genotype plink files!")
  if (check_plink_exist(args.val_gt)):
    pdict['val_gt'] = args.val_gt
  else:
    exit("Cannot find all validation genotype plink files!")

  pdict['coord_out'] = args.coord_out
  pdict['N'] = args.N_sample
  pdict['annotation_flag'] = args.annotation_flag

  if (args.P>0 and args.P<=1):
    pdict['P'] = args.P
  else:
    exit("Tuning parameter needs to be in (0,1)!")

  pdict['num_iter'] = args.num_iter

  if (isdir(args.temp_dir)):
    pdict['temp_dir'] = args.temp_dir
  else:
    exit("Directory for temporary files does not exist!")

  pdict['need_ld_radius'] = args.ld_radius is None
  pdict['ld_radius'] = args.ld_radius

  pdict['local_ld_prefix'] = args.local_ld_prefix

  pdict['need_LDSC'] = args.user_h2 is None
  pdict['user_h2'] = args.user_h2
  if not pdict['need_LDSC'] and not isfile(args.user_h2):
    exit("Per-SNP H2 file does not exist!")

  pdict['out'] = args.out
  return pdict

# Returns the path to the file with name in the temp directory
def tmp(pdict, name):
  return join(pdict['temp_dir'], name)

# Returns pdict used by coord_trimmed
def pdict_coord_trimmed(pdict):
  d = {}
  d['N'] = pdict['N']
  d['gf'] = pdict['ref_gt']
  d['vgf'] = pdict['val_gt']
  d['ssf'] = pdict['sumstats']
  d['ssf_format'] = 'BASIC'
  d['out'] = pdict['coord_out']
  d['gf_format'] = 'PLINK'
  d['skip_coordination'] = False
  d['vbim'] = None
  d['gmdir'] = None
  d['indiv_list'] = None
  return d

# Returns partially filled pdict used by Pred
def pdict_pred_partial(pdict):
  d = {}
  d['coord'] = pdict['coord_out']
  d['ld_radius'] = pdict['ld_radius']
  d['local_ld_file_prefix'] = pdict['local_ld_prefix']
  d['PS'] = pdict['P']
  d['num_iter'] = pdict['num_iter']
  d['N'] = pdict['N']
  d['out'] = pdict['out']
  return d

# Returns pdict used by Pred when using LDSC result
def pdict_pred_ldsc(pdict):
  d = pdict_pred_partial(pdict)
  d['hfile'] = pdict['h2file']
  d['pfile'] = pdict['pTfile']
  d['H2'] = None
  d['user_h2'] = None
  return d

# Returns pdict used by Pred when using user h2
def pdict_pred_user(pdict):
  d = pdict_pred_partial(pdict)
  d['hfile'] = None
  d['pfile'] = None
  d['H2'] = pdict['H2']
  d['user_h2'] = pdict['user_h2_trimmed']
  return d

def main(pdict):
  print(pdict)
  # Filter SNPs
  print('Filtering Summary Stats...')
  org_sumstats = pdict['sumstats']
  sumstats_filtered = tmp(pdict, "sumstats_filtered.txt")
  if not isfile(sumstats_filtered):
    pre_sumstats.get_1000G_snps(pdict['sumstats'], sumstats_filtered)
    pdict['sumstats'] = sumstats_filtered
  else:
    print('Filtered sumstats found, start coordinating genotypes...')

  # Generate coord_genotypes H5 file
  print('Coordinate summary stats and validation/reference genotype data...')
  if not isfile(pdict_coord_trimmed(pdict)['out']):
    coord_trimmed.main(pdict_coord_trimmed(pdict))
  else:
    print('Coord file already exists! Continue calculating priors...')

  if pdict['need_LDSC']:
    print('User-provided heritability file not found. Generating priors...')
#    if isfile()
    ldsc_result = tmp(pdict, pdict['annotation_flag'])
    if not isfile(ldsc_result+'_ldsc.results'):
      LD_PyWrapper.callLDSC(
          org_sumstats, pdict['N'], ldsc_result+'_ldsc', pdict['annotation_flag'])
    else:
      print('LDSC results found! Continue calculating priors ...')
    pdict['h2file'] = tmp(pdict, pdict['annotation_flag'] + "_ldsc_h2.txt")
    pdict['pTfile'] = tmp(pdict, pdict['annotation_flag'] + "_ldsc_pT"+str(pdict['P'])+".txt")
    ld_r = prior_generating.generate_h2_pT(
             pdict['coord_out'], ldsc_result+'_ldsc.results', 
             pdict['h2file'], pdict['P'], pdict['pTfile'], pdict['annotation_flag'])
    if pdict['need_ld_radius']: 
      pdict['ld_radius'] = int(ld_r)
    print('Starting AnnoPred...')
    pred_main.main(pdict_pred_ldsc(pdict))
  else:
    print('User-provided heritability file found. Extracting SNPs in common...')
    pdict['user_h2_trimmed'] = tmp(pdict, "user_h2_trimmed.txt")
    pdict['H2'], ld_r = prior_generating.generate_h2_from_user(
           pdict['user_h2'], pdict['coord_out'], pdict['user_h2_trimmed'])
    if pdict['need_ld_radius']:
      pdict['ld_radius'] = int(ld_r)
    print('Starting AnnoPred...')
    pred_main.main(pdict_pred_user(pdict))


if __name__ == '__main__':
  args = get_argparser().parse_args()
  main(process_args(args))
