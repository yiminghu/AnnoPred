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
  parser.add_argument('--N_case', required=True, type=int,
                      help="Number of cases in GWAS training, for LDSC")
  parser.add_argument('--N_ctrl', required=True, type=int,
                      help="Number of ctrls in GWAS training, for LDSC")

  # For Pred
  parser.add_argument('--P', required=True, type=float,
                      help="Tuning parameter in (0,1)"
                           ", the proportion of causal snps")
  # Local LD file
  parser.add_argument('--local_ld_prefix', default="local_ld",
                      help="A local LD file name prefix"
                           ", will be created if not present")
  # Optional
  parser.add_argument('--num_iter', type=int, default=60
                      help="Number of iterations for mcmc," 
                           " If not specified, use 60")
  parser.add_argument('--ld_radius', type=int,
                      help="If not provided, will use the number of SNPs in" 
                           " common divided by 3000")
  parser.add_argument('--user_h2',
                      help="Path to user-provided heritability."
                           " If not provided, will use LDSC with 53 baseline"
                           " annotations, GenoCanyon and GenoSkyline.")
  ## Temporary file output directory
  parser.add_argument('--temp_dir', default=".",
                      help="Directory to output all temporary files."
                           " If not specified, will use the current directory.")
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
  pdict['N_case'] = args.N_case
  pdict['N_ctrl'] = args.N_ctrl

  if (args.P>0 and args.P<1):
    pdict['P'] = args.P
  else:
    exit("Tuning parameter needs to be in (0,1)!")

  if (isdir(args.temp_dir)):
    pdict['temp_dir'] = args.temp_dir
  else:
    exit("Directory for temporary files does not exist!")

  pdict['need_ld_radius'] = args.ld_radius is None
  pdict['ld_radius'] = args.ld_radius

  pdict['local_ld_prefix'] = args.local_ld_prefix
  
  pdict['num_iter'] = args.num_iter

  pdict['need_LDSC'] = args.user_h2 is None
  pdict['user_h2'] = args.user_h2
  if not pdict['need_LDSC'] and not isfile(args.user_h2):
    exit("User-provided H2 file does not exist!")

  pdict['out'] = args.out
  return pdict

# Returns the path to the file with name in the temp directory
def tmp(pdict, name):
  return join(pdict['temp_dir'], name)

# Returns pdict used by coord_trimmed
def pdict_coord_trimmed(pdict):
  d = {}
  d['N'] = pdict['N_case'] + pdict['N_ctrl']
  d['gf'] = pdict['ref_gt']
  d['vgf'] = pdict['val_gt']
  d['ssf'] = pdict['sumstats']
  d['ssf_format'] = 'BASIC'
  d['out'] = pdict['coord_out']
  d['gf_format'] = 'PLINK'
  d['skip_coordination'] = False
  d['vbim'] = pdict['val_gt']+'.bim'
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
  d['N'] = pdict['N_case'] + pdict['N_ctrl']
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
def pdict_pred_user(pdict, H2):
  d = pdict_pred_partial(pdict)
  d['hfile'] = None
  d['pfile'] = None
  d['H2'] = H2
  d['user_h2'] = pdict['user_hfile']
  return d

# Compute default LD radius
#def compute_default_ld_radius():
#  pass

def main(pdict):
  print(pdict)
  # Filter SNPs
  sumstats_filtered = tmp(pdict, "sumstats_filtered.txt")
  pre_sumstats.get_1000G_snps(pdict['sumstats'], 'ref/1000G_SNP_info.h5', sumstats_filtered)
  pdict['sumstats'] = sumstats_filtered

  # Generate coord_genotypes H5 file
  coord_trimmed.main(pdict_coord_trimmed(pdict))

  if pdict['need_ld_radius']:
    pdict['ld_radius'] = compute_default_ld_radius(???)

  if pdict['need_LDSC']:
    ldsc_result = LD_PyWrapper.callLDSC(
        pdict['sumstats'], pdict['N_case'], pdict[N_ctrl])
    pdict['h2file'] = tmp(pdict, "ldsc_h2_file.txt")
    pdict['pTfile'] = tmp(pdict, "ldsc_pT")
    res = prior_generating.generate_h2_pT(
        'ref/GC1_GS7_Baseline53.h5', 'ref/1000G_SNP_info.h5', pdict['coord_out'], ldsc_result, 
        pdict['h2file'], pdict['P'], pdict['pTfile'])
    pdict['pTfile'] = pdict['pTfile']+'_'+str(pdict['P'])+'_file.txt'
    if pdict['need_ld_radius']:
      pdict['ld_radius'] = res
    pred_main.main(pdict_pred_ldsc(pdict))
  else:
    pdict['user_hfile'] = tmp(pdict, "user_hfile.txt")
    res = prior_generating.generate_h2_from_user(
        pdict['user_h2'], pdict['coord_out'], pdict['user_hfile'])
    if pdict['need_ld_radius']:
      pdict['ld_radius'] = res[1]
    pred_main.main(pdict_pred_user(pdict, res[0]))


if __name__ == '__main__':
  args = get_argparser().parse_args()
  main(process_args(args))
