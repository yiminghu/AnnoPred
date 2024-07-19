## This script was modified from Dr. Bjarni J. Vilhjalmsson's code (https://bitbucket.org/bjarni_vilhjalmsson/ldpred). We thank him for sharing his code

try: 
    import scipy as sp
except Exception:
    print('Using Numpy instead of Scipy.')
    import numpy as sp
    
#from numpy import linalg 
from scipy import linalg 
import pdb
import plinkio
from plinkio import plinkfile
import random

import time
import gzip
import itertools as it
from sklearn import metrics

import getopt
import sys
import traceback
import time
import os
import gzip
import itertools as it

import h5py
import scipy as sp
from scipy import stats
import pickle
from sklearn import metrics


chromosomes_list = ['chrom_%d'%(x) for x in range(1,23)]
chromosomes_list.append('chrom_X')


def pred_accuracy(y_true, y_pred):
    y_true = sp.copy(y_true)
    if len(sp.unique(y_true))==2:
        print('dichotomous trait, calculating AUC')
        y_min = y_true.min()
        y_max = y_true.max()
        if y_min!= 0 or y_max!=1:
            y_true[y_true==y_min]=0
            y_true[y_true==y_max]=1
        fpr, tpr, thresholds = metrics.roc_curve(y_true, y_pred)
        auc = metrics.auc(fpr, tpr)
        return auc
    else:
        print('continuous trait, calculating COR')
        cor = sp.corrcoef(y_true,y_pred)[0,1]
        return cor


def get_LDpred_ld_tables(snps, ld_radius=100, ld_window_size=0):
    """
    Calculates LD tables, and the LD score in one go...
    """
    
    ld_dict = {}
    m,n = snps.shape
    print(m,n)
    ld_scores = sp.ones(m)
    ret_dict = {}
    for snp_i, snp in enumerate(snps):
        # Calculate D
        start_i = max(0, snp_i - ld_radius)
        stop_i = min(m, snp_i + ld_radius + 1)
        X = snps[start_i: stop_i]
        D_i = sp.dot(snp, X.T) / n
        r2s = D_i ** 2
        ld_dict[snp_i] = D_i
        lds_i = sp.sum(r2s - (1-r2s) / (n-2),dtype='float32')
        #lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
        ld_scores[snp_i] =lds_i
    
    ret_dict['ld_dict']=ld_dict
    ret_dict['ld_scores']=ld_scores
    
    if ld_window_size>0:
        ref_ld_matrices = []
        for i, wi in enumerate(range(0, m, ld_window_size)):
            start_i = wi
            stop_i = min(m, wi + ld_window_size)
            curr_window_size = stop_i - start_i
            X = snps[start_i: stop_i]
            D = sp.dot(X, X.T) / n
            ref_ld_matrices.append(D)
        ret_dict['ref_ld_matrices']=ref_ld_matrices
    return ret_dict


def annopred_inf(beta_hats, pr_sigi, n=1000, reference_ld_mats=None, ld_window_size=100):
    """
    infinitesimal model with snp-specific heritability derived from annotation
    used as the initial values for MCMC of non-infinitesimal model
    """
    num_betas = len(beta_hats)
    updated_betas = sp.empty(num_betas)
    m = len(beta_hats)

    for i, wi in enumerate(range(0, num_betas, ld_window_size)):
        start_i = wi
        stop_i = min(num_betas, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        Li = 1.0/pr_sigi[start_i: stop_i]
        D = reference_ld_mats[i]
        A = (n/(1))*D + sp.diag(Li)
        A_inv = linalg.pinv(A)
        updated_betas[start_i: stop_i] = sp.dot(A_inv / (1.0/n) , beta_hats[start_i: stop_i])  # Adjust the beta_hats

    return updated_betas


def annopred_genomewide(data_file=None, ld_radius = None, ld_dict=None, out_file_prefix=None, ps=None, 
               n=None, h2=None, num_iter=None, zero_jump_prob=0.05, burn_in=5, PRF=None):
    """
    Calculate LDpred for a genome
    """    
    prf_chr = PRF['chrom']
    prf_sids = PRF['sids']
    prf_pi = PRF['pi']
    prf_sigi2 = PRF['sigi2']

    df = h5py.File(data_file,'r')
    has_phenotypes=False
    if 'y' in df.keys():
        'Validation phenotypes found.'
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        has_phenotypes=True
        risk_scores_pval_derived = sp.zeros(num_individs)
        risk_scores_pval_derived_inf = sp.zeros(num_individs)

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ld_dict = ld_dict['chrom_ld_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        
    print('LD radius used: %d' % ld_radius)
    results_dict = {}
    num_snps = 0
    sum_beta2s = 0
    cord_data_g = df['cord_data']

    for chrom_str in chromosomes_list: 
        if chrom_str in cord_data_g.keys():
            g = cord_data_g[chrom_str]
            betas = g['betas'][...]
            n_snps = len(betas)
            num_snps += n_snps
            sum_beta2s += sp.sum(betas ** 2)
        
    L = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda = sp.mean(n * sum_beta2s / float(num_snps))
#    print('Genome-wide lambda inflation:', chi_square_lambda,)
    print('Genome-wide mean LD score:', L)
    gw_h2_ld_score_est = max(0.0001, (max(1, chi_square_lambda) - 1) / (n * (L / num_snps)))
    print('Estimated genome-wide heritability:', gw_h2_ld_score_est)
    
    #assert chi_square_lambda>1, 'Check the summary statistic file'
    if h2 is None:
        h2 = gw_h2_ld_score_est
    print(h2)
    h2_new = sp.sum(prf_sigi2)
    sig_12 = (1.0)/n     #######################
    pr_sig = {}
    pr_p = {}
    annopred_inf_chrom_dict = {}
    print('Calculating initial values for MCMC using infinitesimal model')
    for chrom_str in chromosomes_list:
        if chrom_str in cord_data_g.keys():
            print('Calculating posterior betas for Chromosome %s'%((chrom_str.split('_'))[1]))           
            g = cord_data_g[chrom_str]

            #Filter monomorphic SNPs
            snp_stds = g['snp_stds_ref'][...]
            snp_stds = snp_stds.flatten()
            ok_snps_filter = snp_stds>0
            pval_derived_betas = g['betas'][...]
            pval_derived_betas = pval_derived_betas[ok_snps_filter]
            sids = g['sids'][...]
            sids = sids[ok_snps_filter]
            chri = int(chrom_str.split('_')[1])
            prf_sids_chri = prf_sids[prf_chr==chri]
            prf_pi_chri = prf_pi[prf_chr==chri]
            prf_sigi2_chri = prf_sigi2[prf_chr==chri]
            if len(prf_sids_chri)==len(sids):
                if sum(prf_sids_chri==sids)==len(prf_sids_chri):
                    pr_p[chrom_str] = sp.copy(prf_pi_chri)
                    pr_sig[chrom_str] = sp.copy(prf_sigi2_chri)
                else:
                    print('Order of SNPs does not match, sorting prior files')
                    pr_p[chrom_str] = sp.zeros(len(sids))
                    pr_sig[chrom_str] = sp.zeros(len(sids))
                    for i, sid in enumerate(sids):
                        pr_p[chrom_str][i] = prf_pi_chri[prf_sids_chri==sid]
                        pr_sig[chrom_str][i] = prf_sigi2_chri[prf_sids_chri==sid]
            else:
                print('More SNPs found in prior file, extracting SNPs from prior files')
                pr_p[chrom_str] = sp.zeros(len(sids))
                pr_sig[chrom_str] = sp.zeros(len(sids))
                for i, sid in enumerate(sids):
                    pr_p[chrom_str][i] = prf_pi_chri[prf_sids_chri==sid]
                    pr_sig[chrom_str][i] = prf_sigi2_chri[prf_sids_chri==sid]
            pr_sig[chrom_str] = h2*pr_sig[chrom_str]/h2_new
            if h2 is not None:
                h2_chrom = sp.sum(pr_sig[chrom_str])            
            else:
                h2_chrom = gw_h2_ld_score_est * (n_snps / float(num_snps))
            start_betas = annopred_inf(pval_derived_betas, pr_sigi=pr_sig[chrom_str], reference_ld_mats=chrom_ref_ld_mats[chrom_str], n=n, ld_window_size=2*ld_radius)
            annopred_inf_chrom_dict[chrom_str]=start_betas
    
    
    for p in ps:
        print('Starting AnnoPred with ', p)
        p_str = p
        results_dict[p_str]={}
        
        if out_file_prefix:
            #Preparing output files
            raw_effect_sizes = []
            annopred_effect_sizes = []
            annopred_inf_effect_sizes = []
            out_sids = []
            chromosomes = []
            out_positions = []
            out_nts = []
    
        out = []
        out_inf = []
        out.append('The input prior p is '+str(prf_pi[0])+'\n')
        out.append('Estimated Genome-wide heritability: '+str(gw_h2_ld_score_est)+'\n')
        out.append('Posterior variance for each snp: '+str(sig_12)+'\n')
        print('Estimated Genome-wide heritability from Priors:', h2)
        print('Posterior variance for each snp:', sig_12)
        for chrom_str in chromosomes_list:
            if chrom_str in cord_data_g.keys():
                g = cord_data_g[chrom_str]
                if has_phenotypes:
                    if 'raw_snps_val' in g.keys():
                        raw_snps = g['raw_snps_val'][...]
                    else:
                        raw_snps = g['raw_snps_ref'][...]
                    
                #Filter monomorphic SNPs
                snp_stds = g['snp_stds_ref'][...]
                snp_stds = snp_stds.flatten()
                ok_snps_filter = snp_stds>0
                snp_stds = snp_stds[ok_snps_filter]
                pval_derived_betas = g['betas'][...]
                pval_derived_betas = pval_derived_betas[ok_snps_filter]
                positions = g['positions'][...]
                positions = positions[ok_snps_filter]
                sids = g['sids'][...]
                sids = sids[ok_snps_filter]
                log_odds = g['log_odds'][...]
                log_odds = log_odds[ok_snps_filter]
                nts = g['nts'][...]
                nts = nts[ok_snps_filter]
    
                prf_pi_chri_sorted = pr_p[chrom_str]
                prf_sigi2_chri_sorted = pr_sig[chrom_str]
                
                if out_file_prefix:
                    chromosomes.extend([chrom_str]*len(pval_derived_betas))
                    out_positions.extend(positions)
                    out_sids.extend(sids)
                    raw_effect_sizes.extend(log_odds)
                    out_nts.extend(nts)
                
                n_snps = len(pval_derived_betas)
                
                if h2 is not None:
                    h2_chrom = sp.sum(prf_sigi2_chri_sorted)
                    #h2_chrom = h2 * (n_snps / float(num_snps))            
                else:
                    h2_chrom = gw_h2_ld_score_est * (n_snps / float(num_snps))
                #print('Prior parameters: p=%0.3f, n=%d, m=%d, h2_chrom=%0.4f' % (p, n, n_snps, h2_chrom))
                res_dict = non_infinitesimal_mcmc(pval_derived_betas, Pi = prf_pi_chri_sorted, Sigi2=prf_sigi2_chri_sorted, sig_12=sig_12, h2=h2_chrom, n=n, ld_radius=ld_radius,
                                        num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                        start_betas=annopred_inf_chrom_dict[chrom_str], zero_jump_prob=zero_jump_prob)            
                updated_betas = res_dict['betas']
                updated_inf_betas = res_dict['inf_betas']
                sum_sqr_effects = sp.sum(updated_betas ** 2)
                if sum_sqr_effects>gw_h2_ld_score_est:
                    print('Sum of squared updated effects estimates seems too large:', sum_sqr_effects)
                    print('This suggests that the Gibbs sampler did not convergence.')
                
                print('Calculating scores for Chromosome %s'%((chrom_str.split('_'))[1]))
                updated_betas = updated_betas / (snp_stds.flatten())
                updated_inf_betas = updated_inf_betas / (snp_stds.flatten())
                annopred_effect_sizes.extend(updated_betas)
                annopred_inf_effect_sizes.extend(updated_inf_betas)
                #if has_phenotypes:
                if has_phenotypes:
                    prs = sp.dot(updated_betas, raw_snps)
                    prs_inf = sp.dot(updated_inf_betas, raw_snps)
                    risk_scores_pval_derived += prs
                    risk_scores_pval_derived_inf += prs_inf
                    corr = sp.corrcoef(y, prs)[0, 1]
                    r2 = corr ** 2
                    corr_inf = sp.corrcoef(y, prs_inf)[0, 1]
                    r2_inf = corr_inf ** 2
#                    print('The R2 prediction accuracy of PRS using %s was: %0.4f' %(chrom_str, r2))
#                    print('The R2 prediction accuracy of PRS using %s was: %0.4f' %(chrom_str, r2_inf))
                    out.append('The R2 prediction accuracy of PRS using '+chrom_str+' was '+str(r2)+'\n')
                    out_inf.append('The R2 prediction accuracy of PRS using '+chrom_str+' was '+str(r2_inf)+'\n')
        
                    
#        print('There were %d (SNP) effects' % num_snps)
        if has_phenotypes:
            num_indivs = len(y)
            results_dict[p_str]['y']=y
            results_dict[p_str]['risk_scores_pd']=risk_scores_pval_derived
#            print('Prediction accuracy was assessed using %d individuals.'%(num_indivs))
            out.append('Prediction accuracy was assessed using '+str(num_indivs)+' individuals\n')
            
            corr = sp.corrcoef(y, risk_scores_pval_derived)[0, 1]
            r2 = corr ** 2
            results_dict[p_str]['r2_pd']=r2
#            print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1-r2)**2)/num_indivs))
            out.append('The  R2 prediction accuracy (observed scale) for the whole genome was: '+str(r2)+' ('+str(((1-r2)**2)/num_indivs)+')\n')
            
            corr_inf = sp.corrcoef(y, risk_scores_pval_derived_inf)[0, 1]
            r2_inf = corr_inf ** 2
            results_dict[p_str]['r2_pd']=r2_inf
#            print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2_inf, ((1-r2_inf)**2)/num_indivs))
            out_inf.append('The  R2 prediction accuracy (observed scale) for the whole genome was: '+str(r2_inf)+' ('+str(((1-r2_inf)**2)/num_indivs)+')\n')
            
            if corr<0:
                risk_scores_pval_derived = -1* risk_scores_pval_derived
            auc = pred_accuracy(y,risk_scores_pval_derived)
            print('AnnoPred AUC/COR for the whole genome was: %0.4f'%auc)
            out.append('AUC/COR for the whole genome was: '+str(auc)+'\n')
    
            if corr_inf<0:
                risk_scores_pval_derived_inf = -1* risk_scores_pval_derived_inf
            auc_inf = pred_accuracy(y,risk_scores_pval_derived_inf)
            print('AnnoPred-inf AUC/COR for the whole genome was: %0.4f'%auc_inf)
            out_inf.append('AUC/COR for the whole genome was: '+str(auc_inf)+'\n')
    
            sp.savetxt('%s_y_'%(out_file_prefix)+str(p)+'.txt',y)
            #sp.savetxt('%s_prs_'%(out_file_prefix)+str(p)+'.txt',risk_scores_pval_derived)
            #sp.savetxt('%s_prs-inf'%(out_file_prefix)+str(p)+'.txt',risk_scores_pval_derived_inf)
    
            #Now calibration                               
            denominator = sp.dot(risk_scores_pval_derived.T, risk_scores_pval_derived)
            y_norm = (y-sp.mean(y))/sp.std(y)
            numerator = sp.dot(risk_scores_pval_derived.T, y_norm)
            regression_slope = (numerator / denominator)#[0][0]
#            print('The slope for predictions with P-value derived  effects is:',regression_slope)
            out.append('The slope for predictions with P-value derived  effects is: '+str(regression_slope)+'\n')
            results_dict[p_str]['slope_pd']=regression_slope
    
            ff = open('%s_non_inf_auc_'%(out_file_prefix)+str(p)+'.txt',"w")
            ff.writelines(out)
            ff.close()
            
            ff_inf = open('%s_inf_auc_'%(out_file_prefix)+str(p)+'.txt',"w")
            ff_inf.writelines(out_inf)
            ff_inf.close()     
            sp.savetxt('%s_prs_'%(out_file_prefix)+str(p)+'.txt',risk_scores_pval_derived)
            sp.savetxt('%s_prs-inf'%(out_file_prefix)+str(p)+'.txt',risk_scores_pval_derived_inf)
        
        weights_out_file = '%s_non_inf_betas_'%(out_file_prefix)+str(p)+'.txt' ###################################
        with open(weights_out_file,'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta     AnnoPred_beta\n')
            for chrom, pos, sid, nt, raw_beta, annopred_beta in it.izip(chromosomes, out_positions, out_sids, out_nts, raw_effect_sizes, annopred_effect_sizes):
                nt1,nt2 = nt[0],nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n'%(chrom, pos, sid, nt1, nt2, raw_beta, annopred_beta))
    
        weights_out_file = '%s_inf_betas_'%(out_file_prefix)+str(p)+'.txt'
        with open(weights_out_file,'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta    AnnoPred_inf_beta \n')
            for chrom, pos, sid, nt, raw_beta, annopred_inf_beta in it.izip(chromosomes, out_positions, out_sids, out_nts, raw_effect_sizes, annopred_inf_effect_sizes):
                nt1,nt2 = nt[0],nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n'%(chrom, pos, sid, nt1, nt2, raw_beta, annopred_inf_beta))


        
def non_infinitesimal_mcmc(beta_hats, Pi, Sigi2, sig_12, start_betas=None, h2=None, n=1000, ld_radius=100, num_iter=60, burn_in=10, zero_jump_prob=0.05, ld_dict=None):
    """
    MCMC of non-infinitesimal model
    """
    m = len(beta_hats)
    
    curr_betas = sp.copy(start_betas)
    curr_post_means = sp.zeros(m)
    avg_betas = sp.zeros(m)

    # Iterating over effect estimates in sequential order
    iter_order = sp.arange(m)
    
    for k in range(num_iter):  #Big iteration

        #Force an alpha shrink if estimates are way off compared to heritability estimates.  (Improves MCMC convergence.)
        h2_est = max(0.00001,sp.sum(curr_betas ** 2))
        alpha = min(1-zero_jump_prob, 1.0 / h2_est, (h2 + 1 / sp.sqrt(n)) / h2_est)
        rand_ps = sp.random.random(m)

        for i, snp_i in enumerate(iter_order):
            if Sigi2[snp_i]==0:
                curr_post_means[snp_i] = 0
                curr_betas[snp_i] = 0
            else:
                hdmp = (Sigi2[snp_i]/Pi[snp_i])#(h2 / Mp)
                hdmpn = hdmp + sig_12#1.0 / n
                hdmp_hdmpn = (hdmp / hdmpn)
                c_const = (Pi[snp_i] / sp.sqrt(hdmpn))
                d_const = (1 - Pi[snp_i]) / (sp.sqrt(sig_12))
    
                start_i = max(0, snp_i - ld_radius)
                focal_i = min(ld_radius, snp_i)
                stop_i = min(m, snp_i + ld_radius + 1)
                
                #Local LD matrix
                D_i = ld_dict[snp_i]
                
                #Local (most recently updated) effect estimates
                local_betas = curr_betas[start_i: stop_i]
                
                #Calculate the local posterior mean, used when sampling.
                local_betas[focal_i] = 0
                res_beta_hat_i = beta_hats[snp_i] - sp.dot(D_i , local_betas)
                b2 = res_beta_hat_i ** 2
        
                d_const_b2_exp = d_const * sp.exp(-b2 / (2.0*sig_12))
                if sp.isreal(d_const_b2_exp):
                    numerator = c_const * sp.exp(-b2 / (2.0 * hdmpn))
                    if sp.isreal(numerator):
                        if numerator == 0:
                            postp = 0
                        else:
                            postp = numerator / (numerator + d_const_b2_exp)
                            assert sp.isreal(postp), 'Posterior mean is not a real number?' 
                    else:
                        postp = 0
                else:
                    postp = 1
                curr_post_means[snp_i] = hdmp_hdmpn * postp * res_beta_hat_i
        
                if rand_ps[i] < postp * alpha:
                    #Sample from the posterior Gaussian dist.
                    proposed_beta = stats.norm.rvs(0, (hdmp_hdmpn) * sig_12, size=1) + hdmp_hdmpn * res_beta_hat_i
        
                else:
                    #Sample 0
                    proposed_beta = 0
        
                curr_betas[snp_i] = proposed_beta  #UPDATE BETA

        if k >= burn_in:
            avg_betas += curr_post_means #Averaging over the posterior means instead of samples.

    avg_betas = avg_betas/float(num_iter-burn_in)

    return {'betas':avg_betas, 'inf_betas':start_betas}


"""
p_dict = {'coord':None, 'ld_radius':None, 'local_ld_file_prefix':None, 'hfile':None, 'pfile':None, 'PS':None, 'out':None,
          'N':None, 'num_iter': 60, 'H2':None, 'user_h2':None}
"""

def main(p_dict):
    local_ld_dict_file = '%s_ldradius%d.pickled.gz'%(p_dict['local_ld_file_prefix'], p_dict['ld_radius'])
    
    p_dict['PS'] = [p_dict['PS']]
    if not os.path.isfile(local_ld_dict_file):
        df = h5py.File(p_dict['coord'])
                 
        chrom_ld_scores_dict = {}
        chrom_ld_dict = {}
        chrom_ref_ld_mats = {}
        ld_score_sum = 0
        num_snps = 0
        print('Calculating LD information w. radius %d'% p_dict['ld_radius'])

        cord_data_g = df['cord_data']

        for chrom_str in cord_data_g.keys():
            print('Working on %s'%chrom_str)
            g = cord_data_g[chrom_str]
            if 'raw_snps_ref' in g.keys():
                raw_snps = g['raw_snps_ref'][...]
                snp_stds = g['snp_stds_ref'][...]
                snp_means = g['snp_means_ref'][...]
            
            
            #Filter monomorphic SNPs
            ok_snps_filter = snp_stds>0
            ok_snps_filter = ok_snps_filter.flatten()
            raw_snps = raw_snps[ok_snps_filter]
            snp_means = snp_means[ok_snps_filter]
            snp_stds = snp_stds[ok_snps_filter]

            n_snps = len(raw_snps)
            snp_means.shape = (n_snps,1)   
            snp_stds.shape = (n_snps,1)   
            
            
            # Normalize SNPs..
            snps = sp.array((raw_snps - snp_means)/snp_stds,dtype='float32')
            assert snps.shape==raw_snps.shape, 'Array Shape mismatch'
            ret_dict = get_LDpred_ld_tables(snps, ld_radius=p_dict['ld_radius'], ld_window_size=2*p_dict['ld_radius'])
            chrom_ld_dict[chrom_str] = ret_dict['ld_dict']
            chrom_ref_ld_mats[chrom_str] = ret_dict['ref_ld_matrices']
            ld_scores = ret_dict['ld_scores']
            chrom_ld_scores_dict[chrom_str] = {'ld_scores':ld_scores, 'avg_ld_score':sp.mean(ld_scores)}
            ld_score_sum += sp.sum(ld_scores)
            num_snps += n_snps
        avg_gw_ld_score = ld_score_sum / float(num_snps)
        ld_scores_dict = {'avg_gw_ld_score': avg_gw_ld_score, 'chrom_dict':chrom_ld_scores_dict}    
        
        print('Done calculating the LD table and LD score, writing to file:', local_ld_dict_file)
        print('Genome-wide average LD score was:', ld_scores_dict['avg_gw_ld_score'])
        ld_dict = {'ld_scores_dict':ld_scores_dict, 'chrom_ld_dict':chrom_ld_dict, 'chrom_ref_ld_mats':chrom_ref_ld_mats}
        f = gzip.open(local_ld_dict_file, 'wb')
        cPickle.dump(ld_dict, f, protocol=2)
        f.close()
        print('LD information is now pickled.')
    else:
        print('Loading LD information from file: %s'%local_ld_dict_file)
        f = gzip.open(local_ld_dict_file, 'r')
        ld_dict = cPickle.load(f)
        f.close()

    if p_dict['user_h2'] is not None:
        print('Starting calculation using user provided h2 files as priors')
        print('Loading prior information from file: %s'%p_dict['user_h2'])
        with open(p_dict['user_h2']) as f:
            data = f.readlines()
        prf_chr = sp.empty(len(data),dtype='int8')
        prf_sids = []
        prf_pi = sp.zeros(len(data))
        prf_sigi2 = sp.zeros(len(data))
        for i,line in enumerate(data):
            li = line.split()
            prf_chr[i] = int(li[0])
            prf_sids.append(li[1]) 
            prf_pi[i] = p_dict['PS'][0]         
            prf_sigi2[i] = float(li[2]) 
        print('The input prior p is: ', p_dict['PS'])
        prf_sids = sp.array(prf_sids,dtype='str')
        prf = {}
        prf['chrom'] = prf_chr
        prf['sids'] = prf_sids
        prf['pi'] = prf_pi
        prf['sigi2'] = prf_sigi2
        out_user_h2_prefix = p_dict['out']+'_user_h2'
        annopred_genomewide(data_file=p_dict['coord'], out_file_prefix=out_user_h2_prefix, ps=p_dict['PS'], ld_radius=p_dict['ld_radius'], 
                          ld_dict = ld_dict, n=p_dict['N'], num_iter=p_dict['num_iter'], h2=p_dict['H2'], PRF = prf)
    else:
        #### no user-provided heritability files provided ####
        ##################### using hfile as prior #######################
        print('Starting calculation using h2 files as priors')
        print('Loading prior information from file: %s'%p_dict['hfile'])
        with open(p_dict['hfile']) as f:
            data = f.readlines()
        prf_chr = sp.empty(len(data),dtype='int8')
        prf_sids = []
        prf_pi = sp.zeros(len(data))
        prf_sigi2 = sp.zeros(len(data))
        for i,line in enumerate(data):
            li = line.split()
            prf_chr[i] = int(li[0])
            prf_sids.append(li[1]) 
            prf_pi[i] = p_dict['PS'][0]         
            prf_sigi2[i] = float(li[2]) 
        print('The input prior p is: ', p_dict['PS'])
        prf_sids = sp.array(prf_sids,dtype='str')
        prf = {}
        prf['chrom'] = prf_chr
        prf['sids'] = prf_sids
        prf['pi'] = prf_pi
        prf['sigi2'] = prf_sigi2
        out_h2_prefix = p_dict['out']+'_h2'
        #H2 = sp.sum(prf_sigi2)
        annopred_genomewide(data_file=p_dict['coord'], out_file_prefix=out_h2_prefix, ps=p_dict['PS'], ld_radius=p_dict['ld_radius'], 
                              ld_dict = ld_dict, n=p_dict['N'], num_iter=p_dict['num_iter'], PRF = prf)    
               
        ##################### using pfile as prior #######################
        if p_dict['pfile'] is not None:
            print('Starting calculation using p_T files as priors')
            print('Loading prior information from file: %s'%p_dict['pfile'])
            with open(p_dict['pfile']) as f:
                data = f.readlines()
            prf_chr = sp.empty(len(data),dtype='int8')
            prf_sids = []
            prf_pi = sp.zeros(len(data))
            prf_sigi2 = sp.zeros(len(data))
            for i,line in enumerate(data):
                li = line.split()
                prf_chr[i] = int(li[0])
                prf_sids.append(li[1]) 
                prf_pi[i] = float(li[2])         
                prf_sigi2[i] = float(li[3]) 
            print('The input prior p is: ', p_dict['PS'])
            prf_sids = sp.array(prf_sids,dtype='str')
            prf = {}
            prf['chrom'] = prf_chr
            prf['sids'] = prf_sids
            prf['pi'] = prf_pi
            prf['sigi2'] = prf_sigi2
            #H2 = sp.sum(prf_sigi2)
            out_pT_prefix = p_dict['out']+'_pT'
            annopred_genomewide(data_file=p_dict['coord'], out_file_prefix=out_pT_prefix, ps=p_dict['PS'], ld_radius=p_dict['ld_radius'], 
                            ld_dict = ld_dict, n=p_dict['N'], num_iter=p_dict['num_iter'], PRF = prf)
    
    
