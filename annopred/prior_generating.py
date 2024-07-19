#!/usr/bin/env python

### generate prior file from h5py file directly ###
### generate_h2_pT generates two prior files from the results of LDSC and a fixed annotation file ###
### generate_h2_from_user generates one prior file from the user provided prior file ###
import h5py
import numpy as np
import os
from collections import Counter
from collections import defaultdict
import datetime
import math

def generate_h2_pT(h5py_file, LDSC_results_file, output_h2, PS, output_pT, annotation_flag):
    # generate two types of prior files
    ### load the fixed input file ###
    h5f1 = h5py.File('ref/AnnotMatrix/baseline.h5','r')
    baseline = h5f1['annot'][:]
    h5f1.close()

    if annotation_flag=='tier0':
        h5f1 = h5py.File('ref/AnnotMatrix/tier0.h5','r')
        tier = h5f1['annot'][:]
        h5f1.close()
    elif annotation_flag=='tier1':
        h5f1 = h5py.File('ref/AnnotMatrix/tier1.h5','r')
        tier = h5f1['annot'][:]
        h5f1.close()
    elif annotation_flag=='tier2':
        h5f1 = h5py.File('ref/AnnotMatrix/tier2.h5','r')
        tier = h5f1['annot'][:]
        h5f1.close()
    elif annotation_flag=='tier3':
        h5f1 = h5py.File('ref/AnnotMatrix/tier3.h5','r')
        tier = h5f1['annot'][:]
        h5f1.close()
    else:
        exit("Illegal tier name!")
    
    annot = np.concatenate((baseline,tier),axis=1)
    h5f2 = h5py.File('ref/AnnotMatrix/1000G_SNP_info.h5','r')
    snp_chr = h5f2['snp_chr'][:]
    h5f2.close()
    ### get the snp list from h5py ###
    chromosomes_list = ['chrom_%d'%(x) for x in range(1,23)]
    chromosomes_list.append('chrom_X')
    
    df = h5py.File(h5py_file,'r')
    cord_data_g = df['cord_data']
    
    SNPids = []
    for chrom_str in chromosomes_list:
        if chrom_str in cord_data_g.keys():
            g = cord_data_g[chrom_str]
            #Filter monomorphic SNPs (SNPs with variance equal to 0)
            snp_stds = g['snp_stds_ref'][...]
            snp_stds = snp_stds.flatten()
            ok_snps_filter = snp_stds>0
            pval_derived_betas = g['betas'][...]
            pval_derived_betas = pval_derived_betas[ok_snps_filter]
            sids = g['sids'][...]
            SNPids = np.append(SNPids,sids[ok_snps_filter])
    num_snps = len(SNPids)
    ### overlap with SNP in annot files ###
    stt1 = np.in1d(snp_chr[:,2],SNPids)
    ant1 = annot[stt1]
    snp_chr1 = snp_chr[stt1]
    ### check order ###
    if sum(snp_chr1[:,2]==SNPids)==len(SNPids):
        print('Good!')
    else:
        print('Shit happens, sorting ant1 to have the same order as SNPids')
        O1 = np.argsort(snp_chr1[:,2])
        O2 = np.argsort(SNPids)
        O3 = np.argsort(O2)
        ant1 = ant1[O1][O3]

    ### load LDSC results ###
    LD_results = np.genfromtxt(LDSC_results_file,dtype=None,names=True)
    
    tau0 = LD_results['Coefficient']
  
    ### get heritability  ###
    sig2_0 = np.dot(ant1,tau0)
    
    ### adjust for minus terms ###
    sig2_0[sig2_0<0] = np.repeat(min(sig2_0[sig2_0>0]),np.sum(sig2_0<0))
    np.sum(sig2_0)
    
    ### save prior file (h2) ###
    h2_out = []
    for i in range(len(sig2_0)):
        h2_out.append(str(snp_chr1[:,0][i])+' '+str(snp_chr1[:,2][i])+' '+str(sig2_0[i])+'\n')
    #np.savetxt(output_h2,(snp_chr1[:,0],snp_chr1[:,1],sig2_0),fmt="%s")
    ff = open(output_h2,"w")
    ff.writelines(h2_out)
    ff.close()

    ### start calculating p_T ###
    M = np.empty(annot.shape[1])
    for i in range(len(M)):
        M[i] = np.sum(np.logical_and(annot[:,0],annot[:,i]))
    bgt = datetime.datetime.now()
    M_T = defaultdict(int)
    for i in range(annot.shape[0]):
        tup_i = tuple(annot[i])
        M_T[tup_i] += 1
    edt = datetime.datetime.now()
    print(edt-bgt)
    bgt = datetime.datetime.now()
    N_T = defaultdict(int)
    for i in range(ant1.shape[0]):
        tup_i = tuple(ant1[i])
        N_T[tup_i] += 1
    edt = datetime.datetime.now()
    print(edt-bgt)


    H0 = np.dot(M,tau0)
    N0 = float(len(SNPids))
    sig2V = np.dot(ant1,tau0)

    # N_T = {x:annotV1.count(x) for x in annotV1}
    
    M_TV = np.empty(ant1.shape[0])
    N_TV = np.empty(ant1.shape[0])
    for i in range(ant1.shape[0]):
        tup_i = tuple(ant1[i])
        M_TV[i] = M_T[tup_i]
        N_TV[i] = N_T[tup_i]

#    for ps in PS:
    pr_p = (PS*N0/H0)*M_TV*sig2V/N_TV
    sig2 = M_TV*sig2V/N_TV
    m1 = min(pr_p[pr_p>0])
    m2 = min(sig2[sig2>0])
    pr_p[pr_p<0] = np.repeat(m1,np.sum(pr_p<0))
    sig2[sig2<0] = np.repeat(m2,np.sum(sig2<0))
    pr_p[pr_p>1] = np.repeat(1,np.sum(pr_p>1))
    pT_out = []
    for i in range(len(sig2)):
        pT_out.append(str(snp_chr1[:,0][i])+' '+str(snp_chr1[:,2][i])+' '+str(pr_p[i])+' '+str(sig2[i])+'\n')
#    ff = open(output_pT+'_'+str(ps)+'_file.txt',"w")
    ff = open(output_pT,"w")
    ff.writelines(pT_out)
    ff.close()

    return math.ceil(num_snps/3000.0)



def generate_h2_from_user(user_provided_h2, h5py_file, output):
    ### user_provided_h2 format: Chr, SNP_id, heritability ###
    ### load the fixed input file ###
    ## Note: gonna take huge memory!!! Probably need to optimize this part, for example, read in .gz files directly ##
    user_h2 = np.loadtxt(user_provided_h2,dtype=str,skiprows=0)
    ### get the snp list from h5py ###
    chromosomes_list = ['chrom_%d'%(x) for x in range(1,23)]
    chromosomes_list.append('chrom_X')
    
    df = h5py.File(h5py_file,'r')
    cord_data_g = df['cord_data']
    
    SNPids = []
    for chrom_str in chromosomes_list:
        if chrom_str in cord_data_g.keys():
            g = cord_data_g[chrom_str]
            #Filter monomorphic SNPs (SNPs with variance equal to 0)
            snp_stds = g['snp_stds_ref'][...]
            snp_stds = snp_stds.flatten()
            ok_snps_filter = snp_stds>0
            pval_derived_betas = g['betas'][...]
            pval_derived_betas = pval_derived_betas[ok_snps_filter]
            sids = g['sids'][...]
            SNPids = np.append(SNPids,sids[ok_snps_filter])
    num_snps = len(SNPids)
    ### overlap with SNP in annot files ###
    stt1 = np.in1d(user_h2[:,1],SNPids)
    user_h2 = user_h2[stt1]
    
    ### check order ###
    if sum(user_h2[:,1]==SNPids)==len(SNPids):
        print('Good!')
    else:
        print('Shit happens, sorting user_h2 to have the same order as SNPids')
        O1 = np.argsort(user_h2[:,1])
        O2 = np.argsort(SNPids)
        O3 = np.argsort(O2)
        user_h2 = user_h2[O1][O3]

    ### save prior file (h2) ###
    h2_out = []
    for i in range(len(user_h2[:,0])):
        h2_out.append(user_h2[:,0][i]+' '+user_h2[:,1][i]+' '+user_h2[:,2][i]+'\n')
#    np.savetxt(output,(user_h2[:,0],user_h2[:,1],user_h2[:,2]),fmt="%s")
    ff = open(output,"w")
    ff.writelines(h2_out)
    ff.close()
    return [np.sum(user_h2[:,2].astype(np.float)), math.ceil(num_snps/3000.0)]
