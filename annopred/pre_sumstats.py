import numpy as np
import h5py

def get_1000G_snps(sumstats, out_file):
    sf = np.loadtxt(sumstats,dtype=str,skiprows=1)
    h5f = h5py.File('ref/AnnotMatrix/1000G_SNP_info.h5','r')
    rf = h5f['snp_chr'][:]
    h5f.close()
    ind1 = np.in1d(sf[:,1],rf[:,2])
    ind2 = np.in1d(rf[:,2],sf[:,1])
    sf1 = sf[ind1]
    rf1 = rf[ind2]
    ### check order ###
    if sum(sf1[:,1]==rf1[:,2])==len(rf1[:,2]):
        print('Good!')
    else:
        print('Shit happens, sorting sf1 to have the same order as rf1')
        O1 = np.argsort(sf1[:,1])
        O2 = np.argsort(rf1[:,2])
        O3 = np.argsort(O2)
        sf1 = sf1[O1][O3]
    out = ['hg19chrc snpid a1 a2 bp or p'+'\n']
    for i in range(len(sf1[:,1])):
        out.append(sf1[:,0][i]+' '+sf1[:,1][i]+' '+sf1[:,2][i]+' '+sf1[:,3][i]+' '+rf1[:,1][i]+' '+sf1[:,5][i]+' '+sf1[:,6][i]+'\n')
    ff = open(out_file,"w")
    ff.writelines(out)
    ff.close()


