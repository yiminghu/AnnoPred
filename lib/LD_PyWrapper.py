#!/usr/bin/env python

from __future__ import print_function
import subprocess,sys,os,itertools


def formatOptions(optsList):
    command = []
    for f,a in optsList:
        command.append(f)
        if a: command.append(a)
    return(command)

def loadLDPath():
    configPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + "/LDSC.config"
    if not os.path.isfile(configPath):
        raise Exception("Please provide LDSC.config file in the AnnoPred directory")
    with open(configPath, "r") as fp:
        option = fp.readline().split()
        ldPath = option[1]
        return ldPath

#Which of these options should be user supplied?
def callMunge(sumstats, n_case, n_ctrl, ldPath, refPath):
    print("Calling munge_sumstats.py...")
    mungeOut = "/SS_Formatted/"
    #outName = os.path.basename(sumstats).split("_")[1]
    mungeFlags = ["--" + flag for flag in ["N-cas", "N-con", "merge-alleles", "sumstats", "out"]]
    mungeArgs = [n_case, n_ctrl, refPath + "Misc/w_hm3.snplist", sumstats, "Curated_GWAS"]
    mungeOptsList = [(f,a) for f,a in zip(mungeFlags, mungeArgs)]
    mungeOpts = formatOptions(mungeOptsList)
    subprocess.call(["python", ldPath + "/munge_sumstats.py"] + mungeOpts)
    return mungeArgs[4]


def callLDSC(sumstats, n_case, n_ctrl):
    ldPath = loadLDPath()
    refPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + "/ref/"
    GenoAnnotations = ["Brain", "GI", "Lung", "Heart", "Blood", "Muscle", "Epithelial"]
    AnnotationPaths = ["/Annotations/Baseline/baseline.", "/Annotations/GenoCanyon/GenoCanyon_Func."] + ["/Annotations/GenoSkyline/" + g + "." for g in GenoAnnotations] 
    mungeFile = callMunge(sumstats, str(n_case), str(n_ctrl), ldPath, refPath)
    refFiles = [refPath + a for a in AnnotationPaths]
    print("Running LD Score calculation...")
    ldscFlags = ["--" + flag for flag in ["h2", "ref-ld-chr", "w-ld-chr", "frqfile-chr", "overlap-annot", "print-coefficients", "out"]]
    ldscArgs = [mungeFile + '.sumstats.gz', ','.join(refFiles), refPath+'/Misc/weights.', refPath+'/Misc/1000G.mac5eur.', '', '', 'SNP_Heritability']
    ldscOptsList = [(f,a) for f,a in zip(ldscFlags, ldscArgs)]
    ldscOpts = formatOptions(ldscOptsList)
    subprocess.call(["python", ldPath + "/ldsc.py"] + ldscOpts)
    print(subprocess.list2cmdline(["python", ldPath + "/ldsc.py"] + ldscOpts))

if __name__ == "__main__":
    callLDSC("/net/zhao/ql68/GWAS/PD/Data/Data_SimonSanchez_Lite.txt", 1178, 5213)
    #callLDSC(sys.argv[1], sys.argv[2], sys.argv[3])


