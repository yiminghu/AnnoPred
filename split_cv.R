args = commandArgs(trailingOnly=TRUE)
orig_gt_file = args[1] ## original genotype file (plink bed/bim/fam)
output_gt_fam_file1 = args[2] ## list of the first half of individuals
output_gt_fam_file2 = args[3] ## list of the second half of individuals

options(stringsAsFactors=F)
fam = read.table(orig_gt_file, header=F)
N = nrow(fam)
pt1 = sample(1:N, floor(N/2))
pt2 = (1:N)[-pt1]
write.table(fam[pt1,1:2], output_gt_fam_file1, quote=F, row.names=F, col.names=F, sep='\t')
write.table(fam[pt2,1:2], output_gt_fam_file2, quote=F, row.names=F, col.names=F, sep='\t')
