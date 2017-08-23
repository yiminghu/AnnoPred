args = commandArgs(trailingOnly=TRUE)
cv1_output_prefix = args[1] ## output file prefix of cv1
cv2_output_prefix = args[2] ## output file prefix of cv2
cv_range = args[3]

get_res <- function(strV){
	N = length(strV)
	as.numeric(unlist(strsplit(strV[N-1], ": "))[2])
}
pv = unlist(strsplit(cv_range, ' '))

res_cv1 = res_cv2 = matrix(0,2,length(pv))
for(i in 1:length(pv)){
	h2_cv1 = readLines(paste0(cv1_output_prefix, "_h2_non_inf_auc_", pv[i], ".txt"))
	pT_cv1 = readLines(paste0(cv1_output_prefix, "_pT_non_inf_auc_", pv[i], ".txt"))
	h2_cv2 = readLines(paste0(cv2_output_prefix, "_h2_non_inf_auc_", pv[i], ".txt"))
	pT_cv2 = readLines(paste0(cv2_output_prefix, "_pT_non_inf_auc_", pv[i], ".txt"))
	res_cv1[1,i] = get_res(h2_cv1)
	res_cv1[2,i] = get_res(pT_cv1)
	res_cv2[1,i] = get_res(h2_cv2)
	res_cv2[2,i] = get_res(pT_cv2)
}

best1 = which(res_cv1 == max(res_cv1), arr.ind = TRUE)
best2 = which(res_cv2 == max(res_cv2), arr.ind = TRUE)

avg_cv = (res_cv1[best2[1],best2[2]]+res_cv2[best1[1],best1[2]])/2

cat("Average CV AUC/COR is ", avg_cv, '\n')
