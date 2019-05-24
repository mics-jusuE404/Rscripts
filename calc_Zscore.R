## Calculate row-wise z-score: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
 
data_subset_norm <- t(apply(dataframe, 1, cal_z_score))
