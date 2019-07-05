resample_avg_nmi <- function(x, 
                             group_membership, 
                             iterations = 10000){
  
  avgcor <- function(x){mean(abs(x[lower.tri(x)]))}
  
  avg_nmi = rep(0, iterations)
  
  for(i in 1:iterations){
    pseudo_group <- sample(group_membership)
    group_avg <- avgcor(x[pseudo_group, pseudo_group])
    nongroup_avg <- avgcor(x[!pseudo_group, !pseudo_group])
    avg_nmi[i] <- mean(group_avg, nongroup_avg)
  }  
  
  resample_nmi_avg <- c(mean(avg_nmi), 
                        quantile(avg_nmi, c(.05, .95)) )
  
  orig_group_avg <- rep(0, 3)
  orig_group_avg[1] <- avgcor(x[group_membership, group_membership])
  orig_group_avg[2] <- avgcor(x[!group_membership, !group_membership])
  orig_group_avg[3] <- mean(orig_group_avg[1:2])
  
  pval <- sum(avg_nmi >  orig_group_avg[3] )/ iterations
  
  output = c(orig_group_avg, resample_nmi_avg, pval)
  
  names(output) = c("mean_group1", "mean_group2", "mean_grouped",
                    "resample_mean", "CI_Lower", "CI_Upper",
                    "pval")
    
  return(output)
  
}
