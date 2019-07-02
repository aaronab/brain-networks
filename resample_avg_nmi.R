resample_avg_nmi <- function(x, 
                             group_membership, 
                             iterations = 10000, 
                             output_samples = F){
  
  avgcor <- function(x){mean(abs(x[lower.tri(x)]))}
  
  avg_nmi = rep(0, iterations)
  
  for(i in 1:iterations){
    pseudo_group = sample(group_membership)
    group_avg = avgcor(x[pseudo_group, pseudo_group])
    nongroup_avg = avgcor(x[!pseudo_group, !pseudo_group])
    avg_nmi[i] = mean(group_avg, nongroup_avg)
  }  
  
  if(output_samples == F) outputs <- c(mean(avg_nmi), 
                                       quantile(avg_nmi, c(.05, .95)) )
  else outputs <- avg_nmi
  return(outputs)
  
}