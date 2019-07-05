# For a given node X, the other 277 nodes are relabeled (86 for us)
# to reflect simply whether or not they belong to X’s module. 
# 
# These labels can then be compared across subjects: 
# The similarity of two subjects, in terms of node X’s functional community, 
# can be quantified as Pearson’s phi, − 1≤φ≤1, a statistic that is essentially 
# the Pearson correlation of a dichotomous variable (Pearson, 1900). 
# 
# The phi coefficient can then serve the same role for node X that 
# NMI served for the whole network. If there is a genuine difference in node X’s 
# functional community between the patients and the controls, then the average 
# within group phi coefficient should be higher than the average between group phi coefficients. 

make_comb_mat <- function(vec, CombLength, Perms = F, Repeats = F){
  # vec # vector of numbers from which to create perms
  # PermLength # length of each perm
  # Perms = F # if T, return permutations, if F, return combinations
  # Repeats = F # if T allow repetition of values (e.g., 1-1-1), if F don't allow
  
  if(!length(vec) > CombLength) print("vector of combinates shorter than Combination Length")
  
  lst = lapply(numeric(CombLength), function(x) vec)
  mat = as.matrix(expand.grid(lst))
  
  if(Repeats == F) mat = mat[colSums(apply(mat,1,duplicated))==0,] 
  # remove rows including repeats of the same value
  
  if(Perms == F){
    sort.rows = apply(mat,1,sort)
    row.content = apply(sort.rows,2,paste, collapse = "")
    mat = mat[!duplicated(row.content), ]
  } # perms False
  return(mat)
}

calc_node_phi <- function(x){
  # x is matrix with partitions as rows, nodes as columns and module labels as the data
  require(psych) # make sure psych is installed for phi
  
  # Step 1. set up 3d matrix to hold phis for each node for each pair
  node_phi_array <- array(0, dim=c(nrow(x), nrow(x), ncol(x)))
  
  # Step 2. generate a list of all pairs - more than halves the number of calculations
  pairs_mat <- make_comb_mat(1:nrow(x), CombLength = 2)
  
  for(node in 1:ncol(x)){
    # For a given node X, the other nodes are relabeled
    # to reflect simply whether or not they belong to X’s module. 
    # node <- 1
    print(paste("Calculating pairs phi for node", node, sep = " ") )
    relabelled_mat <- all_mat[, -node] == all_mat[, node]
  
    # The similarity of two subjects, in terms of node X’s functional community, 
    # can be quantified as Pearson’s phi, − 1≤φ≤1, a statistic that is essentially 
    # the Pearson correlation of a dichotomous variable (Pearson, 1900). 
    
    # Step 3. calculate phi for every pair in this node
    for(i in 1:nrow(pairs_mat) ){

      # check for lonely nodes
      if(sum((relabelled_mat[pairs_mat[i,1],])) + 
         sum((relabelled_mat[pairs_mat[i,2],])) == 0) phi_pair = 0 
      else {

      phi_pair <- phi( matrix ( rbind( table(relabelled_mat[pairs_mat[i,1],]), 
                                       table(relabelled_mat[pairs_mat[i,2],]) ), nrow = 2 ) )
      node_phi_array[pairs_mat[i,1], pairs_mat[i,2], node] <- phi_pair
      node_phi_array[pairs_mat[i,1], pairs_mat[i,2], node] <- phi_pair
      } # lonely nodes
      
    } # i
    
  }# node
  
  return(node_phi_array)

}

#  currently only works for two groups
resample_avg_node_phi <- function(x, 
                             group_membership, 
                             iterations = 10000){
  
  # x is a 3d aray produced by calc_node_phi
  # group_membership is a logical vector indicating group membership 
  # currently only works for two groups

  avgcor <- function(x){mean(abs(x[lower.tri(x)]))}
  
  num_nodes <- dim(x)[3]
  
  avg_node_phi = matrix(0, nrow = num_nodes, ncol = iterations)
  
  for(i in 1:iterations){
    node_avg = cbind(1:86, 0, 0)
    # get average phi for each node for pseudo groups
    pseudo_group = sample(group_membership)
    for(node in 1:86) {
      node_avg[node, 2] = avgcor(x[pseudo_group, pseudo_group, node])
      node_avg[node, 3] = avgcor(x[!pseudo_group, !pseudo_group, node])
    }
    
    # average across groups and add to matrix
    avg_node_phi[, i] <-  rowMeans(node_avg[,2:3])

  }  
  
  
  CI_Lower <- apply(avg_node_phi, 1, quantile, probs = .05)
  CI_Upper <- apply(avg_node_phi, 1, quantile, probs = .95)
  
  resample_node_avg <- cbind(1:86, rowMeans(avg_node_phi), CI_Lower, CI_Upper)
  colnames(resample_node_avg)[1:2] = c("node", "resample_mean")
  
  group_node_avg = cbind(1:86, 0, 0, 0)
  # get average phi for each node for original groups
  
  for(node in 1:86) {
    group_node_avg[node, 2] = avgcor(x[group_membership, group_membership, node])
    group_node_avg[node, 3] = avgcor(x[!group_membership, !group_membership, node])
  }
  
  group_node_avg[,4] = rowMeans(group_node_avg[,2:3])
  colnames(group_node_avg) = c("node", "mean_group1", 
                                    "mean_group2", "mean_grouped")
  
  grouped_pval = cbind(1:86, 0, 0)
  for(node in 1:86) {
    grouped_pval[node, 2] = sum(avg_node_phi[node,] >group_node_avg[node,4] )/ iterations
  }
  grouped_pval[,3] = round(p.adjust(grouped_pval[,2], method = "fdr"),4)
  colnames(grouped_pval) = c("node", "pval", "fdr_pval")
  
  output = cbind(group_node_avg, resample_node_avg[,2:4], grouped_pval[,2:3])
  
  return(output)
  
}

avgcor <- function(x){mean(abs(x[lower.tri(x)]))}

## -------------------- Import data ---------------------------------

# read in data
# to calculate node_phi
# x is matrix with partitions as rows, nodes as columns and module labels as the data
# example below assumes data includes more than the partition info

## -------------------- Working ---------------------------------

# put module assignments in matrix
all_mat <- as.matrix(data[, 3:88]) # extract module assignments

node_phi_array = calc_node_phi(all_mat)

# generate group_membership logical vectors
control_rows = data$Group=="C"
# patient_rows = data$Group=="P"

rnd_avg_phi <- resample_avg_node_phi(node_phi_array, 
                                group_membership = control_rows, 
                                iterations = 10000)

hist(rnd_avg_phi[,4] )
plot(rnd_avg_phi[,1], rnd_avg_phi[,4], xlab = "node", ylab = "average phi" )
points(rnd_avg_phi[,1], rnd_avg_phi[,5], col = "red") # red dots for resampled mean



