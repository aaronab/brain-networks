#R code for graph construction methods used in Alexander-Bloch et al. Frontiers in Systems Neuroscience. 2010.
#For questions or comments, or to upload code of your own, email aalexanderbloch@gmail.com

############################## command sequence for each subject would be something like as follows ##########################
#where corTs2 is the functional connectivity matrix
#need igraph (www.igraph.sourceforge.net) for the mst function to work
#need to have defined the functions at the bottom of this document

library(igraph)
corTs2 <- abs(corTs2)
diag(corTs2) <- 0
knn.2 <- knn(corTs2)
mst.2 <- mst(corTs2)
lt.2 <- combine.mst(mst.2, knn.2)
gt.2 <- gt(corTs2)

#at this point lt.2 contains the edges in the order in which they should be added for local thresholding
#gt.2 contains the edges in order for global thresholding
#you use extract to look at different numbers of edges i, e.g. below for a local threshold

ic <- extract(corTs2, lt.2, i)
icb <- ic; icb[,] <- 0; icb[ic > 0] <- 1
G <- graph.adjacency(icb, mode = 'lower', diag = F)

#ic is a weighted graph
#icb is a binary graph, in the format that brainwaver (http://cran.r-project.org/web/packages/brainwaver/) wants
#G is the graph in igraph format

############################################################ functions ########################################################
#copy and paste below into R

knn  <- function(graph){
#graph=functional connectivity matrix
#returns the order in which edges should be added to make k nearest neighbors graphs (k-NNG)
#i.e. first, edges of 1-NNG, then 2-NNG, etc until all edges have been added
#within each k-NNG cycle, edges are ordered in terms of their connectivity
	graphi <- graph
	thcg <- graph[,]
	thcg[,] <- 0
	diag(graphi) <- 0
	diag(graph) <- 0
	thcgi <- 1
	k <- 0
	as <- vector()
	while(sum(graph) > 0){
		k <- k + 1
		a <- cbind(sort(apply(graphi, 1, function(x) sort(x, decreasing = T))[k,], 
						index.return = T, decreasing = T)$ix, (apply(graphi, 1, function(x) (sort(x, 
						decreasing = T, index.return = T)$ix))[k,])[sort(apply(graphi, 1, function(x) 
						sort(x, decreasing = T))[k,], index.return = T, decreasing = T)$ix])[ sort(apply(graphi, 1,
						function(x) sort(x, decreasing = T))[k,], index.return = T, decreasing = T)$x  > 0 ,]
		if(is.matrix(a)){
			aa <- rbind(a, cbind(a[,2], a[,1]))
		} else{
			aa <- rbind(aa, aa[c(2,1)])
		}		
		if(length(as) > 1){
			ad <- rbind(aa, as, cbind(as[,2], as[,1]))
		} else{
			ad <- aa
		}		
		if(is.matrix(aa)){
			iii <- c(1:(length(aa[,1])/2), (1:(length(aa[,1])/2) + 
						.5))[(!duplicated(ad, fromLast = T))[1:length(aa[,1])]]
			aa <- ad[(!duplicated(ad, fromLast = T))[1:length(aa[,1])],]
			iii <- sort(iii, index.return = T)$ix
			aa <- aa[iii,]			
		} else{
			aa <- ad[(!duplicated(ad, fromLast = T))[1,],]
		}	
		if(length(aa) > 2){
			aa <- aa[aa[,1] > aa[,2],]
		}
		as <- rbind(as, aa)
		if(length(aa) == 2){
			thcg[aa[1], aa[2]] <- thcgi
			thcg[aa[2], aa[1]] <- thcgi
			thcgi <- thcgi + 1
		} 
		if(length(aa) > 2){		
			for(i in 1:nrow(aa)){
				thcg[aa[i,1], aa[i,2]] <- thcgi
				thcg[aa[i,2], aa[i,1]] <- thcgi
				thcgi <- thcgi + 1
		}  	}
		graph[thcg > 0] <- 0
	}
	return(thcg)
}

# renamed the fn below to avoid clashes with igraph::mst
fcm_mst <- function(fcm){
  #fcm = functional connectivity matrix
  #computes the minimum spanning tree of the graph (using igraph function)
  #returns these edges ordered in terms of their connectivity
  corTs2 <- abs(fcm)
  diag(corTs2) <- 0
  thc <- corTs2
  thc[,] <- 0
  g3 <- graph.adjacency(1-corTs2, weighted = T, mode = 'lower', diag = F)
  G4 <- minimum.spanning.tree(g3)
  for(i in 1:(nrow(thc))){
    print(i)
    thc[i, neighbors(G4, i)] <- 1
  }
  corTs2 <- corTs2*thc
  diag(corTs2) <- 0
  corTs2[upper.tri(corTs2)] <- 0
  index <- sort(corTs2, index.return = T, decreasing = T)
  ttt <- corTs2
  ttt[,] <- 0
  ttt[index[[2]][1:sum(corTs2 > 0)]] <- 1:sum(corTs2 > 0)
  ttt[upper.tri(ttt)] <- t(ttt)[upper.tri(ttt)]
  return(ttt)
}

combine.mst <- function(g.mst = tm, g.o = ata){
#g.mst = output from mst(); g.o = output from knn()
#reorders so as to start with mst and then add edges of the k-NNG
	g.r <- g.mst
	i.mst <- g.mst > 0
	g.o[i.mst] <- Inf
	diag(g.o) <- Inf
	g.o[upper.tri(g.o)] <- Inf
	index <- sort(g.o, index.return = T)
	g.ra <- g.r
	g.ra[,] <- 0
	g.ra[index$ix[1:sum(!is.infinite(g.o))]] <- (((sum(g.mst > 0))/2) + 1):(((sum(g.mst > 0))/2) + 
					length(index$ix[1:sum(!is.infinite(g.o))]))
	g.ra[upper.tri(g.ra)] <- t(g.ra)[upper.tri(g.ra)]
	g.r <- g.r + g.ra
	return(g.r)
}

extract <- function(graph = abs(corTs2), index = g.r, edge = i){
#graph=functional connectivity matrix
#index=order of edges, e.g. output from combine.mst()
#edge=number of edges to include in graph
#returns a weighted adjacency matrix
	diag(graph) <- 0
	newindex <- index < (edge + 1)
	thc <- graph
	thc[,] <- 0
	thc[newindex] <- 1
	graphnew <- graph*thc
	return(graphnew)
}

gt <- function(graph){
#graph=functional connectivity matrix
#returns the order in which edges should be added for standard global threhsolding
	corTs2 <- abs(graph)
	corTs2[upper.tri(corTs2)] <- 0
	diag(corTs2) <- 0
	index <- sort(corTs2, index.return = T, decreasing = T)
	ttt <- corTs2
	ttt[,] <- 0
	ttt[index[[2]][1:sum(corTs2 > 0)]] <- 1:sum(corTs2 > 0)
	ttt[upper.tri(ttt)] <- t(ttt)[upper.tri(ttt)]
	return(ttt)
}

