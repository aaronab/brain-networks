#R code for the analysis of graph community structure, from 
#Alexander-Bloch et al. Neuroimage. 2011. PMID: 22119652
#Alexander-Bloch et al. Cereberal Cortex. 2012. PMID: 22275481
#See graph_construction.R for graph construction methods
#email aalexanderbloch@gmail.com

pc <- function(G2, com){
#calculate participation coefficients of every node
#G2 is igraph graph, com is partition
#e.g. com<-spinglass.community(G2)$membership
	require(igraph)
	degs <- degree(G2)
	mat.id <- matrix(0, length(unique(com)), length(com))
	for(i in unique(com) + 1){
		for(j in 1:length(com)){
			mat.id[i,j] <- sum(com[(V(G2) [ nei( j-1 ) ]) +1] == i-1)
		}
		mat.id[i,] <- mat.id[i,]/degs
	}
	mat.id <- mat.id**2
	mat.id <- apply(mat.id,2,sum)
	mat.id <- 1-mat.id
	return(mat.id)
}

nid <- function(G2, com){
#calculate normalized intramodular degree of very node
#G2 is igraph graph, com is partition
	require(igraph)
	id <- rep(0, length(com)) #intramodular degree
	mid <- rep(0, length(com)) #mean intramodular degree
	sid <- rep(0, length(com)) #sd intramodular degree
	for(i in 1:length(com)){
		id[i] <- sum(com[(V(G2) [ nei( i-1 ) ]) +1] == com[i])
	}
	for(i in 1:length(com)){
		mid[i] <- mean(id[which(com == com[i])])
		sid[i] <- sd(id[which(com == com[i])])
	}
	nid <- (id-mid)/sid
	nid[is.na(nid)] <- 0
	return(nid)
}

part.nmi <- function(p1,p2){
#normalized mutal information between 2 partitions
	pp <- p1
	cc <- 1
	for(i in unique(p1)){
		pp[p1 == i] <- cc
		cc <- cc + 1
	}
	p1 <- pp
	pp <- p2
	cc <- 1
	for(i in unique(p2)){
		pp[p2 == i] <- cc
		cc <- cc + 1
	}
	p2 <- pp
	nn <- length(p1)
	n1 <- length(unique(p1))
	n2 <- length(unique(p2))
	if((n1 == 1) & (n2 == 1) ){frac <- 0} else {
		N <- matrix(0, nrow = n1, ncol = n2)
		for(i in 1:n1 ){
			for(j in 1:n2 ){
				N[i,j] <- sum ( (p1 == i) & (p2 == j) )
			}
		}
		num <- 0
		for(i in 1:n1 ){
			for(j in 1:n2 ){
				if(N[i,j] > 0){
					num <- num + (N[i,j] * log( (  N[i,j] * nn  )/(  sum(   N[i,]   ) *  sum(   N[,j]   )  ) ))
				}
			}
		}
		num <- -2 * num
		denomi <- 0
		for(i in 1:n1 ){
			denomi <- denomi + (sum(N[i,]) * log( sum(N[i,])/nn ) )
		}
		denomj <- 0
		for(j in 1:n2 ){
			denomj <- denomj + (sum(N[,j]) * log( sum(N[,j])/nn ) )
		}
		frac <- num/(denomi + denomj)
	}
	return(frac)
}

create.nmi.matrix.group <- function(x){
#find NMI between every pair of partitions, for a group of partitions
#x is a matrix of partitions for a group of subjects
	combined <- x
	bigmat <- matrix(0, nrow(x),nrow(x))
	for(i in 1:nrow(bigmat)){
		for(j in 1:nrow(bigmat)){
			bigmat[i,j] <- part.nmi(combined[i,], combined[j,])
		}
	}
	return(bigmat)
}

hungarianmatch <- function(g, d){
#match 2 partitions using the Hungarian algorithm
#registers d to g, i.e. first partition determines labels
#note that algorithm needs the communities to go 1:n to work propertly
require(clue)
if(all(g==d)){
return(g)} else{
if(length(unique(g)) <= length(unique(d))){
x <- g; y <- d
N <- T
} else{
y <- g; x <- d
N <- F
}
if(min(x) == 0){
x <- x+1}
if(min(y) == 0){
y <- y+1}
xx <- length(unique(x))
yy <- length(unique(y))
lu <- min(xx, yy)
lp <- max(xx, yy)
simmatij <- matrix(rep(0, xx*yy), nrow = xx)
for(i in 1:xx){for(j in 1:yy){simmatij[i, j] <-  sum((x == (i)) & (y == (j)))  }}
a <- solve_LSAP(simmatij, maximum = T)
newx <- x
newy <- y
for(i in 1:lu){
newx[x == i] <- a[[i]]
}
b <- list()
length(b) <- lp
bvec <-1:lp
for(i in 1:lu){
b[[ a[[i]] ]] <- i
bvec <- bvec[-which(bvec == a[[i]])]
}
if(length(bvec) > 0){
for(i in 1:length(bvec)){
b[[bvec[i]]] <- lu + i
}
}
for(i in 1:lp){
newy[y == i] <- b[[i]]
}
if(N){
return(newy)
} else{
return(newx)
}
}
}

median.match <- function(x){
#match a group of partitions to the median in terms of NMI
#for greedy optimization step, email aalexanderbloch@gmail.com
#x is a matrix of partitions for a group of subjects
#includes step to make all community labels go 1:n
	for(i in 1:nrow(x)){
		switch <- x[i,]
		uniquem <- unique(x[i,])
		for(j in 1:length(unique(x[i,]))){
			switch[x[i,] == uniquem[j]] <- j
		}
		x[i,] <- switch
	}
	findmed <- create.nmi.matrix.group(x)
	diag(findmed) <- 0
	findmed <- apply(findmed, 1, sum)
	findmed <- which.max(findmed)
	medianpart <- x[findmed,]
	bigtogether <- vector()
	for(i in 1:nrow(x)){
		bigtogether <- rbind(bigtogether, hungarianmatch(medianpart, x[i,]))
	}
	return(bigtogether)
}

