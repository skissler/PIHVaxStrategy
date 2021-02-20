library(tidyverse) 


makeSymm <- function(m) {
   m[upper.tri(m)] <- t(m)[upper.tri(m)]
   return(m)
}

# make_adj <- function(N=4000, nbin.size=1.68, nbin.mean=10.47){
# 	Wi <- rgamma(N, shape=nbin.size, scale=nbin.mean/nbin.size)/sqrt(nbin.mean)
# 	ui <- Wi/sqrt(N)
# 	pmat <- ui%o%ui / (1 + ui%o%ui)
# 	pmat[upper.tri(pmat,diag=TRUE)] <- 0
# 	drawmat <- matrix(runif(N^2),nrow=N)
# 	adj <- 1*(drawmat < pmat)
# 	adj <- makeSymm(adj)
# }

make_contact_weights <- function(N=4000, nbin.size=1.68, nbin.mean=10.47){
	Wi <- rgamma(N, shape=nbin.size, scale=nbin.mean/nbin.size)/sqrt(nbin.mean)
	return(Wi)
}

make_adj <- function(weights){
	Wi <- weights
	N <- length(Wi)
	ui <- Wi/sqrt(N)
	pmat <- ui%o%ui / (1 + ui%o%ui)
	pmat[upper.tri(pmat,diag=TRUE)] <- 0
	drawmat <- matrix(runif(N^2),nrow=N)
	adj <- 1*(drawmat < pmat)
	adj <- makeSymm(adj)
	return(adj)
}
