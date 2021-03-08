library(tidyverse)
library(statnet)

# Examine ?'ergm-terms' for degree distribution statistic (is there one?)

# Simulatig a network with arbitrary degree: 
g <- network.initialize(5, directed=FALSE)
g.fit <- ergm(g~degree(1:5), target.stats=c(0,5,0,0,0))
g.sim <- simulate(g.fit)
plot(g.sim)
degreedist(g.sim)


# cut off at 50
popsize <- 100
g <- network.initialize(popsize, directed=FALSE)
sampled.degrees <- rnbinom(n=popsize, size=1.68, mu=10.47)
sampled.degree.dist <- unlist(lapply(1:50, function(x){sum(sampled.degrees==x)}))
g.fit <- ergm(g~degree(1:50), target.stats=sampled.degree.dist)

temp <- data.frame(draw=rnbinom(n=10000, size=1.68, mu=10.47))
ggplot(temp, aes(x=draw)) + 
	geom_density()

# =============================================================================
# A hard about-face to generate some of my own networks, using theory from Britton et al (https://staff.math.su.se/tom.britton/papers/fulltext.pdf) - "The generalized random graph" 

N <- 4000
nbin.size <- 1.68
nbin.mean <- 10.47
Wi <- rgamma(N, shape=nbin.size, scale=nbin.mean/nbin.size) / sqrt(nbin.mean)
ui <- Wi/sqrt(N)
pmat <- ui%o%ui / (1 + ui%o%ui)
pmat[upper.tri(pmat,diag=TRUE)] <- 0

makeSymm <- function(m) {
   m[upper.tri(m)] <- t(m)[upper.tri(m)]
   return(m)
}

drawmat <- matrix(runif(N^2),nrow=N)
adj <- 1*(drawmat < pmat)
adj <- makeSymm(adj)

temp <- tibble(sim=colSums(adj), ref=rnbinom(n=N,size=nbin.size,mu=nbin.mean))

temp %>% 
	pivot_longer(everything()) %>%
	ggplot(aes(x=value,col=name)) + 
		geom_density(adjust=2) + 
		theme_minimal()