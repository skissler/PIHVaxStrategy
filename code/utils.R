library(tidyverse) 

makeSymm <- function(m) {
   m[upper.tri(m)] <- t(m)[upper.tri(m)]
   return(m)
}

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

# https://data.census.gov/cedsci/table?q=population&tid=ACSST1Y2019.S0101&hidePreview=false
popdist2019 <- tibble(
	agegrp=c("<5",
		"5-9",
		"10-14",
		"15-19",
		"20-24",
		"25-29",
		"30-34",
		"35-39",
		"40-44",
		"45-49",
		"50-54",
		"55-59",
		"60-64",
		"65-69",
		"70-74",
		"75-79",
		"80+"),
	pop=c(19404835,
		19690437,
		21423479,
		21353524,
		21468680,
		23233299,
		22345176,
		21728259,
		20186586,
		20398226,
		20464881,
		21484060,
		20984053,
		17427013,
		14148548,
		9759764,
		12738703)) %>%
	mutate(prop=pop/sum(pop))

ifr <- tibble(
	agegrp=c("<5",
		"5-9",
		"10-14",
		"15-19",
		"20-24",
		"25-29",
		"30-34",
		"35-39",
		"40-44",
		"45-49",
		"50-54",
		"55-59",
		"60-64",
		"65-69",
		"70-74",
		"75-79",
		"80+"),
	ifr=c(2.912324e-05,
		6.214030e-06,
		9.250515e-06,
		2.543196e-05,
		6.444111e-05,
		1.320592e-04,
		2.426553e-04,
		4.027768e-04,
		7.533561e-04,
		1.207733e-03,
		2.073839e-03,
		3.231949e-03,
		4.571147e-03,
		1.076168e-02,
		1.676234e-02,
		3.206417e-02,
		8.301160e-02))







