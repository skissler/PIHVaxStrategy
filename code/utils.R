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


make_toprow <- function(rho_l, rho_h, N){
	eval_pts <- seq(from=0, to=1, length.out=N)
	out <- rep(0, N)
	out[eval_pts<=0.5] <- rho_l - (rho_l-1)/0.5*eval_pts[eval_pts<=0.5]
	out[eval_pts>0.5] <- (1-(rho_h-1)) + (rho_h-1)/0.5*(eval_pts[eval_pts>0.5])
	return(out)
}

make_RRmat <- function(rho_l, rho_h, N){
	toprow <- make_toprow(rho_l, rho_h, N)
	bottomrow <- make_toprow(1/rho_l, 1/rho_h, N)
	out <- matrix(unlist(map2(toprow, bottomrow, function(x,y){seq(from=x,to=y,length.out=N)})),nrow=N)
	return(out)
}

assign_mortality <- function(rho_l, rho_h, N){

	mcpairs <- tibble(m=rep(0,N), c=rep(0,N))

	RRmat <- make_RRmat(rho_l, rho_h, N)
	for(m in sample(1:N)){
		c <- sample(1:N, 1, prob=RRmat[m,])
		RRmat[,c] <- 0
		mcpairs[m,] <- tibble(m=m, c=c)
	}

	# Need to reverse order of 'c' since RRmat[1,1] is high-mortality but low-contact. Want 1 to correspond to both high mortality and high contact.
	mcpairs$c <- (N+1)-mcpairs$c

	return(mcpairs)

}

# temp <- assign_mortality(1000,1000,101)

# plot(temp$c, temp$m)

# N <- 100
# ggplot(assign_mortality(1,1,N), aes(x=c/N, y=m/N)) + geom_point()


# 1:10




