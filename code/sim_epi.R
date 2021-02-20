# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
source('code/utils.R')

# =============================================================================
# Set key parameters
# =============================================================================

R0 <- 3 			# Basic reproduction number
N <- 4000 			# Population size
nbin.size <- 1.68	# Neg.Bin. contact distribution size parameter
nbin.mean <- 10.47	# Neg.Bin. contact distribution mean

Emean <- 3 			# Mean exposure period (days)
Imean <- 4 			# Mean infectious period (days) 

tmax <- 90			# Max days to run the simulation

# A list that contains the start time for vaccination and the percent of the population at which to switch strategies. Possible strategies are "risk" (prioritize those at highest risk), "contact" (prioritize those with the most contacts), or "anyone" (vacinate at random), or "none" if vaccinattion is to stop. Also sets the amountt of infecttion blocking and transmission blocking that the vaccine gives. I think this makes it a 'leaky' vaccine. 

vax_strategy <- list(
	tstart=1,
	daily_nvax=floor(0.002*N),
	# daily_pvax=0.002,
	efficacy=0.95,
	infblock=0,
	transblock=0,
	pswitch=tibble(
		pstart=c(0,0.5,0.6), 
		pend=c(0.5,0.6,1.0), 
		prioritize=c("mortality","contact","none")))

# =============================================================================
# Calculate key quantities
# =============================================================================

# Infection probability per contact-day:
pinf <- R0/(nbin.mean*Imean)	

# Define contact weights for specifying adjacency matrix: 
Wi <- make_contact_weights(N=N, nbin.size=nbin.size, nbin.mean=nbin.mean)

# Define individual ifr
mortrisk <- sample(ifr$ifr, size=N, replace=TRUE, prob=popdist2019$prop)

# Set each person's vaccination priority by contacts and by mortality risk:
priority <- tibble(id=1:N, contacts=Wi, mortality=mortrisk) %>%
	sample_n(N) %>%
	arrange(desc(contacts)) %>%
	mutate(contact_priority=1:N) %>%
	sample_n(N) %>%
	arrange(desc(mortality)) %>%
	mutate(mortality_priority=1:N) %>%
	arrange(id) %>%
	select(-contacts, -mortality)

# =============================================================================
# Run simulation
# =============================================================================

# Initialize output tibble: 
casecounts <- tibble(t=1:tmax,
	E=rep(0,tmax),	# Exposed
	I=rep(0,tmax),	# Infectious
	R=rep(0,tmax),	# Recovered
	X=rep(0,tmax),	# Deceased
	V=rep(0,tmax),	# Vaccinated
	strategy="none" # Vaccination strategy
	)
priority_axis <- "none"

# Initialize compartment vectors:
Evec <- matrix(rep(0,N),ncol=1)
Ivec <- matrix(rep(0,N),ncol=1)
Ivec[1] <- 1
Rvec <- matrix(rep(0,N),ncol=1)
Xvec <- matrix(rep(0,N),ncol=1)
Vvec <- matrix(rep(0,N),ncol=1)

# Fill output tibble with initial compartment sums:
casecounts$E[1] <- sum(Evec)
casecounts$I[1] <- sum(Ivec)
casecounts$R[1] <- sum(Rvec)
casecounts$X[1] <- sum(Xvec)
casecounts$V[1] <- sum(Vvec)

# Increment time:
t <- 2

# Define adjacency matrix here if it's meant to be static:
A <- make_adj(weights=Wi)

while(t<=tmax){

	# Define adjacency matrix here if it's meant to be dynamic:
	# A <- make_adj(weights=Wi) 

	# Set vaccine priority
	if(t>=(vax_strategy$tstart)){
		priority_axis <- vax_strategy$pswitch %>% 
			filter((sum(Vvec)/N >= pstart) & (sum(Vvec)/N < pend)) %>%
			slice(1) %>%
			pull(prioritize)
		if(priority_axis == "mortality"){
			priority <- arrange(priority, mortality_priority)
		} else if(priority_axis == "contact"){
			priority <- arrange(priority, contact_priority)
		} else if(priority_axis == "anyone"){
			priority <- sample_n(priority, nrow(priority))
		}
	}

	# Propose transitions	
	newE <- (runif(N)<(1-(1-pinf)^(A%*%Ivec)))*(1-Evec)*(1-Ivec)*(1-Rvec)*(1-Xvec) # New S -> E
	newI <- (runif(N)<(1/Emean*Evec))*1 # New E -> I
	newR <- (runif(N)<(1/Imean*Ivec))*1	# New I -> R (though some will split off to X)
	newX <- (runif(N)<(newR*mortrisk*(1-Vvec*vax_strategy$efficacy)))*1 # New I -> X
	newR <- newR - newX # Revised I -> R, accounting for mortality

	# Propose new vaccinated: 
	if((t>=(vax_strategy$tstart)) & (priority_axis != "none") & (nrow(priority)>0)){
		newV <- matrix(rep(0,N),ncol=1)
		newV[priority$id[1:min(vax_strategy$daily_nvax, nrow(priority))]] <- 1
	} else {
		newV <- matrix(rep(0,N),ncol=1)
	}

	# Run transitions
	Evec <- Evec + newE - newI
	Ivec <- Ivec + newI - newR
	Rvec <- Rvec + newR
	Xvec <- Xvec + newX
	Vvec <- Vvec + newV

	# Remove vaccinated and deceased people from the vaccine queue: 
	priority <- filter(priority, !(id %in% which(newV==1)))
	priority <- filter(priority, !(id %in% which(newX==1)))

	# Fill in output array: 
	casecounts$E[t] <- sum(Evec)
	casecounts$I[t] <- sum(Ivec)
	casecounts$R[t] <- sum(Rvec)
	casecounts$X[t] <- sum(Xvec)
	casecounts$V[t] <- sum(Vvec)
	casecounts$strategy[t] <- priority_axis

	# Increment time: 
	t <- t+1

}

# Plot output: 
fig_casecounts <- casecounts %>% 
	pivot_longer(c("E","I","R","X","V")) %>%
	ggplot(aes(x=t, y=value, col=name)) + 
		geom_point(size=0.5) + 
		geom_line() + 
		scale_y_continuous(limits=c(0,N)) + 
		theme_minimal() + 
		labs(x="Day", y="Cases",col="Compartment")
ggsave(fig_casecounts, file="figures/casecounts.png", width=8, height=5)






