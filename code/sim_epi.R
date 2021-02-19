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

# =============================================================================
# Calculate key quantities
# =============================================================================

pinf <- R0/(nbin.mean*Imean)	# Infection probability per contact-day

# =============================================================================
# Run simulation
# =============================================================================

casecounts <- tibble(t=1:tmax,	# Initialize output tibble
	E=rep(0,tmax),
	I=rep(0,tmax),
	R=rep(0,tmax))	

Evec <- matrix(rep(0,N),ncol=1)		# Initialize exposed vector	
Ivec <- matrix(rep(0,N),ncol=1)		# Initialize infectious vector
Ivec[2] <- 1						# Initialize initial infected
Rvec <- matrix(rep(0,N),ncol=1)		# Initialize recovered vector

casecounts$E[1] <- sum(Evec)
casecounts$I[1] <- sum(Ivec)
casecounts$R[1] <- sum(Rvec)
t <- 2

# Define adjacency matrix here if it's meant to be static:
A <- make_adj(N=N, nbin.size=nbin.size, nbin.mean=nbin.mean)

while(t<=tmax & (sum(Evec)+sum(Ivec))>0){

	# Define adjacency matrix here if it's meant to be dynamic:
	# A <- make_adj(N=N, nbin.size=nbin.size, nbin.mean=nbin.mean) 

	newE <- (runif(N)<(pinf*A%*%Ivec))*(1-Evec)*(1-Ivec)*(1-Rvec)
	newI <- (runif(N)<(1/Emean*Evec))*1
	newR <- (runif(N)<(1/Imean*Ivec))*1

	Evec <- Evec + newE - newI
	Ivec <- Ivec + newI - newR
	Rvec <- Rvec + newR

	casecounts$E[t] <- sum(Evec)
	casecounts$I[t] <- sum(Ivec)
	casecounts$R[t] <- sum(Rvec)

	t <- t+1

}
casecounts[(t-1):tmax,-1] <- casecounts[t-1,-1]

fig_casecounts <- casecounts %>% 
	pivot_longer(c("E","I","R")) %>%
	ggplot(aes(x=t, y=value, col=name)) + 
		geom_point() + 
		geom_line() + 
		scale_y_continuous(limits=c(0,N)) + 
		theme_minimal() + 
		labs(x="Day", y="Cases",col="Compartment")


