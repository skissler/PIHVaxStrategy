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

# Infection probability per contact-day:
pinf <- R0/(nbin.mean*Imean)	

# Define contact weights for specifying adjacency matrix: 
Wi <- make_contact_weights(N=N, nbin.size=nbin.size, nbin.mean=nbin.mean)

# =============================================================================
# Run simulation
# =============================================================================

# Initialize output tibble: 
casecounts <- tibble(t=1:tmax,	
	E=rep(0,tmax),
	I=rep(0,tmax),
	R=rep(0,tmax))	

# Initialize compartment vectors:
Evec <- matrix(rep(0,N),ncol=1)
Ivec <- matrix(rep(0,N),ncol=1)
Ivec[1] <- 1
Rvec <- matrix(rep(0,N),ncol=1)

# Fill output tibble with initial compartment sums:
casecounts$E[1] <- sum(Evec)
casecounts$I[1] <- sum(Ivec)
casecounts$R[1] <- sum(Rvec)

# Increment time:
t <- 2

# Define adjacency matrix here if it's meant to be static:
A <- make_adj(weights=Wi)

while(t<=tmax & (sum(Evec)+sum(Ivec))>0){

	# Define adjacency matrix here if it's meant to be dynamic:
	# A <- make_adj(weights=Wi) 

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
# Fill in remaining time if necessary:
casecounts[(t-1):tmax,-1] <- casecounts[t-1,-1]

fig_casecounts <- casecounts %>% 
	pivot_longer(c("E","I","R")) %>%
	ggplot(aes(x=t, y=value, col=name)) + 
		geom_point(size=0.5) + 
		geom_line() + 
		scale_y_continuous(limits=c(0,N)) + 
		theme_minimal() + 
		labs(x="Day", y="Cases",col="Compartment")


