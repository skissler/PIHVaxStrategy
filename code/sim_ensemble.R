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

nsims <- 500		# Set number of simulations in ensemble

# A list that contains the start time for vaccination and the percent of the population at which to switch strategies. Possible strategies are "risk" (prioritize those at highest risk), "contact" (prioritize those with the most contacts), or "anyone" (vacinate at random), or "none" if vaccinattion is to stop. Also sets the amountt of infecttion blocking and transmission blocking that the vaccine gives. I think this makes it a 'leaky' vaccine. 

# vax_strategy <- list(
# 	tstart=30,
# 	daily_nvax=floor(0.002*N),
# 	efficacy=0.95,
# 	infblock=0.5,
# 	transblock=0.5,
# 	pswitch=tibble(
# 		pstart= c(0.0, 0.1, 0.6), 
# 		pend  = c(0.1, 0.6, 1.0), 
# 		prioritize=c("mortality","contact","none")))

vax_strategy_list <- list(
	# ----------------------------------------------
	earlycontact_notrans=list(
	tstart=15,
	daily_nvax=floor(0.002*N),
	efficacy=0.95,
	infblock=0.0,
	transblock=0.0,
	pswitch=tibble(
		pstart= c(0.0, 0.1, 0.6), 
		pend  = c(0.1, 0.6, 1.0), 
		prioritize=c("mortality","contact","none"))),
	# ----------------------------------------------
	latecontact_notrans=list(
	tstart=15,
	daily_nvax=floor(0.002*N),
	efficacy=0.95,
	infblock=0.0,
	transblock=0.0,
	pswitch=tibble(
		pstart= c(0.0, 0.5, 0.6), 
		pend  = c(0.5, 0.6, 1.0), 
		prioritize=c("mortality","contact","none"))),
	# ----------------------------------------------	
	earlycontact_withtrans=list(
	tstart=15,
	daily_nvax=floor(0.002*N),
	efficacy=0.95,
	infblock=0.5,
	transblock=0.5,
	pswitch=tibble(
		pstart= c(0.0, 0.1, 0.6), 
		pend  = c(0.1, 0.6, 1.0), 
		prioritize=c("mortality","contact","none"))),
	# ----------------------------------------------
	latecontact_withtrans=list(
	tstart=15,
	daily_nvax=floor(0.002*N),
	efficacy=0.95,
	infblock=0.5,
	transblock=0.5,
	pswitch=tibble(
		pstart= c(0.0, 0.5, 0.6), 
		pend  = c(0.5, 0.6, 1.0), 
		prioritize=c("mortality","contact","none")))
	)

# =============================================================================
# Run sims
# =============================================================================

casecounts_full <- tibble()

# casecountslist <- list()
for(vs in 1:length(vax_strategy_list)){

	vax_strategy <- vax_strategy_list[[vs]]
	print(paste0("Vax Strategy ", vs))

	for(indexA in 1:nsims){
		source('code/sim_epi.R')
		casecounts_full <- bind_rows(casecounts_full, 
			mutate(casecounts, vax_strategy=vs, sim=indexA))
		# casecountslist[[indexA]] <- casecounts
		print(paste0("   Sim ",indexA))
	}
}

# casecounts <- bind_rows(casecountslist,.id="sim") %>% mutate(sim=as.numeric(sim))
# rm(casecountslist)

# =============================================================================
# Summarize output
# =============================================================================

casecounts_summary <- summarize_casecounts(casecounts_full)
write_csv(casecounts_summary, file="output/casecounts_summary.csv")

casecounts_summary %>% 
	filter(finalsize>50) %>%
	ggplot(aes(x=finalsize)) + 
		geom_histogram(bins=30)