library(tidyverse) 
library(lubridate) 

# Seattle: 
seattlecovid <- read_csv("data/raw/seattlecovid.csv") %>% 
	select(Location_Name, Positive_Rate) %>%
	mutate(Positive_Rate = case_when(Positive_Rate=="."~"0",TRUE~Positive_Rate)) %>%
	mutate(Positive_Rate = as.numeric(Positive_Rate)) %>% 
	mutate(Positive_Rate = Positive_Rate*100) %>%
	rename(zip=Location_Name) %>%
	rename(covidrate=Positive_Rate) %>%
	filter(zip!="All King County")

write_csv(seattlecovid, "data/seattlecovid.csv")


# Chicago: 
chicagocovid <- read_csv("data/raw/chicagocovid.csv") %>%
	select(`ZIP Code`, `Week End`, `Case Rate - Cumulative`) %>%
	rename(zip="ZIP Code") %>%
	rename(weekend="Week End") %>%
	rename(covidrate="Case Rate - Cumulative") %>% 
	mutate(weekend=mdy(weekend)) %>%
	filter(weekend==ymd("2021-02-27")) %>%
	select(-weekend) 

write_csv(chicagocovid, "data/chicagocovid.csv")

# Houston: 
houstoncovid <- read_csv("data/raw/houstoncovid.csv") %>%
	select(ZIP, TotalPop, TotalConfirmedCases) %>%
	mutate(zip=as.character(ZIP)) %>% 
	mutate(covidrate = TotalConfirmedCases/TotalPop*100000) %>%
	select(zip, covidrate) 

write_csv(houstoncovid, "data/houstoncovid.csv")


# New York:
nyccovid <- read_csv("data/raw/nyccovid.csv") %>%
	mutate(zip=as.character(MODIFIED_ZCTA), covidrate=COVID_CASE_RATE) %>%
	select(zip, covidrate)
write_csv(nyccovid, "data/nyccovid.csv")

