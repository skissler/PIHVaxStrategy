library(tidyverse) 
library(sf) 
library(scales)
library(tidycensus)

source("code/utils.R")
source("code/utils_private.R")

census_api_key(my_censuskey)

bostoncovid <- read_csv("data/bostoncovid.csv")
chicagocovid <- read_csv("data/chicagocovid.csv")
houstoncovid <- read_csv("data/houstoncovid.csv")
nyccovid <- read_csv("data/nyccovid.csv")
seattlecovid <- read_csv("data/seattlecovid.csv")