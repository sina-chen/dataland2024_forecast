#-------------------------------------------------------------------------------
#
# Prepare data for forecasting model for provoince-level elections 1985-2024
#
# Author: Sina Chen
#
#
#-------------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------

{
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
}


# Data --------------------------------------------------------------------

setwd("your working directory")

# polls
polls1984_2023 <- read_csv("dataland/dataland_polls_1984_2023.csv")
polls2024 <- read_csv("dataland/dataland_polls_2024_scenarios.csv")

# election results
res1984_2023 <- read_csv("dataland/dataland_election_results_1984_2023.csv")

# calendar
calendar <- read_csv("dataland/dataland_electoral_calendar.csv")

# 2024 poll scenario, select A to E
sc <- "A"


# Preparation -------------------------------------------------------------

# filter province level polls &
#   remove 1984 polls since we have no previous election result as start value for random walk
polls1985_2023 <- polls1984_2023 %>% 
  filter(!(geography %in% c("National", "Metaflux Realm", 
                            "Circuit Confederation",
                            "Synapse Territories")) &
           year >= 1985)

# filter 2024 scenario & province level polls
polls2024 <- polls2024 %>% 
  filter(!(geography %in% c("National", "Metaflux Realm", 
                            "Circuit Confederation",
                            "Synapse Territories")) &
           scenario == sc)

# add 2024 polls
polls1985_2024 <- rbind(polls1985_2023, polls2024 %>% 
                          select(-scenario))
rm(polls1984_2023, polls1985_2023, polls2024)

# add election date
polls1985_2024 <- merge(polls1985_2024, calendar %>% 
                          select(election_cycle, election_day), 
                        by.x = "year", by.y = "election_cycle")
rm(calendar)

# order data frame by year and geography
polls1985_2024 <- polls1985_2024 %>% 
  arrange(year, geography)

# compute days to election of each poll
polls1985_2024 <- polls1985_2024 %>% 
  mutate(dte = difftime(election_day, date_conducted),
         election = paste0(year, geography),
         r_id = as.integer(as.factor(election)),
         st_id = if_else(geography %in% c("Cerebrica", "Cortexia", "Neuronia"),
                         1, 0),
         t_id =  as.numeric(max(dte) - dte) + 1)

# polls in long
polls1985_2024_long <- polls1985_2024 %>% 
  mutate(poll_id = paste0(year, geography, date_conducted, pollster)) %>% 
  pivot_longer(cols = c(cc_poll_share, dgm_poll_share, pdal_poll_share, 
                        ssp_poll_share),
               names_to = "party", values_to = "poll_share") %>% 
  mutate(party = str_remove_all(party, "_poll_share"))  %>% 
  rename("province" = geography)

# add previous election result
prev_res1985_2024_long <- res1984_2023 %>% 
  select(year, province, cc_share, dgm_share, pdal_share, ssp_share) %>% 
  pivot_longer(cols = c(cc_share, dgm_share, pdal_share, ssp_share), 
               names_to = "party", values_to = "prev_vote_share") %>% 
  mutate(year = year + 1,
         party = str_remove_all(party, "_share"))

# add previous election vote share
polls1985_2024_long <- merge(polls1985_2024_long, prev_res1985_2024_long,
                             by = c("year", "province", "party"), all.x = T)
rm(prev_res1985_2024_long, res1984_2023)

# order polls
polls1985_2024_long <- polls1985_2024_long %>% 
  arrange(year, province, poll_id, party, -dte)

# add party and poll id
polls1985_2024_long <- polls1985_2024_long  %>% 
  mutate(k_id = as.numeric(as.factor(party)),
         p_id = as.numeric(as.factor(poll_id)))

# election-party-level data
election_party_data <- polls1985_2024_long %>% 
  group_by(year, province, party, prev_vote_share) %>% 
  summarise(n=n()) %>% 
  ungroup() 

# order election-party data
election_party_data <- election_party_data %>% 
  arrange(year, province, party)

# election-level data
election_data <- polls1985_2024_long %>% 
  group_by(year, province, election, st_id) %>% 
  summarise() %>% 
  ungroup()

# order election data
election_data <- election_data %>% 
  arrange(year, province)

# add year and state id
election_data <- election_data %>%  
  mutate(y_id = as.integer(as.factor(year)),
         s_id = as.integer(as.factor(province)))

# add Synapse Territories id
election_data <- election_data  %>% 
  group_by(st_id) %>% 
  mutate(st_r_id = sequence(n())) %>% 
  ungroup()

# generate 0 matrices for redundant parametrization of mvn distributed party coefficients
zero_r <- matrix(nrow = length(unique(polls1985_2024_long$election)), 
                 ncol = length(unique(polls1985_2024_long$party))-1, 0)
zero_rK <- matrix(nrow = length(unique(polls1985_2024_long$election)), 
                  ncol = length(unique(polls1985_2024_long$party)), 0)


# stan input data
stan_dat <- list(
  
  N = nrow(polls1985_2024_long),
  K = length(unique(polls1985_2024_long$party)),
  R = length(unique(polls1985_2024_long$election)),
  KR = nrow(unique(polls1985_2024_long[,c("party", "election")])),
  P = length(unique(polls1985_2024_long$poll_id)),
  ST = sum(election_data$st_id),
  
  poll = polls1985_2024_long$poll_share,
  prev_vote = election_party_data$prev_vote_share,
  n = polls1985_2024$sample_size,
  max_D = as.numeric(max(polls1985_2024_long$dte)),
  
  k_id = polls1985_2024_long$k_id,
  p_id = polls1985_2024_long$p_id,
  r_id = polls1985_2024$r_id,
  t_id = polls1985_2024$t_id,
  st_id = election_data$st_id,
  st_r_id = election_data$st_r_id,
  
  zero_r = zero_r,
  zero_rK = zero_rK
  
)

# check model input
sapply(stan_dat, length)
sapply(stan_dat, range)

save(polls1985_2024_long, election_data, stan_dat, 
     file = paste0("forecast_data_", sc, ".RData"))


