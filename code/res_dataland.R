#-------------------------------------------------------------------------------
#
# Estimates mean vote shares and win probabilities elections 1985-2024
#
# Author: Sina Chen
#
#-------------------------------------------------------------------------------


# Libraries ---------------------------------------------------------------

{
  library(dplyr)
  library(tidyr)
  library(rv)
  library(stringr)
  library(readr)
}



# Data --------------------------------------------------------------------

setwd("your working directory")

# scenario: enter A to F
sc <- "A"

# poll data
load(paste0("~/economist/forecast_data_", sc, ".RData"))

# filter 2024 elections
election_data <- election_data %>% 
  filter(str_detect(election, "2024"))

# filter 2024 polls
polls2024 <- polls1985_2024_long %>% 
  filter(year == 2024)
     
# posterior simulations for estimated party support
theta <- readRDS(paste0("theta_", sc, ".rds"))

# demographics
demographics <- read_csv("dataland/dataland_demographics.csv")



# Preparation -------------------------------------------------------------

# summary of estimated support
theta_summary <- lapply(1:dim(theta)[1], function(x) summary(theta[x,,]))

# add party
theta_summary <- lapply(theta_summary, function(x) mutate(
  x, party = rep(c("cc", "dgm", "pdal", "ssp"), each = stan_dat$max_D+1)))

# combine single elections
theta_election_df <- do.call(rbind, theta_summary)

# add election, day of campaign (d) and Synapse Territories id (st_id)
theta_election_df <- theta_election_df %>% 
  mutate(election = rep(election_data$election, 
                        each = (stan_dat$max_D+1) * stan_dat$K),
         d = rep(1:(stan_dat$max_D+1), stan_dat$K*length(theta_summary)),
         st_id =  rep(election_data$st_id, 
                      each=(stan_dat$max_D+1) * stan_dat$K))

# reshape to wide
theta_election_df_wide <- theta_election_df %>% 
  select("election", "party", "mean", "d", "st_id") %>% 
  pivot_wider(names_from = party, values_from = mean)



# Mean vote share & win probability ---------------------------------------

#### Province level ####

# mean vote share
province_forecast <- theta_election_df_wide %>% 
  filter(str_detect(election, "2024")) %>% 
  mutate(date = rep(seq(from = as.Date("2024-04-02", "%Y-%m-%d"), 
                        to = as.Date("2024-06-16", "%Y-%m-%d"),
                        by = "day"), nrow(election_data)),
         province = str_remove_all(election, "\\d+")) %>% 
  rename("cc_mean_vote_share" = cc,
         "dgm_mean_vote_share" = dgm,
         "pdal_mean_vote_share" = pdal,
         "ssp_mean_vote_share" = ssp)

# election win probability
province_win_prob <- list()
for(r in 1:nrow(election_data)){           # loop through elections
  election_sim <- theta[r,,]               # subset election
  win_prob <- data.frame(cc = c(), dgm = c(), pdal = c(), ssp = c())
  for(d in 1:nrow(election_sim)) {         # loop through days of campaign
    day <- election_sim[d,]                # subset day
    day_sim <- day %>% unclass()%>% as.data.frame()
    colnames(day_sim) <- c("cc", "dgm", "pdal", "ssp")
    winner <- t(apply(day_sim, 1, function(x) ifelse(x == max(x), 1, 0)))   # identify winning party 
    win_prob <- rbind(win_prob, apply(winner, 2, function(x)                # compute win probability as simulation share winning
      sum(x)/nrow(winner)))
    colnames(win_prob) <- c("cc_win_probability", "dgm_win_probability", 
                            "pdal_win_probability", "ssp_win_probability")
  }
  province_win_prob[[(r)]] <- win_prob
}
province_win_prop_df <- do.call(rbind, province_win_prob)
rm(province_win_prob, election_sim, day, day_sim, winner, win_prob, r, d)

# add date and election
province_win_prop_df <- province_win_prop_df %>% 
  mutate(date = rep(seq(from = as.Date("2024-04-02", "%Y-%m-%d"), 
                        to = as.Date("2024-06-16", "%Y-%m-%d"),
                        by = "day"), nrow(election_data)),
         election = rep(election_data$election, each = stan_dat$max_D+1))

# merge mean vote share with win probability
province_df <- merge(province_forecast, province_win_prop_df, 
                     by = c("election", "date")) %>% 
  select(date, province, cc_mean_vote_share, cc_win_probability, 
         dgm_mean_vote_share, dgm_win_probability, pdal_mean_vote_share, 
         pdal_win_probability, ssp_mean_vote_share, ssp_win_probability)
write.csv(province_df, paste0("provincial_forecast_", sc, ".csv"))
rm(province_win_prop_df, province_df)


#### National level ####

# compute population share for each province for national popular vote 
demographics <- demographics %>% 
  mutate(pop_share = population/sum(population))

# add population share to province forecast
province_forecast <- merge(province_forecast, demographics %>% 
                             select(province, pop_share), 
                           by = "province", all.x = T)

# aggregate to national forecast
national_forecast <- province_forecast %>% 
  group_by(date) %>% 
  summarise(cc_mean_vote_share = sum(cc_mean_vote_share*pop_share),
            dgm_mean_vote_share = sum(dgm_mean_vote_share*pop_share),
            pdal_mean_vote_share = sum(pdal_mean_vote_share*pop_share),
            ssp_mean_vote_share = sum(ssp_mean_vote_share*pop_share))

# electoral votes (ev) for winning party over every simulation
ev_sim_df <-list()
for(r in 1: nrow(election_data)){          # loop through elections
  election_sim <- theta[r,,]               # subset election
  ev_sim <- list()
  for(i in 1:nrow(election_sim)) {         # loop through days of campaign
    day <- election_sim[i,]                # subset day
    day_sim <- day %>% unclass() %>% 
      as.data.frame()
    colnames(day_sim) <- c("cc", "dgm", "pdal", "ssp")
    winner <- t(apply(day_sim, 1, function(x) 
      ifelse(x == max(x), demographics$electoral_college_votes[r], 0))) # assign ev to winning party
    ev_sim[[i]] <- winner
  }
  ev_sim_df[[r]] <- ev_sim
}

# compute national win probability for every day based on simulations
nat_win_day_df <- data.frame()
for(d in 1:76){                            # loop through campaign days
  day_sim_ev <- data.frame()
  for(s in 1:750){                         # loop through simulations
    election_day_sim_df <- data.frame()
    for(r in 1: nrow(election_data)){      # loop through elections
      election_day_sim_df <- rbind(election_day_sim_df, ev_sim_df[[r]][[d]][s,])
    }
    sim_ev <- apply(election_day_sim_df, 2, function(x) sum(x)) %>% unname()    # sum ev single simulation
    day_sim_ev <- rbind(day_sim_ev, sim_ev)                                     # sum ev all simulations 
  }
  
  winner_ev <- t(apply(day_sim_ev, 1, function(x) ifelse(x == max(x), 1, 0)))   # identify national winner
  win_ev_prob <- apply(winner_ev, 2, function(x) sum(x)/nrow(winner_ev)) %>%    # compute win probability based on simulations
    unname()                
  
  nat_win_day_df <- rbind(nat_win_day_df, win_ev_prob)
}
rm(winner, day_sim_ev, election_day_sim_df, sim_ev, winner_ev, win_ev_prob, d, 
   s, r, i)
colnames(nat_win_day_df) <- c("cc_win_probability", "dgm_win_probability", 
                               "pdal_win_probability", "ssp_win_probability")

# cbind national mean vote share with win probability
national_df <- cbind(national_forecast, nat_win_day_df) %>% 
  select(date, cc_mean_vote_share, cc_win_probability, dgm_mean_vote_share, 
          dgm_win_probability, pdal_mean_vote_share, pdal_win_probability, 
          ssp_mean_vote_share, ssp_win_probability)
write.csv(national_df, paste0("national_forecast_", sc, ".csv"))
















