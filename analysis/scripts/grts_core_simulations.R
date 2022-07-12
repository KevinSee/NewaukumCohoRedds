# Author: Kevin See
# Purpose: test some GRTS draws on coho redd data
# Created: 6/24/22
# Last Modified: 7/11/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(lubridate)
library(magrittr)
library(sf)
library(spsurvey)
library(lwgeom)
library(measurements)
library(sfnetworks)
library(units)
library(here)

theme_set(theme_bw())

#-----------------------------------------------------------------
# needed functions
devtools::load_all()

#-----------------------------------------------------------------
# read in frame
coho_frame = st_read(here("analysis/data/raw_data/GIS",
                          "Newaukum_CohoGRTS_SurveyFrame.shp"))

# index areas only
index_frame = st_read(here("analysis/data/raw_data/GIS",
                           "IndexOnlyCohoSurveyFrame.shp")) %>%
  st_transform(st_crs(coho_frame))

# read in core areas
core_sf = st_read(here("analysis/data/raw_data/GIS",
                       "CoreNewaukumCohoSurveys.shp")) %>%
  st_transform(st_crs(coho_frame))

# read in non-core areas
noncore_sf = st_read(here("analysis/data/raw_data/GIS",
                          "GRTS_NewaukumLayerNoCore.shp")) %>%
  st_transform(st_crs(coho_frame))


# clip the core and non-core areas to the index areas
core_sf %<>%
  st_intersection(index_frame %>%
                    st_buffer(10) %>%
                    select(geometry))
noncore_sf %<>%
  st_intersection(index_frame %>%
                    st_buffer(10) %>%
                    select(geometry))



ggplot() +
  geom_sf(data = coho_frame,
          aes(color = "Coho Frame")) +
  geom_sf(data = index_frame,
          aes(color = "Index Areas"),
          size = 2) +
  geom_sf(data = core_sf,
          aes(color = "Core")) +
  geom_sf(data = noncore_sf,
          aes(color = "GRTS")) +
  scale_color_manual(name = "Category",
                     values = c("Coho Frame" = "lightblue",
                                "Index Areas" = "gray",
                                "Core" = "red",
                                "GRTS" = "purple"))

# focus on the index frame only
strm_sf = core_sf %>%
  mutate(type = "core") %>%
  rbind(noncore_sf %>%
          mutate(type = "grts"))

# read in all redds
redd_sf <- st_read(here("analysis/data/raw_data/GIS",
                        "CohoRedds_2019.shp")) %>%
  rename(survey_dat = Date,
         Reach_Name = Reach,
         redd_latit = Latitude,
         redd_longi = Longitude,
         sgs_redd_n = Redd_name) %>%
  select(-Stream) %>%
  add_column(Year = 2019,
             .before = 0) %>%
  st_transform(st_crs(coho_frame))

redd_sf %<>%
  rbind(st_read(here("analysis/data/raw_data/GIS",
                     "CohoRedds_2020.shp")) %>%
          add_column(Year = 2020,
                     .before = 0) %>%
          select(any_of(names(redd_sf))) %>%
          st_transform(st_crs(coho_frame)))

redd_2021 <- st_read(here("analysis/data/raw_data/GIS",
                          "CohoRedds_2021.shp")) %>%
  add_column(Year = 2021,
             .before = 0)
names(redd_2021)[which(names(redd_2021) == "reach")] = "Reach_Name"
names(redd_2021)[which(names(redd_2021) == "redd_name")] = "sgs_redd_n"
redd_2021 %<>%
  select(any_of(names(redd_sf))) %>%
  st_transform(st_crs(coho_frame))

redd_sf %<>%
  rbind(redd_2021) %>%
  mutate(id = 1:n()) %>%
  relocate(id)

rm(redd_2021)

# which redds were in core areas, and which in non-core areas?
redd_sf %<>%
  st_join(strm_sf %>%
            select(type) %>%
            st_buffer(5)) %>%
  distinct() %>%
  arrange(id,
          type) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

# drop redds too far away
other_redds <- redd_sf %>%
  filter(is.na(type))

redd_sf %<>%
  filter(!is.na(type))

# pull out total number of redds observed each year
obs_redds <- redd_sf %>%
  st_drop_geometry() %>%
  group_by(Year, type) %>%
  summarize(obs = n(),
            .groups = "drop")
obs_redds %<>%
  bind_rows(obs_redds %>%
              group_by(Year) %>%
              summarize(across(obs,
                               sum)) %>%
              mutate(type = "total")) %>%
  arrange(Year, type)

#------------------------------------------------
# quick plot of points
ggplot() +
  geom_sf(data = coho_frame,
          color = "lightblue") +
  geom_sf(data = strm_sf,
          color = "blue") +
  # geom_sf(data = other_redds,
  #         aes(color = as.factor(Year)),
  #         shape = 2) +
  geom_sf(data = redd_sf,
          aes(color = as.factor(Year))) +
  scale_color_brewer(palette = "Set1",
                     name = "Year") +
  theme(axis.text = element_blank()) +
  facet_wrap(~ Year,
             nrow = 2,
             ncol = 2)


#------------------------------------------------
# GRTS draw
# how long is the core & non-core stream layer?
coho_lngth <- strm_sf %>%
  mutate(lngth = st_length(geometry)) %>%
  group_by(type) %>%
  summarize(tot_lngth = sum(lngth),
            .groups = "drop") %>%
  st_drop_geometry() %>%
  mutate(across(tot_lngth,
                as.numeric),
         across(tot_lngth,
                measurements::conv_unit,
                "ft",
                "mi")) %>%
  rename(tot_lngth_mi = tot_lngth)


frame_lngth <- coho_lngth %>%
  filter(type == "grts") %>%
  pull(tot_lngth_mi) %>%
  round()

# how many points?
coverage_df <- tibble(n_sites = c(15, 20, 25, 30)) %>%
  mutate(n_miles = n_sites,
         perc_cov = n_miles / frame_lngth)

coverage_df


# # aim for 20 sites
# n_pts = 20

#-------------------------------------------
# set up loop of simulations
# simulate GRTS points, test what estimates
# would be if those points were used in
# 2019 - 2021

# how long should the reach extend on either side of each point?
# half-mile
lngth_req = 0.5

# how many simulations should we do?
n_sim = 100

# loop over different numbers of points
for(n_pts in coverage_df$n_sites[c(3,4)]) {

  cat(paste("\n\t Starting with", n_pts, "sites \n\n\n"))

  sim_res = NULL
  set.seed(6)
  for(s in 1:n_sim) {
    cat(paste("Starting sim", s, "with", n_pts, "points\n"))

    # draw the points
    grts_draw <- grts(noncore_sf,
                      n_base = n_pts,
                      mindis = 4000,
                      maxtry = 50)

    # pull out GRTS points object
    grts_pts <- grts_draw$sites_base

    cat("Creating GRTS reaches\n\n")

    # generate the reaches around each point
    grts_rchs = NULL
    while(class(grts_rchs)[1] != "sf") {
      grts_rchs = try(createGRTSreaches(grts_pts,
                                        noncore_sf,
                                        length_buffer_mi = lngth_req))
      if(class(grts_rchs)[1] == "try-error") {
        cat("Reach creation failed. Trying again.\n")
      }
    }

    # remove any overlap between reaches
    grts_rchs %<>%
      st_difference() %>%
      mutate(rch_lngth_mi = st_length(geometry),
             across(rch_lngth_mi,
                    conv_unit,
                    "ft", "mi"),
             across(rch_lngth_mi,
                    as.numeric))

    # create small buffer around each GRTS reach
    grts_buff <- grts_rchs %>%
      st_buffer(10)

    # determine which observed redds were found in each GRTS reach
    grts_redds <- redd_sf %>%
      filter(type == "grts") %>%
      st_intersection(grts_buff %>%
                        select(siteID, rch_lngth_mi))

    mod_df <- grts_pts %>%
      select(siteID,
             ends_with("WGS84"),
             stratum, wgt) %>%
      mutate(pop = "Newaukum") %>%
      right_join(grts_rchs %>%
                   st_drop_geometry() %>%
                   select(siteID, rch_lngth_mi) %>%
                   tidyr::crossing(Year = unique(redd_sf$Year)),
                 by = "siteID") %>%
      left_join(grts_redds %>%
                  st_drop_geometry() %>%
                  as_tibble() %>%
                  group_by(siteID, Year) %>%
                  summarize(n_redds = sum(!is.na(redd_longi)),
                            .groups = "drop"),
                by = c("siteID", "Year")) %>%
      mutate(across(n_redds,
                    replace_na,
                    0)) %>%
      arrange(Year, siteID) %>%
      mutate(rch_lngth_ft = conv_unit(rch_lngth_mi,
                                      from = "mi",
                                      to = "ft")) %>%
      mutate(redd_dens = n_redds / rch_lngth_ft) %>%
      mutate(across(wgt,
                    as.numeric))

    # # do the weights add up to the total length of the frame?
    # identical(sum(mod_df$wgt[mod_df$Year == 2019]),
    #           sum(mod_df$wgt[mod_df$Year == 2020]))
    # identical(sum(mod_df$wgt[mod_df$Year == 2019]),
    #           as.numeric(sum(st_length(strm_sf))))

    # generate estimates for each year
    est_df <- mod_df %>%
      group_by(Year) %>%
      nest() %>%
      mutate(spsurv_list = map(data,
                               .f = function(x) {
                                 cont_analysis(dframe = x,
                                               vars = "redd_dens",
                                               subpops = "pop",
                                               siteID = "siteID",
                                               weight = "wgt")
                               }),
             tot = map(spsurv_list,
                       "Total")) %>%
      ungroup() %>%
      unnest(tot) %>%
      left_join(obs_redds %>%
                  filter(type == "grts"),
                by = "Year") %>%
      add_column(sim = s,
                 .before = 0)

    if(s == 1) {
      sim_res = est_df
    } else {
      sim_res <- sim_res %>%
        bind_rows(est_df)
    }
    rm(grts_draw,
       grts_pts,
       grts_redds,
       grts_rchs,
       grts_buff,
       mod_df)
  }

  # save results
  write_rds(sim_res,
            file = here("analysis/data/derived_data",
                        paste0("grts_core_sims_", n_pts, ".rds")))
}

# get all simulation results
sim_res <- coverage_df %>%
  select(perc_cov,
         n_sites) %>%
  mutate(sim_pts = map(n_sites,
                       .f = function(n_pts) {
                         read_rds(here("analysis/data/derived_data",
                                       paste0("grts_core_sims_", n_pts, ".rds")))
                       })) %>%
  unnest(sim_pts)

write_rds(sim_res,
          file = here("analysis/data/derived_data",
                      "grts_core_sims_all.rds"))

# write to a csv sheet
write_csv(sim_res,
          file = here("analysis/data/derived_data",
                      "Newaukum_GRTS_core_sims.csv"))
#------------------------------------------------------------
sim_res <- read_rds(here("analysis/data/derived_data",
                         "grts_core_sims_all.rds"))

save(strm_sf,
     redd_sf,
     obs_redds,
     sim_res,
     file = here("analysis/data/derived_data",
                 "grts_data_core_sims.rda"))


#------------------------------------------------------------
load(here("analysis/data/derived_data",
          "grts_data_core_sims.rda"))

# some quick analyses on the results
sim_res %>%
  rowwise() %>%
  mutate(in_ci = between(obs, LCB95Pct, UCB95Pct)) %>%
  ungroup() %>%
  janitor::tabyl(perc_cov, in_ci) %>%
  janitor::adorn_percentages() %>%
  janitor::adorn_pct_formatting()

sim_res %<>%
  select(-type,
         -obs) %>%
  left_join(obs_redds %>%
              pivot_wider(names_from = type,
                          values_from = obs)) %>%
  mutate(across(c(Estimate,
                  ends_with("95Pct")),
                ~ . + core)) %>%
  rename(obs = total)


dw = 0.8
sim_res %>%
  mutate(ci_width = UCB95Pct - LCB95Pct,
         est_cv = StdError / Estimate) %>%
  mutate(across(c(Year,
                  perc_cov,
                  n_sites),
                as_factor)) %>%
  ggplot(aes(x = Year,
             y = est_cv,
             fill = n_sites)) +
  geom_boxplot(position = position_dodge(dw)) +
  geom_hline(yintercept = 0.15,
             linetype = 2) +
  scale_fill_viridis_d(name = "n Sites") +
  labs(y = "CV of Estimate")

ggsave(filename = here("analysis/figures",
                       "CV_Newaukum_GRTS.png"),
       width = 6,
       height = 6)



dw = 0.8
sim_res %>%
  mutate(across(c(Year,
                  perc_cov,
                  n_sites),
                as_factor)) %>%
  ggplot(aes(x = Year,
             y = Estimate,
             fill = n_sites)) +
  geom_boxplot(position = position_dodge(dw)) +
  geom_point(aes(y = obs),
             position = position_dodge(dw),
             size = 4) +
  scale_fill_viridis_d(name = "n Sites")

sim_res %>%
  mutate(across(c(Year,
                  perc_cov,
                  n_sites),
                as_factor)) %>%
  ggplot(aes(x = sim,
             y = Estimate,
             color = Year)) +
  geom_errorbar(aes(ymin = LCB95Pct,
                    ymax = UCB95Pct),
                width = 0) +
  geom_point() +
  geom_hline(aes(yintercept = obs),
             linetype = 2) +
  facet_wrap(~ Year,
             scales = "free_y")


sim_res %>%
  mutate(across(c(perc_cov,
                  Year,
                  n_sites),
                as_factor)) %>%
  mutate(bias = Estimate - obs,
         rel_bias = bias / obs) %>%
  ggplot(aes(x = Year,
             y = rel_bias * 100,
             fill = n_sites)) +
  geom_boxplot() +
  geom_hline(yintercept = 0,
             linetype = 2) +
  labs(y = "Relative Bias (%)") +
  scale_y_continuous(limits = c(-20, 30)) +
  scale_fill_viridis_d(name = "n Sites")

ggsave(filename = here("analysis/figures",
                       "Relative_Bias_Newaukum_GRTS.png"),
       width = 6,
       height = 6)



sim_res %>%
  mutate(across(c(perc_cov,
                  Year,
                  n_sites),
                as_factor)) %>%
  mutate(ci_width = UCB95Pct - LCB95Pct,
         est_cv = StdError / Estimate,
         bias = Estimate - obs,
         rel_bias = bias / obs) %>%
  group_by(Year,
           n_sites) %>%
  summarize(MB = median(bias),
            MAE = mean(abs(bias), na.rm = T),
            MAPE = median(abs(rel_bias) * 100, na.rm = T),
            MSA = exp(median(abs(log(Estimate/obs)), na.rm = T)) * 100,
            RMSE = sqrt(mean(bias^2)),
            MCI = median(ci_width / Estimate),
            MCV = median(est_cv),
            .groups = "drop")
