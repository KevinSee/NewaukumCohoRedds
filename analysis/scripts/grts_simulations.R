# Author: Kevin See
# Purpose: test some GRTS draws on coho redd data
# Created: 6/24/22
# Last Modified: 6/29/2022
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
                           "IndexOnlyCohoSurveyFrame.shp"))

ggplot() +
  geom_sf(data = coho_frame,
          color = "blue") +
  geom_sf(data = index_frame,
          color = "red")

# focus on the index frame only
strm_sf = index_frame

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
  st_transform(st_crs(strm_sf))

redd_sf %<>%
  rbind(st_read(here("analysis/data/raw_data/GIS",
                     "CohoRedds_2020.shp")) %>%
          add_column(Year = 2020,
                     .before = 0) %>%
          select(any_of(names(redd_sf))) %>%
          st_transform(st_crs(strm_sf)))

redd_2021 <- st_read(here("analysis/data/raw_data/GIS",
                          "CohoRedds_2021.shp")) %>%
  add_column(Year = 2021,
             .before = 0)
names(redd_2021)[which(names(redd_2021) == "reach")] = "Reach_Name"
names(redd_2021)[which(names(redd_2021) == "redd_name")] = "sgs_redd_n"
redd_2021 %<>%
  select(any_of(names(redd_sf))) %>%
  st_transform(st_crs(strm_sf))

redd_sf %<>%
  rbind(redd_2021)

rm(redd_2021)
# rm(coho_frame,
#    index_frame)

# drop a couple redd points WAY outside stream layer
# strm_bb <- st_bbox(strm_sf) %>%
#   st_as_sfc()
#
# redd_sf %<>%
#   st_intersection(strm_bb)

strm_buff <- st_buffer(strm_sf,
                       30)

other_redds <- redd_sf %>%
  anti_join(redd_sf %>%
              st_intersection(strm_buff) %>%
              select(sgs_redd_n) %>%
              st_drop_geometry()) %>%
  filter(NEAR_X != -1,
         NEAR_Y != -1)

redd_sf %<>%
  st_intersection(strm_buff)


# pull out total number of redds observed each year
obs_redds <- redd_sf %>%
  st_drop_geometry() %>%
  group_by(Year) %>%
  summarize(obs = n(),
            .groups = "drop")

#------------------------------------------------
# quick plot of points
ggplot() +
  geom_sf(data = coho_frame,
          color = "lightblue") +
  geom_sf(data = strm_sf,
          color = "blue") +
  geom_sf(data = other_redds,
          aes(color = as.factor(Year)),
          shape = 2) +
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
# how long is the whole stream layer?
sum(st_length(strm_sf))

# convert meters to miles
sum(st_length(strm_sf)) |>
  as.numeric() |>
  measurements::conv_unit("ft", "mi")

frame_lngth <- sum(st_length(strm_sf)) |>
  as.numeric() |>
  measurements::conv_unit("ft", "mi") |>
  round()


# how many points?
coverage_df <- tibble(perc_cov = seq(.2, .5, by = .1)) %>%
  mutate(n_miles = frame_lngth * perc_cov) %>%
  # mutate(n_sites = plyr::round_any(n_miles, 5))
  mutate(n_sites = plyr::round_any(n_miles, 1))
coverage_df


# # aim for 30% coverage
# n_pts = 45

# loop over different numbers of points
for(n_pts in coverage_df$n_sites) {

  cat(paste("\n\t Starting with", n_pts, "sites \n\n\n"))

  #-------------------------------------------
  # set up loop of simulations
  # simulate GRTS points, test what estimates
  # would be if those points were used in
  # 2019 and 2020

  # how long should the reach extend on either side of each point?
  # half-mile
  lngth_req = 0.5

  # set buffer size around that point (in miles)
  buff_size = 1


  # how many simulations should we do?
  n_sim = 50

  sim_res = NULL
  set.seed(6)
  for(s in 1:n_sim) {
    cat(paste("Starting sim", s, "with", n_pts, "points\n"))

    # draw the points
    grts_draw <- grts(strm_sf,
                      n_base = n_pts,
                      mindis = 4000,
                      maxtry = 50)

    # pull out GRTS points object
    grts_pts <- grts_draw$sites_base

    cat(paste("Creating GRTS reaches\n\n"))

    # generate the reaches around each point
    grts_rchs = createGRTSreaches(grts_pts,
                                  strm_sf,
                                  length_buffer_mi = lngth_req)

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
      st_intersection(grts_buff %>%
                        select(siteID, rch_lngth_mi))

    mod_df <- grts_pts %>%
      select(siteID,
             ends_with("WGS84"),
             stratum, wgt) %>%
      mutate(pop = "Newaukum") %>%
      right_join(tidyr::expand(grts_rchs %>%
                                 st_drop_geometry(),
                               nesting(siteID, rch_lngth_mi),
                               tidyr::crossing(Year = unique(redd_sf$Year)))) %>%
      left_join(grts_redds %>%
                  st_drop_geometry() %>%
                  as_tibble() %>%
                  group_by(Year, siteID, rch_lngth_mi) %>%
                  summarize(n_redds = sum(!is.na(redd_longi)),
                            .groups = "drop")) %>%
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
      left_join(obs_redds) %>%
      add_column(sim = s,
                 .before = 0)

    if(s == 1) {
      sim_res = est_df
    } else {
      sim_res <- sim_res %>%
        bind_rows(est_df)
    }
  }

  # save results
  write_rds(sim_res,
            file = here("analysis/data/derived_data",
                        paste0("grts_sims_", n_pts, ".rds")))
}

# get all simulation results
sim_res <- coverage_df %>%
  select(perc_cov,
         n_sites) %>%
  mutate(sim_pts = map(n_sites,
                       .f = function(n_pts) {
                         read_rds(here("analysis/data/derived_data",
                                       paste0("grts_sims_", n_pts, ".rds")))
                       })) %>%
  unnest(sim_pts)

write_rds(sim_res,
          file = here("analysis/data/derived_data",
                      "grts_sims_all.rds"))

# write to a csv sheet
write_csv(sim_res,
          file = here("analysis/data/derived_data",
                      "Newaukum_GRTS_sims.csv"))
#------------------------------------------------------------
sim_res <- read_rds(here("analysis/data/derived_data",
                         "grts_sims_all.rds"))

save(strm_sf,
     redd_sf,
     obs_redds,
     sim_res,
     file = here("analysis/data/derived_data",
                 "grts_data_sims.rda"))


#------------------------------------------------------------
load(here("analysis/data/derived_data",
          "grts_data_sims.rda"))

# some quick analyses on the results
sim_res %>%
  rowwise() %>%
  mutate(in_ci = between(obs, LCB95Pct, UCB95Pct)) %>%
  ungroup() %>%
  janitor::tabyl(perc_cov, in_ci) %>%
  janitor::adorn_percentages() %>%
  janitor::adorn_pct_formatting()

dw = 0.8
sim_res %>%
  mutate(ci_width = UCB95Pct - LCB95Pct,
         est_cv = StdError / Estimate) %>%
  mutate(across(c(Year,
                  perc_cov),
                as_factor)) %>%
  ggplot(aes(x = Year,
             y = est_cv,
             fill = perc_cov)) +
  geom_boxplot(position = position_dodge(dw)) +
  geom_hline(yintercept = 0.15,
             linetype = 2) +
  scale_fill_viridis_d(name = "% Surveyed") +
  labs(y = "CV of Estimate")

ggsave(filename = here("analysis/figures",
                       "CV_Newaukum_GRTS.png"),
       width = 6,
       height = 6)



dw = 0.8
sim_res %>%
  mutate(across(c(Year,
                  perc_cov),
                as_factor)) %>%
  ggplot(aes(x = Year,
             y = Estimate,
             fill = perc_cov)) +
  geom_boxplot(position = position_dodge(dw)) +
  geom_point(aes(y = obs),
             position = position_dodge(dw),
             size = 4) +
  scale_fill_viridis_d(name = "% Surveyed")

sim_res %>%
  mutate(across(c(Year,
                  perc_cov),
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
                  Year),
                as_factor)) %>%
  mutate(bias = Estimate - obs,
         rel_bias = bias / obs) %>%
  ggplot(aes(x = Year,
             y = rel_bias * 100,
             fill = perc_cov)) +
  geom_boxplot() +
  geom_hline(yintercept = 0,
             linetype = 2) +
  labs(y = "Relative Bias (%)") +
  scale_fill_viridis_d(name = "% Surveyed")

ggsave(filename = here("analysis/figures",
                       "Relative_Bias_Newaukum_GRTS.png"),
       width = 6,
       height = 6)



sim_res %>%
  mutate(across(c(perc_cov,
                  Year),
                as_factor)) %>%
  mutate(ci_width = UCB95Pct - LCB95Pct,
         est_cv = StdError / Estimate,
         bias = Estimate - obs,
         rel_bias = bias / obs) %>%
  group_by(Year,
           perc_cov) %>%
  summarize(MB = median(bias),
            MAE = mean(abs(bias), na.rm = T),
            MAPE = mean(abs(rel_bias) * 100, na.rm = T),
            MSA = exp(median(abs(log(Estimate/obs)), na.rm = T)) * 100,
            RMSE = sqrt(mean(bias^2)),
            MCI = median(ci_width / Estimate),
            MCV = median(est_cv),
            .groups = "drop")
