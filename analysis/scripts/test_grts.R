# Author: Kevin See
# Purpose: test some GRTS draws on coho redd data
# Created: 6/24/22
# Last Modified: 6/27/2022
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

theme_set(theme_bw())

#-----------------------------------------------------------------
# needed functions
source("grts_functions.R")

#-----------------------------------------------------------------
# read in frame
strm_sf = st_read("data/GIS/Newaukum_CohoGRTS_SurveyFrame.shp")

strm_sf %>%
  ggplot() +
  geom_sf(color = "blue")

# read in all redds
redd_sf <- st_read("data/GIS/CohoRedds_2019.shp") %>%
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
  rbind(st_read("data/GIS/CohoRedds_2020.shp") %>%
          add_column(Year = 2020,
                     .before = 0) %>%
          select(any_of(names(redd_sf))) %>%
          st_transform(st_crs(strm_sf)))

# drop one redd point WAY outside stream layer
strm_bb <- st_bbox(strm_sf) %>%
  st_as_sfc()

redd_sf %<>%
  st_intersection(strm_bb)


#------------------------------------------------
# quick plot of points
ggplot() +
  geom_sf(data = strm_sf,
          color = "blue") +
  geom_sf(data = redd_sf,
          aes(color = as.factor(Year))) +
  labs(color = "Year") +
  theme(axis.text = element_blank())


#------------------------------------------------
# GRTS draw
# how long is the whole stream layer?
sum(st_length(strm_sf))

# convert meters to miles
sum(st_length(strm_sf)) |>
  as.numeric() |>
  measurements::conv_unit("ft", "mi")

# how many points?
143 * 0.3
n_pts = 45

#-------------------------------------------
# draw the points
# use set.seed() to make it consistent

set.seed(5)
grts_draw <- grts(strm_sf,
                  n_base = n_pts,
                  mindis = 4000,
                  maxtry = 50)

# quick plot of points
sp_plot(grts_draw, strm_sf)

# pull out points object
grts_pts <- grts_draw$sites_base

# how long should the reach extend on either side of each point?
# half-mile
lngth_req = 0.5

# set buffer size around that point
buff_size = 1

grts_rchs = create_grts_reaches(grts_pts,
                                strm_sf,
                                length_buffer_mi = lngth_req,
                                pt_buffer_mi = buff_size)

xtabs( ~ rch_lngth_mi <= 0.8, grts_rchs)

#--------------------------------------------------------
# any GRTS reaches overlap each other?
grts_overlap <- grts_rchs %>%
  st_overlaps() %>%
  enframe(name = "rch_1",
          value = "rch_2") %>%
  unnest(rch_2) %>%
  rowwise() %>%
  mutate(id = paste(min(rch_1, rch_2),
                    max(rch_1, rch_2),
                    sep = ".")) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-id)

grts_overlap




i = sample(1:nrow(grts_overlap), 1)
grts_rchs %>%
  slice(grts_overlap %>%
          slice(i) %>%
          as.numeric()) %>%
  ggplot() +
  geom_sf(aes(color = siteID,
              fill = siteID),
          size = 4,
          alpha = 0.5) +
  scale_color_brewer(palette = "Set1")


# remove the overlap
grts_rchs %>%
  st_difference() %>%
  rename(old_rch_lngth = rch_lngth_mi) %>%
  mutate(rch_lngth_mi = st_length(geometry),
         across(rch_lngth_mi,
                conv_unit,
                "ft", "mi"),
         across(rch_lngth_mi,
                as.numeric)) %>%
  mutate(diff = old_rch_lngth - rch_lngth_mi) %>%
  st_drop_geometry() %>%
  select(siteID, contains("lngth"), diff) %>%
  arrange(desc(diff))

grts_rchs %<>%
  st_difference() %>%
  mutate(rch_lngth_mi = st_length(geometry),
         across(rch_lngth_mi,
                conv_unit,
                "ft", "mi"),
         across(rch_lngth_mi,
                as.numeric))





ggplot() +
  geom_sf(data = strm_sf,
          color = "blue") +
  geom_sf(data = grts_pts,
          color = "red",
          size = 3,
          shape = 1) +
  geom_sf(data = grts_rchs,
          color = "firebrick",
          size = 1) +
  # geom_sf_text(data = grts_pts,
  #              aes(label = siteID),
  #              nudge_x = 5000,
  #              nudge_y = 1000) +
  theme(axis.text = element_blank())


#---------------------------------------------------------
# create small buffer around each GRTS reach
grts_buff <- grts_rchs %>%
  st_buffer(10)

grts_redds <- redd_sf %>%
  filter(Year == 2019) %>%
  st_intersection(grts_buff %>%
                    select(siteID, rch_lngth_mi))


# redd_sf %>%
#   filter(Year == 2019) %>%
#   mutate(redd_id = 1:n()) %>%
#   st_distance(grts_rchs) %>%
#   as_tibble() %>%
#   mutate(redd_id = 1:n()) %>%
#   pivot_longer(starts_with("V"),
#                names_to = "siteID",
#                values_to = "dist") %>%
#   group_by(redd_id) %>%
#   filter(dist == min(dist)) %>%
#   ungroup() %>%
#   mutate(siteID = paste0("Site-", str_extract(siteID, "[:digit:]+"))) %>%
#   filter(dist <= units::set_units(10, "ft")) %>%
#   left_join(redd_sf %>%
#               filter(Year == 2019) %>%
#               mutate(redd_id = 1:n())) %>%
#   st_sf()

  # filter(dist > units::set_units(6, "ft"),
  #        dist < units::set_units(100, "ft")) %>%
  # arrange(dist)
  # qplot(dist, data = .)


ggplot() +
  theme(axis.text = element_blank()) +
  geom_sf(data = strm_sf,
          color = "lightblue") +
  geom_sf(data = grts_buff,
          # size = 1,
          color = 'purple',
          fill = "purple",
          alpha = 0.5) +
  geom_sf(data = grts_redds,
          color = 'red',
          size = 1) +
  geom_sf(data = redd_sf %>%
            filter(Year == 2019) %>%
            anti_join(grts_redds %>%
                        st_drop_geometry() %>%
                        select(Year:NEAR_Y)),
          shape = 3,
          color = 'black')


mod_df <- grts_pts %>%
  select(siteID,
         ends_with("WGS84"),
         stratum, wgt) %>%
  mutate(pop = "Newaukum") %>%
  right_join(tidyr::expand(grts_rchs %>%
                             st_drop_geometry(),
                           nesting(siteID, rch_lngth_mi))) %>%
  left_join(grts_redds %>%
              st_drop_geometry() %>%
              as_tibble() %>%
              group_by(Year, siteID, rch_lngth_mi) %>%
              summarize(n_redds = sum(!is.na(redd_longi)),
                        .groups = "drop")) %>%
  mutate(across(Year,
                replace_na,
                2019),
         across(n_redds,
                replace_na,
                0)) %>%
  arrange(Year, siteID) %>%
  mutate(rch_lngth_ft = conv_unit(rch_lngth_mi,
                                  from = "mi",
                                  to = "ft")) %>%
  mutate(redd_dens = n_redds / rch_lngth_ft) %>%
  mutate(across(wgt,
                as.numeric))

# do the weights add up to the total length of the frame?
identical(sum(mod_df$wgt),
          as.numeric(sum(st_length(strm_sf))))

est_list <- spsurvey::cont_analysis(dframe = mod_df,
                                    vars = "redd_dens",
                                    subpops = "pop",
                                    siteID = "siteID",
                                    weight = "wgt")

est_list$Total

redd_sf %>%
  filter(Year == 2019) %>%
  st_drop_geometry() %>%
  summarize(n_redds = n())
