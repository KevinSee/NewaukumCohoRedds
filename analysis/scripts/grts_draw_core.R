# Author: Kevin See
# Purpose: Draw GRTS sites in non-core areas of the Newaukum for Coho redd surveys
# Created: 7/14/22
# Last Modified: 7/14/2022
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
                          "Newaukum_CohoGRTS_SurveyFrame.shp"),
                     quiet = T)

# read in core areas
core_sf = st_read(here("analysis/data/raw_data/GIS",
                       "CoreNewaukumCohoSurveys.shp"),
                  quiet = T) %>%
  st_transform(st_crs(coho_frame))

# read in non-core areas
noncore_sf = st_read(here("analysis/data/raw_data/GIS",
                          "GRTS_NewaukumLayerNoCore.shp"),
                     quiet = T) %>%
  st_transform(st_crs(coho_frame))

#-----------------------------------------------------------------
# set number of desired reaches
n_pts = 25

# draw the points
# set the seed to ensure reproducibility
# for a different set of sites, change the seed
set.seed(5)
grts_draw <- grts(noncore_sf,
                  n_base = n_pts,
                  n_over = 2*n_pts,
                  mindis = 4000,
                  maxtry = 50)

names(grts_draw)

# pull out GRTS points object
grts_pts <- grts_draw$sites_base

# grab the replacement points object
grts_replace <- grts_draw$sites_over

# generate 1 mile reaches near each GRTS point
grts_rchs = createGRTSreaches(grts_pts,
                              noncore_sf,
                              length_buffer_mi = 1)


# save some objects
st_write(grts_pts,
         dsn = here("analysis/data/derived_data",
                    "Newaukum_Coho_NonCore_GRTS_base.gpkg"))
st_write(grts_replace,
         dsn = here("analysis/data/derived_data",
                    "Newaukum_Coho_NonCore_GRTS_replace.gpkg"))
st_write(grts_rchs,
         dsn = here("analysis/data/derived_data",
                    "Newaukum_Coho_NonCore_GRTS_reach_base.gpkg"))


# quick plot showing sites
ggplot() +
  geom_sf(data = core_sf,
          aes(color = "Core")) +
  geom_sf(data = noncore_sf,
          aes(color = "Non-Core")) +
  geom_sf(data = grts_rchs,
          # size = 1,
          aes(color = 'GRTS')) +
  geom_sf(data = grts_pts,
          shape = 1,
          size = 2,
          aes(color = "GRTS")) +
  scale_color_manual(name = "Classificiation",
                     values = c("Core" = "darkblue",
                                "Non-Core" = "lightblue",
                                "GRTS" = "red"))

#-----------------------------------------------------------------
# if need a replacement point/reach, use this to generate 1 mile reach
site_id = "Site-26"
grts_rch_replace <- grts_replace %>%
  filter(siteID == site_id) %>%
  createGRTSreaches(noncore_sf,
                    length_buffer_mi = 1)
