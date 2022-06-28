# Author: Kevin See
# Purpose: prep data for Coho in the Newaukum
# Created: 3/18/2022
# Last Modified: 4/26/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(janitor)
library(mgcv)
library(gratia)

#------------------------------------------------------------
# read in data

# read in new redds
coho_new_redds <- read_excel("2021CohoNewRedds_25Apr2022.xlsx",
                             sheet = 1) %>%
  mutate(across(-c(basin, reach),
                as.numeric)) %>%
  pivot_longer(cols = -c(basin, reach),
               names_to = "stat_wk",
               values_to = "new_redds") %>%
  mutate(across(stat_wk,
                ordered,
                levels = c(30:52,
                           "53-1",
                           1:10))) %>%
  mutate(across(reach,
                as_factor)) %>%
  arrange(reach, stat_wk) %>%
  group_by(reach) %>%
  mutate(first_wk = min(stat_wk[!is.na(new_redds)], na.rm = T),
         last_wk = max(stat_wk[!is.na(new_redds)], na.rm = T)) %>%
  mutate(in_surv_season = between(stat_wk, unique(first_wk), unique(last_wk))) %>%
  ungroup()

# read in visible redds
coho_vis_redds <- read_excel("VisibleCohoRedds25Apr2022.xlsx",
                             sheet = 1) %>%
  clean_names() %>%
  mutate(across(-c(basin:survey_type),
                as.numeric)) %>%
  pivot_longer(cols = starts_with("x"),
               names_to = "stat_wk",
               values_to = "vis_redds") %>%
  mutate(across(stat_wk,
                str_remove,
                "^x")) %>%
  mutate(across(stat_wk,
                recode,
                "53_1" = "53-1")) %>%
  mutate(across(stat_wk,
                ordered,
                levels = c(30:52,
                           "53-1",
                           1:10))) %>%
  mutate(across(reach,
                as_factor)) %>%
  arrange(reach, stat_wk) %>%
  group_by(reach) %>%
  mutate(first_wk = min(stat_wk[!is.na(vis_redds)], na.rm = T),
         last_wk = max(stat_wk[!is.na(vis_redds)], na.rm = T)) %>%
  mutate(in_surv_season = between(stat_wk, unique(first_wk), unique(last_wk))) %>%
  ungroup()

# combine new and visible redds
coho_redds <- full_join(coho_vis_redds,
                        coho_new_redds) %>%
  select(basin:survey_type,
         first_wk, last_wk,
         stat_wk,
         in_surv_season,
         new_redds, vis_redds,
         everything()) %>%
  arrange(reach, stat_wk) %>%
  filter(in_surv_season) %>%
  group_by(basin, reach, survey_type) %>%
  mutate(cum_redds = cumsum(replace_na(new_redds, 0))) %>%
  ungroup()

# summary stats for each reach
rch_summ <- coho_redds %>%
  group_by(reach,
           first_wk,
           last_wk) %>%
  summarize(n_surv_poss = sum(in_surv_season),
            n_surv = sum(!is.na(new_redds)),
            n_non_zero = sum(new_redds > 0,
                             na.rm = T),
            n_zero = sum(new_redds == 0,
                         na.rm = T),
            n_missing = sum(in_surv_season & is.na(new_redds)),
            .groups = "drop")

rch_summ %>%
  # arrange(n_non_zero)
  # arrange(desc(n_non_zero))
  arrange(desc(n_missing),
          desc(n_non_zero))

rch_summ %>%
  # arrange(n_non_zero)
  # arrange(desc(n_non_zero))
  arrange(desc(n_missing),
          desc(n_non_zero)) %>%
  slice(5) %>%
  select(reach) %>%
  left_join(coho_redds) %>%
  filter(in_surv_season) %>%
  ggplot(aes(x = as.numeric(stat_wk),
             y = new_redds)) +
  geom_point() +
  geom_line() +
  labs(x = "Stat Week",
       y = "New Redds")


set.seed(5)
coho_redds %>%
  filter(in_surv_season) %>%
  filter(reach %in% sample(levels(coho_redds$reach), 6)) %>%
  ggplot(aes(x = as.numeric(stat_wk),
             y = new_redds)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ reach,
             scales = "free_y")

#---------------------------------------------------------------
# prepare a model dataset
mod_df <- rch_summ %>%
  filter(n_non_zero > 0) %>%
  left_join(coho_redds) %>%
  mutate(across(reach,
                fct_drop)) %>%
  # rename(redds = new_redds) %>%
  # select(-vis_redds, -cum_redds) %>%
  rename(redds = vis_redds) %>%
  select(-new_redds, -cum_redds) %>%
  filter(!is.na(redds)) %>%
  mutate(surv_day = as.numeric(stat_wk))

# fit a bunch of models

# single common smoother for all observations (G)
G <- gam(redds ~ s(surv_day, k = 5, m = 2, bs = "tp") +
           s(reach, k = nlevels(mod_df$reach), bs = "re"),
         data = mod_df,
         family = nb(theta = NULL,
                     link = "log"),
         method = "REML")
summary(G)
draw(G)

# global smoother with group level smoothers with the same wiggliness (GS)
GS <- gam(redds ~ s(surv_day, k = 5, m = 2, bs = "tp") +
            s(surv_day, reach, k = 5, bs = "fs", m = 2),
          data = mod_df,
          family = nb(theta = NULL,
                      link = "log"),
          method = "REML")

draw(GS)

# global smoother with group level smoothers with different wiggliness (GI)
GI <- gam(redds ~ s(surv_day, k = 5, m = 2, bs = "tp") +
            s(surv_day, by = reach, k = 5, m = 1, bs = "tp") +
            s(reach, bs = "re", k = nlevels(mod_df$reach)),
          data = mod_df,
          family = nb(theta = NULL,
                      link = "log"),
          method = "REML")

draw(GI)

# Group specific smoothers, no global smoother. All have the same wiggliness (S)
S <- gam(redds ~ s(surv_day, reach, k = 5, m = 2, bs = "fs"),
         data = mod_df,
         family = nb(theta = NULL,
                     link = "log"),
         method = "REML")
draw(S)

# Group specific smoothers, no global smoother. Different wiggliness (I)
I <- gam(redds ~ s(surv_day, by = reach, k = 5, m = 1, bs = "tp") +
           s(reach, bs = "re", k = nlevels(mod_df$reach)),
         data = mod_df,
         family = nb(theta = NULL,
                     link = "log"),
         method = "REML")

draw(I)


# build a model selection type table
tibble(mod_nm = c("G", "GS", 'GI', 'S', "I"),
       model = list(G, GS, GI, S, I)) %>%
  mutate(AIC = map_dbl(model,
                       AIC),
         AICc = map_dbl(model,
                        MuMIn::AICc)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         delta_AICc = AICc - min(AICc)) %>%
  arrange(AICc)
