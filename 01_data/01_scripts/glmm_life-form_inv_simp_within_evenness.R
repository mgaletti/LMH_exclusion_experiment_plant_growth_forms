
# load packages and table ----------------------------------------------------------


library(tidyverse)
library(hillR)
library(hablar)
library(stringr)
library(lme4) 
library(car)
library(naniar)
library(ggpubr)
library(grid)
library(codyn)
library(rcompanion)
library(MASS)
library(lmtest)



data_biota <- read_csv("00_tables/life_form_2020_for_abudance.csv")
glimpse(data_biota)


# table with inverse simpsons calculated  -------------------------------------------------------------------


bt.invsim.mtz <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", "..")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  mutate(site=Site, treatment = Treatment, plot = Plot, month = Month, time = Time, lifeform = life_form) %>% 
  unite(Site, Site, Treatment, Plot, Month, Time, life_form) %>% rename(PlotID = Site) %>% 
  dplyr::group_by(PlotID, site, treatment, plot, month, time, Species) %>%
  dplyr::summarise(abundances = sum(value)) %>% 
  ungroup(PlotID, site, treatment, plot, month, time, Species) %>%
  spread(Species, abundances) %>% 
  replace(is.na(.), 0) %>% 
  remove_rownames %>% 
  dplyr::select(-c(2:6)) %>% 
  column_to_rownames(var="PlotID") %>% 
  hillR::hill_taxa(q = 2, MARGIN = 1) %>% 
  as.data.frame() %>% rownames_to_column(var="PlotID") %>% 
  separate(PlotID, c("Site", "Treatment", "Plot", "Month", "Time", "Life_form"), convert = TRUE) %>% 
  rename(simp.inv.div = ".") %>% 
  hablar::rationalize() %>% 
  spread(Life_form, simp.inv.div) %>%
  hablar::rationalize() %>% 
  mutate(Month = as.factor(Month)) %>%
  mutate(Treatment = as.factor(Treatment)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open", "closed"))) %>% 
  mutate(site = Site,
         treatment = Treatment,
         plot = Plot,
         month = Month,
         time = Time) %>%
  unite(site, site, treatment, plot, month, time) %>% rename(PlotID = site) %>% 
  rename(bamboo.invs = bamboo,
         herb.invs = herb,
         liana.invs = liana,
         palm.invs = palm,
         shrub.invs = shrub,
         tree.invs = tree)
bt.invsim.mtz



# total abundance table  -------------------------------------------------------------------


bt.abu.mtz <- data_biota %>% 
  dplyr::select(-c(4, 24:26,28)) %>%
  gather(key = "Month", value = "value", 4:22) %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::group_by(Site, Plot, Treatment, Month, Time, life_form) %>% 
  summarise(abundances = sum(value)) %>% 
  dplyr::mutate(site = Site, plot = Plot, month = Month,
                time = Time, treatment = Treatment) %>%
  unite(Site, Site, Treatment, Plot, Month, Time) %>% rename(PlotID = Site) %>% 
  dplyr::select(PlotID, site, plot, month,
                time, treatment, abundances,life_form) %>% 
  spread(life_form, abundances) %>% 
  rename(bamboo.abn = bamboo,
         herb.abn = herb,
         liana.abn = liana,
         palm.abn = palm,
         shrub.abn = shrub,
         tree.abn = tree) %>% 
  dplyr::select(1,7:12)

bt.abu.mtz

# join inverse simpsons, abundances and richness table --------------------

bt.invsim.abu.mtz <-left_join(bt.invsim.mtz, bt.abu.mtz, by = "PlotID") %>% 
  dplyr::select(12, 1:11, 13:18) %>% 
  hablar::rationalize() %>%
  mutate_all(replace_na, 0) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% 
  hablar::rationalize() %>%
  mutate_all(replace_na, 0)


bt.invsim.abu.mtz


# glmm life-form species evenness -----------------------------------------



# trees -------------------------------------------------------------------

sp.lf.even.time.tree <- glmer(tree.invs~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                          family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.time.tree)
car::Anova(sp.lf.even.time.tree, type = "III")

sp.lf.even.tree <- glmer(tree.invs~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.tree)
car::Anova(sp.lf.even.tree, type = "III")

anova(sp.lf.even.time.tree, sp.lf.even.tree)



# palms -------------------------------------------------------------------

sp.lf.even.time.palm <- glmer(palm.invs~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.time.palm)
car::Anova(sp.lf.even.time.palm, type = "III")


sp.lf.even.palm <- glmer(palm.invs~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.palm)
car::Anova(sp.lf.even.palm, type = "III")

anova(sp.lf.even.time.palm,sp.lf.even.palm)



# lianas -------------------------------------------------------------------


sp.lf.even.time.liana.nb <- glmer.nb(liana.invs~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                    data = bt.invsim.abu.mtz)
summary(sp.lf.even.time.liana.nb)
car::Anova(sp.lf.even.time.liana.nb, type = "III")

sp.lf.even.liana.nb <- glmer.nb(liana.invs~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                                     data = bt.invsim.abu.mtz)
summary(sp.lf.even.liana.nb)
car::Anova(sp.lf.even.liana.nb, type = "III")

anova(sp.lf.even.time.liana.nb,sp.lf.even.liana.nb)



# shrubs -------------------------------------------------------------------


sp.lf.even.time.shrub <- glmer(shrub.invs~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                               family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.time.shrub)
car::Anova(sp.lf.even.time.shrub, type = "III")

sp.lf.even.shrub <- glmer(shrub.invs~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                               family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.shrub)
car::Anova(sp.lf.even.shrub, type = "III")


anova(sp.lf.even.time.shrub,sp.lf.even.shrub)


# herbs -------------------------------------------------------------------

sp.lf.even.time.herb <- glmer(herb.invs~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                               family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.time.herb)
car::Anova(sp.lf.even.time.herb, type = "III")


sp.lf.even.herb <- glmer(herb.invs ~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.herb)
car::Anova(sp.lf.even.herb, type = "III")


anova(sp.lf.even.time.herb,sp.lf.even.herb)



# bamboos -------------------------------------------------------------------


sp.lf.even.time.bamboo <- glmer(bamboo.invs~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.time.bamboo)
car::Anova(sp.lf.even.time.bamboo, type = "III")

sp.lf.even.bamboo <- glmer(bamboo.invs~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                                family = gaussian, data = bt.invsim.abu.mtz)
summary(sp.lf.even.bamboo)
car::Anova(sp.lf.even.bamboo, type = "III")


anova(sp.lf.even.time.bamboo, sp.lf.even.bamboo)
