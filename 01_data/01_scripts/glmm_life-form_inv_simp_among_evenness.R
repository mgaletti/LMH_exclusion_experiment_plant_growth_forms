
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
  mutate(site=Site, treatment = Treatment, plot = Plot, month = Month, time = Time, lifeform = life_form) %>% 
  unite(Site, Site, Treatment, Plot, Month, Time) %>% rename(PlotID = Site) %>% 
  dplyr::group_by(PlotID, site, treatment, plot, month, time, lifeform) %>%
  summarise(abundances = sum(value)) %>% 
  ungroup(PlotID, site, treatment, plot, month, time, lifeform) %>%
  spread(lifeform, abundances) %>% 
  replace(is.na(.), 0) %>% 
  remove_rownames %>% 
  dplyr::select(-c(2:6)) %>% 
  column_to_rownames(var="PlotID") %>% 
  hillR::hill_taxa(q = 2, MARGIN = 1) %>% 
  as.data.frame() %>% rownames_to_column(var="PlotID") %>% 
  separate(PlotID, c("Site", "Treatment", "Plot", "Month", "Time"), convert = FALSE) %>% 
  rename(Inv_simp = ".") %>% 
  hablar::rationalize() %>%
  mutate(site = Site,
         treatment = Treatment,
         plot = Plot,
         month = Month,
         time = Time) %>% 
  unite(site, site, treatment, plot, month, time) %>% rename(PlotID = site)

bt.invsim.mtz  


# table to abundance ------------------------------------------------------



bt.abund.mtz <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::group_by(Site, Treatment, Plot, Month, Time) %>%
  summarise(abundances = sum(value)) %>% 
  mutate(site = Site,
         treatment = Treatment,
         plot = Plot,
         month = Month,
         time = Time) %>% 
  unite(site, site, treatment, plot, month, time) %>% rename(PlotID = site)

bt.abund.mtz 



# join table --------------------------------------------------------------



bt.invsim.abund.mtz <- left_join(bt.invsim.mtz, bt.abund.mtz, by = "PlotID") %>% 
  dplyr::select(-c(7:12)) %>% 
  rename(Site = Site.x,
         Treatment = Treatment.x,
         Plot = Plot.x,
         Month = Month.x,
         Time = Time.x,
         Abundance = abundances) %>% 
  mutate(Time = as.numeric(Time)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open","closed")))
bt.invsim.abund.mtz



# glmm --------------------------------------------------------------------


ttl.lf.even.time.int <- glmer(Inv_simp~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                         family = gaussian, data = bt.invsim.abund.mtz)
summary(ttl.lf.even.time.int)
car::Anova(ttl.lf.even.time.int, type = "III")

ttl.lf.even.time <- glmer(Inv_simp~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                          family = gaussian, data = bt.invsim.abund.mtz)
summary(ttl.lf.even.time)
car::Anova(ttl.lf.even.time, type = "III")


anova(ttl.lf.even.time.int, ttl.lf.even.time)


