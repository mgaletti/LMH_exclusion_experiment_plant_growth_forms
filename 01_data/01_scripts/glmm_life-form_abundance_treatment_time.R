
# load packages and table ----------------------------------------------------------


library(tidyverse)
library(rcompanion)
library(stringr)
library(MASS)
library(lmtest)
library(lme4) 
library(car)


data_biota <- read_csv("00_tables/life_form_2020.csv")
glimpse(data_biota)


# life form abundance table
bt.lf.abn <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4,24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", "..")) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  dplyr::group_by(Site, Treatment, Plot, Month, Time, life_form) %>%
  summarise(Abundances = sum(value)) %>% 
  ungroup(Site, Treatment, Plot, Month, Time,life_form) %>%
  spread(life_form, Abundances) %>% 
  replace(is.na(.), 0) %>% 
  rename(bamboo.abn = bamboo,
         herb.abn = herb,
         liana.abn = liana,
         palm.abn = palm,
         shrub.abn = shrub,
         tree.abn = tree) %>% 
  mutate(Treatment = as.factor(Treatment)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open", "closed")))
bt.lf.abn


# GLMM abundance ----------------------------------------------------------



# trees -------------------------------------------------------------------


bt.div.sim.tree.glmm.poisson02<- glmer(tree.abn~ Treatment*log(Time+1) + (1 |Site/Plot/Month), family=poisson, data = bt.lf.abn)
summary(bt.div.sim.tree.glmm.poisson02)
car::Anova(bt.div.sim.tree.glmm.poisson02, type = "III")

mo1lf.tree <- glmer(tree.abn~ Treatment+log(Time+1) + (1 |Site/Plot/Month), family=poisson, data = bt.lf.abn)
summary(mo1lf.tree)

anova(bt.div.sim.tree.glmm.poisson02,mo1lf.tree)

# palms -------------------------------------------------------------------


bt.div.sim.palm.glmm.poisson02<- glmer(palm.abn~ Treatment*log(Time+1) + (1 |Site/Plot/Month), family=poisson, data = bt.lf.abn)
summary(bt.div.sim.palm.glmm.poisson02)
car::Anova(bt.div.sim.palm.glmm.poisson02, type = "III")

mo1lf.palm<- glmer(palm.abn~ Treatment+log(Time+1) + (1 |Site/Plot/Month), family=poisson, data = bt.lf.abn)
summary(mo1lf.palm)

anova(bt.div.sim.palm.glmm.poisson02,mo1lf.palm)


# lianas ------------------------------------------------------------------


nbglmer.liana <- glmer.nb(liana.abn~ Treatment*log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(nbglmer.liana)
car::Anova(nbglmer.liana, type = "III")

mo1lf.liana <- glmer.nb(liana.abn~ Treatment+log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(mo1lf.liana)

anova(nbglmer.liana,mo1lf.liana)

lianas.t108 <- bt.lf.abn %>% 
  filter(!Month == "T108")
lianas.t108

nbglmer.liana.t108 <- glmer.nb(liana.abn~ Treatment*log(Time+1)+ (1 |Site/Plot/Month), data = lianas.t108)
summary(nbglmer.liana.t108)
car::Anova(nbglmer.liana.t108, type = "III")

nbglmer.liana.t108.ni <- glmer.nb(liana.abn~ Treatment+log(Time+1)+ (1 |Site/Plot/Month), data = lianas.t108)
summary(nbglmer.liana.t108.ni)

anova(nbglmer.liana.t108, nbglmer.liana.t108.ni)




# shrubs ------------------------------------------------------------------


nbglmer.shrub<- glmer.nb(shrub.abn~ Treatment*log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(nbglmer.shrub)
car::Anova(nbglmer.shrub, type = "III")

mo1lf.shrub<- glmer.nb(shrub.abn~ Treatment+log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(mo1lf.shrub)

anova(nbglmer.shrub,mo1lf.shrub)


# herbs ------------------------------------------------------------------


nbglmer.herb<- glmer.nb(herb.abn~ Treatment*log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(nbglmer.herb)
car::Anova(nbglmer.herb, type = "III")

mo1lf.herb<- glmer.nb(herb.abn~ Treatment+log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(mo1lf.herb)

anova(nbglmer.herb,mo1lf.herb)


# bamboos -----------------------------------------------------------------


nbglmer.bamboo<- glmer.nb(bamboo.abn~ Treatment*log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(nbglmer.bamboo)
car::Anova(nbglmer.bamboo, type = "III")

mo1lf.bamboo<- glmer.nb(bamboo.abn~ Treatment+log(Time+1)+ (1 |Site/Plot/Month), data = bt.lf.abn)
summary(mo1lf.bamboo)

anova(nbglmer.bamboo,mo1lf.bamboo)

