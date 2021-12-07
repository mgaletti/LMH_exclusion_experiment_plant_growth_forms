
# load packages and table -------------------------------------------------



library(tidyverse)
library(stringr)
library(lme4) 
library(car) 
library(naniar)
library(ggpubr)
library(rcompanion)
library(MASS)
library(lmtest)

### Fitting the data

data_biota <- read_csv("00_tables/life_form_2020.csv")
glimpse(data_biota)


lf.ttl.abun <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4,24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>%
  group_by(Site, Treatment, Plot, Month, Time, life_form) %>% 
  summarise(Abundance = sum(value)) %>%
  spread(life_form, Abundance) %>% 
  replace(is.na(.), 0) %>% 
  ungroup() %>% 
  mutate(Treatment = as.factor(Treatment)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open", "closed")))
lf.ttl.abun



# trees -----------------------------------------------------------------

## trees x treatments


mo1lf.tree.treat <- glmer(tree ~ Treatment + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.tree.treat)

car::Anova(mo1lf.tree.treat, type = "III")


## trees x grupos x treatments

# palms
mo1lf.tree.palm <- glmer(tree ~ Treatment*log(palm+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.tree.palm)
car::Anova(mo1lf.tree.palm, type = "III")


# lianas
mo1lf.tree.liana <- glmer(tree ~ Treatment*log(liana+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.tree.liana)

car::Anova(mo1lf.tree.liana, type = "III")

# shrub
mo1lf.tree.shrub <- glmer(tree ~ Treatment*log(shrub+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.tree.shrub)
car::Anova(mo1lf.tree.shrub, type = "III")


# herbs
mo1lf.tree.herbs <- glmer(tree ~ Treatment*log(herb+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.tree.herbs)
car::Anova(mo1lf.tree.herbs, type = "III")


# bamboo
mo1lf.tree.bamboo <- glmer(tree ~ Treatment*log(bamboo+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.tree.bamboo)
car::Anova(mo1lf.tree.bamboo, type = "III")


# palms ---------------------------------------------------------------



## palms x treatments

mo1lf.palm.treat <- glmer(palm ~ Treatment + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.palm.treat)
car::Anova(mo1lf.palm.treat, type = "III")



# lianas
mo1lf.palm.liana <- glmer(palm ~  Treatment*log(liana+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.palm.liana)

car::Anova(mo1lf.palm.liana, type = "III")

mo1lf.palm.no.int <- glmer(palm ~ (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
anova(mo1lf.palm.liana, mo1lf.palm.no.int)


# shrub
mo1lf.palm.shrub <- glmer(palm ~  Treatment*log(shrub+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.palm.shrub)

car::Anova(mo1lf.palm.shrub, type = "III")

mo1lf.palm.no.int <- glmer(palm ~ (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
anova(mo1lf.palm.shrub, mo1lf.palm.no.int)



# herbs
mo1lf.palm.herb <- glmer(palm ~  Treatment*log(herb+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.palm.herb)

car::Anova(mo1lf.palm.herb, type = "III")

mo1lf.palm.no.int <- glmer(palm ~ (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
anova(mo1lf.palm.herb, mo1lf.palm.no.int)



# bamboo
mo1lf.palm.bamboo <- glmer(palm ~  Treatment*log(bamboo+1) + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.palm.bamboo)

car::Anova(mo1lf.palm.bamboo, type = "III")

mo1lf.palm.no.int <- glmer(palm ~ (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
anova(mo1lf.palm.bamboo, mo1lf.palm.no.int)


# lianas ------------------------------------------------------------------


mo1lf.liana.treat.nb <- glmer.nb(liana ~ Treatment + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.liana.treat.nb)

car::Anova(mo1lf.liana.treat.nb, type = "III")


# shrub

mo1lf.liana.shrub.nb <- glmer.nb(liana ~  Treatment*log(shrub+1) + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.liana.shrub.nb)
car::Anova(mo1lf.liana.shrub.nb, type = "III")

mo1lf.liana.no.int <- glmer.nb(liana ~ (1 |Site/Plot/Month), data = lf.ttl.abun)
anova(mo1lf.liana.shrub.nb, mo1lf.liana.no.int)


# herbs

mo1lf.liana.herb.nb <- glmer.nb(liana ~  Treatment*log(herb+1) + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.liana.herb.nb)
car::Anova(mo1lf.liana.herb.nb, type = "III")

mo1lf.liana.no.int <- glmer.nb(liana ~ (1 |Site/Plot/Month), data = lf.ttl.abun)
anova(mo1lf.liana.herb.nb, mo1lf.liana.no.int)


# bamboo

mo1lf.liana.bamboo.nb <- glmer.nb(liana ~  Treatment*log(bamboo+1) + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.liana.bamboo.nb)
car::Anova(mo1lf.liana.bamboo.nb, type = "III")

mo1lf.liana.no.int <- glmer.nb(liana ~ (1 |Site/Plot/Month), data = lf.ttl.abun)
anova(mo1lf.liana.bamboo.nb, mo1lf.liana.no.int)


# shrubs ----------------------------------------------------------------


## shrubs x treatments

mo1lf.shrub.treat.nb <- glmer.nb(shrub ~ Treatment + (1 |Site/Plot/Month),  data = lf.ttl.abun)
summary(mo1lf.shrub.treat.nb)

car::Anova(mo1lf.shrub.treat.nb, type = "III")


# herbs

mo1lf.shrub.herb.nb <- glmer.nb(shrub ~  Treatment*log(herb+1) + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.shrub.herb.nb)
car::Anova(mo1lf.shrub.herb.nb, type = "III")

mo1lf.shrub.no.int <- glmer.nb(shrub ~ (1 |Site/Plot/Month), data = lf.ttl.abun)
anova(mo1lf.shrub.herb.nb, mo1lf.shrub.no.int)

# bamboo

mo1lf.shrub.bamboo.nb <- glmer.nb(shrub ~  Treatment*log(bamboo+1) + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.shrub.bamboo.nb)
car::Anova(mo1lf.shrub.bamboo.nb, type = "III")


mo1lf.shrub.no.int <- glmer.nb(shrub ~ (1 |Site/Plot/Month), data = lf.ttl.abun)
anova(mo1lf.shrub.bamboo.nb, mo1lf.shrub.no.int)


# herbs -------------------------------------------------------------------

## herbs x treatments

mo1lf.herb.treat.nb <- glmer(herb ~ Treatment + (1 |Site/Plot/Month), family = poisson, data = lf.ttl.abun)
summary(mo1lf.herb.treat.nb)

car::Anova(mo1lf.herb.treat.nb, type = "III")


# bamboo

mo1lf.herb.bamboo.nb <- glmer.nb(herb ~  Treatment*log(bamboo+1) + (1 |Site/Plot/Month), data = lf.ttl.abun)
summary(mo1lf.herb.bamboo.nb)
car::Anova(mo1lf.herb.bamboo.nb, type = "III")

mo1lf.herb.no.int <- glmer.nb(herb ~ (1 |Site/Plot/Month), data = lf.ttl.abun)
anova(mo1lf.herb.bamboo.nb, mo1lf.herb.no.int)

