
# load packages  ----------------------------------------------------------


library(tidyverse)
library(textclean)
library(ggpubr)
library(rcompanion)
library(lme4)


#################### Within ------------------------------------------------------------------

# load table --------------------------------------------------------------


data_biota <- read_csv("00_tables/life_form_2020_for_spp_rich.csv")
glimpse(data_biota)


# species richness and abundance table ----------------------------------------------

bt.spp.rich.abun <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", " ")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" aff. ", " ")) %>%
  textclean::drop_row("Species", c(" cf.")) %>% 
  textclean::drop_row("Species", c(" sp.")) %>% 
  textclean::drop_row("Species", c("sp.")) %>%
  textclean::drop_row("Species", c("eae")) %>% 
  dplyr::mutate(life_form = life_form %>% stringr::str_to_title()) %>%
  dplyr::group_by(Species, life_form) %>%
  dplyr::summarise(abundance  = sum(value))

bt.spp.rich.abun




#write_csv(bt.spp.rich.abun, ".../life_form_2020_spp.abundance.csv")

# calculating species richness ----------------------------------------------

# Species richness
bt.spp.rich <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", ".")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  textclean::drop_row("Species", c(" cf.")) %>% 
  textclean::drop_row("value", c(0)) %>% 
  dplyr::group_by(Site, Plot, Treatment, life_form, Month, Time) %>%
  dplyr::count(Species) %>% dplyr::select(!c(n)) %>% 
  dplyr::count(Species) %>% dplyr::select(!c(Species)) %>% 
  dplyr::summarise(n_spp = sum(n)) %>% 
  dplyr::mutate(Treatment = factor(Treatment)) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "open", "closed")) %>% 
  dplyr::mutate(Treatment = recode_factor(Treatment, "open" = "Open", "closed" = "Closed")) %>% 
  dplyr::mutate(life_form = factor(life_form)) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "tree", "palm", "liana", "shrub", "herb", "bamboo")) %>%
  dplyr::mutate(life_form = recode_factor(life_form, "tree" = "Tree", "palm" = "Palm","liana" = "Liana",
                                   "shrub" = "Shrub","herb" = "Herb","bamboo" = "Bamboo")) 
  
  
bt.spp.rich

######## glmms

bt.spp.rich.mtz <- bt.spp.rich %>% 
  spread(life_form, n_spp) %>% 
  replace(is.na(.), 0)
bt.spp.rich.mtz


#### trees

# with interaction
lf.sp.rich.int.time.tree <- glmer(Tree~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                              family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.int.time.tree)
car::Anova(lf.sp.rich.int.time.tree , type = "III")

# no interaction
lf.sp.rich.time.tree <- glmer(Tree~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.time.tree)
car::Anova(lf.sp.rich.time.tree, type = "III")
  
# best model
anova(lf.sp.rich.int.time.tree, lf.sp.rich.time.tree)



#### palms

# with interaction
lf.sp.rich.int.time.palm <- glmer(Palm~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                  family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.int.time.palm)
car::Anova(lf.sp.rich.int.time.palm , type = "III")

# no interaction
lf.sp.rich.time.palm <- glmer(Palm~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.time.palm)
car::Anova(lf.sp.rich.time.palm, type = "III")

# best model
anova(lf.sp.rich.int.time.palm, lf.sp.rich.time.palm)



#### liana


# with interaction
lf.sp.rich.int.time.liana <- glmer(Liana~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                  family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.int.time.liana)
car::Anova(lf.sp.rich.int.time.liana , type = "III")

# no interaction
lf.sp.rich.time.liana <- glmer(Liana~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.time.liana)
car::Anova(lf.sp.rich.time.liana, type = "III")

# best model
anova(lf.sp.rich.int.time.liana, lf.sp.rich.time.liana)


#### shrub


# with interaction
lf.sp.rich.int.time.shrub <- glmer(Shrub~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                   family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.int.time.shrub)
car::Anova(lf.sp.rich.int.time.shrub , type = "III")

# no interaction
lf.sp.rich.time.shrub <- glmer(Shrub~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                               family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.time.shrub)
car::Anova(lf.sp.rich.time.shrub, type = "III")

# best model
anova(lf.sp.rich.int.time.shrub, lf.sp.rich.time.shrub)



#### herb


# with interaction
lf.sp.rich.int.time.herb <- glmer(Herb~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                   family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.int.time.herb)
car::Anova(lf.sp.rich.int.time.herb , type = "III")

# no interaction
lf.sp.rich.time.herb <- glmer(Herb~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                               family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.time.herb)
car::Anova(lf.sp.rich.time.herb, type = "III")

# best model
anova(lf.sp.rich.int.time.herb, lf.sp.rich.time.herb)



#### bamboo


# with interaction
lf.sp.rich.int.time.bamboo <- glmer(Bamboo~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                  family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.int.time.bamboo)
car::Anova(lf.sp.rich.int.time.bamboo , type = "III")

# no interaction
lf.sp.rich.time.bamboo <- glmer(Bamboo~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = poisson, data = bt.spp.rich.mtz)
summary(lf.sp.rich.time.bamboo)
car::Anova(lf.sp.rich.time.bamboo, type = "III")

# best model
anova(lf.sp.rich.int.time.bamboo, lf.sp.rich.time.bamboo)


# calculating Shannon diversity -------------------------------------------

#shannon diversity 
bt.shan.div.mtz <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", "..")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  textclean::drop_row("Species", c(" cf.")) %>% 
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
  vegan::diversity(., "shannon") %>% 
  as.data.frame() %>% rownames_to_column(var="PlotID") %>%
  separate(PlotID, c("Site", "Treatment", "Plot", "Month", "Time", "Life_form"), convert = TRUE) %>% 
  rename(shannon.div = ".") %>% 
  hablar::rationalize() %>% 
  spread(Life_form, shannon.div) %>%
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
  rename(bamboo.shan = bamboo,
         herb.shan = herb,
         liana.shan = liana,
         palm.shan = palm,
         shrub.shan = shrub,
         tree.shan = tree)

bt.shan.div.mtz
  
# abundance

bt.abun.mtz <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", "..")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  textclean::drop_row("Species", c(" cf.")) %>% 
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

bt.abun.mtz


# join the Shannon diversity and abundance table 

bt.invsim.abu.mtz <-left_join(bt.shan.div.mtz, bt.abun.mtz, by = "PlotID") %>% 
  dplyr::select(12, 1:11, 13:18) %>% 
  hablar::rationalize() %>%
  mutate_all(replace_na, 0) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% 
  hablar::rationalize() %>%
  mutate_all(replace_na, 0)


str(bt.invsim.abu.mtz)



######## glmms
#### trees

# with interaction
lf.sp.shan.int.time.tree <- glmer(tree.shan~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                  family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.int.time.tree)
car::Anova(lf.sp.shan.int.time.tree , type = "III")

# no interaction
lf.sp.shan.time.tree <- glmer(tree.shan~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.time.tree)
car::Anova(lf.sp.shan.time.tree, type = "III")

# best model
anova(lf.sp.shan.int.time.tree, lf.sp.shan.time.tree)



#### palms

# with interaction
lf.sp.shan.int.time.palm <- glmer(palm.shan~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                  family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.int.time.palm)
car::Anova(lf.sp.shan.int.time.palm , type = "III")

# no interaction
lf.sp.shan.time.palm <- glmer(palm.shan~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.time.palm)
car::Anova(lf.sp.shan.time.palm, type = "III")

# best model
anova(lf.sp.shan.int.time.palm, lf.sp.shan.time.palm)



#### liana


# with interaction
lf.sp.shan.int.time.liana <- glmer(liana.shan~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                   family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.int.time.liana)
car::Anova(lf.sp.shan.int.time.liana , type = "III")

# no interaction
lf.sp.shan.time.liana <- glmer(liana.shan~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                               family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.time.liana)
car::Anova(lf.sp.shan.time.liana, type = "III")

# best model
anova(lf.sp.shan.int.time.liana, lf.sp.shan.time.liana)


#### shrub


# with interaction
lf.sp.shan.int.time.shrub <- glmer(shrub.shan~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                   family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.int.time.shrub)
car::Anova(lf.sp.shan.int.time.shrub , type = "III")

# no interaction
lf.sp.shan.time.shrub <- glmer(shrub.shan~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                               family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.time.shrub)
car::Anova(lf.sp.shan.time.shrub, type = "III")

# best model
anova(lf.sp.shan.int.time.shrub, lf.sp.shan.time.shrub)



#### herb


# with interaction
lf.sp.shan.int.time.herb <- glmer(herb.shan~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                  family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.int.time.herb)
car::Anova(lf.sp.shan.int.time.herb , type = "III")

# no interaction
lf.sp.shan.time.herb <- glmer(herb.shan~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                              family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.time.herb)
car::Anova(lf.sp.shan.time.herb, type = "III")

# best model
anova(lf.sp.shan.int.time.herb, lf.sp.shan.time.herb)



#### bamboo


# with interaction
lf.sp.shan.int.time.bamboo <- glmer(bamboo.shan~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                    family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.int.time.bamboo)
car::Anova(lf.sp.shan.int.time.bamboo , type = "III")

# no interaction
lf.sp.shan.time.bamboo <- glmer(bamboo.shan~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                                family = gaussian, data = bt.shan.div.mtz)
summary(lf.sp.shan.time.bamboo)
car::Anova(lf.sp.shan.time.bamboo, type = "III")

# best model
anova(lf.sp.shan.int.time.bamboo, lf.sp.shan.time.bamboo)



#################### among ------------------------------------------------------------------


# load table --------------------------------------------------------------


data_biota <- read_csv("C:/Users/Yuri/Google Drive/Yuri/Mestrado/parcelas_biota/projeto/PPG/tese/publicacoes/JOE/02_revision/01_data/00_tables/life_form_2020_for_abudance.csv")
glimpse(data_biota)


# calculating life-form richness ----------------------------------------------

# life-form richness
bt.lf.rich <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", ".")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  textclean::drop_row("Species", c(" cf.")) %>% 
  textclean::drop_row("value", c(0)) %>% 
  dplyr::group_by(Site, Plot, Treatment, Month, Time) %>%
  dplyr::count(life_form) %>% dplyr::select(!c(n)) %>% 
  dplyr::count(life_form) %>% dplyr::select(!c(life_form)) %>% 
  dplyr::summarise(n_lf = sum(n)) %>% 
  dplyr::mutate(Treatment = factor(Treatment)) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "open", "closed")) %>% 
  dplyr::mutate(Treatment = recode_factor(Treatment, "open" = "Open", "closed" = "Closed")) 
bt.lf.rich





####### glmm

# with interaction
lf.rich.int.time <- glmer(n_lf~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                                    family = poisson, data = bt.lf.rich)
summary(lf.rich.int.time)
car::Anova(lf.rich.int.time , type = "III")

# no interaction
lf.rich.time <- glmer(n_lf~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                        family = poisson, data = bt.lf.rich)
summary(lf.rich.time)
car::Anova(lf.rich.time , type = "III")

# best model
anova(lf.rich.int.time, lf.rich.time)




# calculating Shannon diversity -------------------------------------------

#shannon diversity 
bt.shan.div.mtz <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", "..")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace("..", ".")) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  textclean::drop_row("Species", c(" cf.")) %>% 
  mutate(site=Site, treatment = Treatment, plot = Plot, month = Month, time = Time, lifeform = life_form) %>% 
  unite(Site, Site, Treatment, Plot, Month, Time) %>% rename(PlotID = Site) %>% 
  dplyr::group_by(PlotID, site, treatment, plot, month, time, lifeform) %>%
  dplyr::summarise(abundances = sum(value)) %>% 
  ungroup(PlotID, site, treatment, plot, month, time) %>%
  spread(lifeform, abundances) %>% 
  replace(is.na(.), 0) %>% 
  remove_rownames %>% 
  dplyr::select(-c(2:6)) %>% 
  column_to_rownames(var="PlotID") %>% 
  vegan::diversity(., "shannon") %>% 
  as.data.frame() %>% rownames_to_column(var="PlotID") %>%
  separate(PlotID, c("Site", "Treatment", "Plot", "Month", "Time"), convert = TRUE) %>% 
  rename(shannon.div = ".") %>% 
  hablar::rationalize() %>% 
  mutate(Month = as.factor(Month)) %>%
  mutate(Treatment = as.factor(Treatment)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open", "closed"))) %>% 
  dplyr::mutate(Treatment = recode_factor(Treatment, "open" = "Open", "closed" = "Closed"))
bt.shan.div.mtz



####### glmm

# with interaction
lf.shan.int.time <- glmer(shannon.div~ Treatment*log(Time+1) + (1 |Site/Plot/Month), 
                        family = gaussian, data = bt.shan.div.mtz)
summary(lf.shan.int.time )
car::Anova(lf.shan.int.time  , type = "III")

# no interaction
lf.shan.time <- glmer(shannon.div~ Treatment+log(Time+1) + (1 |Site/Plot/Month), 
                          family = gaussian, data = bt.shan.div.mtz)
summary(lf.shan.time )
car::Anova(lf.shan.time  , type = "III")

# best model
anova(lf.shan.int.time, lf.shan.time)

