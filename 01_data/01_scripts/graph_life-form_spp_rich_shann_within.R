# load packages  ----------------------------------------------------------


library(tidyverse)
library(textclean)

#################### Within ------------------------------------------------------------------

# load table --------------------------------------------------------------


data_biota <- read_csv("00_tables/life_form_2020_for_spp_rich.csv")
glimpse(data_biota)

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

#graph

spp.rich.graph <- ggplot(bt.spp.rich, aes(x = Month , y = n_spp, color = life_form, group = life_form)) + 
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = life_form),
               color = "0.12",
               alpha = 0.2) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               stroke = 0.5,
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = life_form),
               size = 1.5) +
  theme_bw() +
  facet_grid(Treatment~., space="free", switch="both") +
  labs(y = "Species Richness", x = "Sampled period (months)", color = "Life forms") +
  theme(strip.background = element_rect(color="grey50", fill="gray90"),
        axis.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.y = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position="top", legend.box = "horizontal", 
        plot.background=element_blank(),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(1,1,1,1), "lines")) + #top, right, botton, left
  scale_color_manual(values = c("springgreen4","yellowgreen",
                                "#fb81bf","tan1","sienna3","brown")) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96", "", "108")) +
  expand_limits(y=c(0, 10))

spp.rich.graph


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


# graph

bt.invsim.abu.mtz.graph <- bt.invsim.abu.mtz %>% 
  dplyr::select(!c(13:18)) %>% 
  gather(key = "life_form", value = "shannon.div", 7:12) %>% 
  dplyr::mutate(Treatment = factor(Treatment)) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "open", "closed")) %>% 
  dplyr::mutate(Treatment = recode_factor(Treatment, "open" = "Open", "closed" = "Closed")) %>% 
  dplyr::mutate(life_form = factor(life_form)) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "tree.shan", "palm.shan","liana.shan",
                                        "shrub.shan","herb.shan","bamboo.shan")) %>%
  dplyr::mutate(life_form = recode_factor(life_form, "tree.shan" = "Tree", "palm.shan" = "Palm","liana.shan" = "Liana",
                                          "shrub.shan" = "Shrub","herb.shan" = "Herb","bamboo.shan" = "Bamboo"))


bt.invsim.abu.mtz.graph



spp.shannon.graph <- ggplot(bt.invsim.abu.mtz.graph, aes(x = Month , y = shannon.div, color = life_form, group = life_form)) + 
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = life_form),
               color = "0.12",
               alpha = 0.2) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               stroke = 0.5,
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = life_form),
               size = 1.5) +
  theme_bw() +
  facet_grid(Treatment~., space="free", switch="both") +
  labs(y = "Shannon diversity", x = "Sampled period (months)", color = "Life forms") +
  theme(strip.background = element_rect(color="grey50", fill="gray90"),
        axis.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.y = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position="top", legend.box = "horizontal", 
        plot.background=element_blank(),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(1,1,1,1), "lines")) + #top, right, botton, left
  scale_color_manual(values = c("springgreen4","yellowgreen",
                                "#fb81bf","tan1","sienna3","brown")) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96", "", "108"))


spp.shannon.graph

# ggarrange ---------------------------------------------------------------

ggpubr::ggarrange(spp.rich.graph, spp.shannon.graph,  common.legend = TRUE)

#ggsave("02_figures/spp_rich_shannon_lf_geom_line.jpeg", w = 30, h = 20, unit = "cm", dpi = 1000)



