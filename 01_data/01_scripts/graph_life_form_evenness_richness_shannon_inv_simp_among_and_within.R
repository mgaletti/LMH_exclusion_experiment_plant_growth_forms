

rm(list = ls())


# packages -----------------------------------------------------------------


library(tidyverse)
library(textclean)
library(vegan)
library(ggpubr)
library(gridExtra)
library(gtable)
library(grid)
library(codyn)
library(ggthemes)

# load table ------------------------------------------------

s.d <- read_csv2("00_tables/species_diversity.csv")



# inverse simpson within life-forms ------------------------------------------------------

s.d <- s.d %>% 
  mutate(Estimates = as.numeric(Estimates)) %>% 
  mutate(se = as.numeric(se)) %>% 
  mutate(min = Estimates - se) %>% 
  mutate(max = Estimates + se) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "Open", "Closed")) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "Tree", "Palm","Liana","Shrub","Herb","Bamboo"))
s.d


# graph -------------------------------------------------------------------

within <- ggplot(s.d, aes(x=life_form, y=Estimates, color = Treatment, fill = life_form)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9), 
                show.legend = FALSE)  + 
  geom_rect(aes(xmin = 0.40, xmax = 1.5, ymin = -Inf, ymax = Inf),fill = "springgreen4", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),fill = "yellowgreen", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),fill = "#fb81bf", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf),fill = "tan1", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),fill = "sienna3", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 5.5, xmax = 6.6, ymin = -Inf, ymax = Inf),fill = "brown", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_hline(yintercept = 0, colour="grey60", linetype = "longdash" ) +
  geom_pointrange(aes(ymin=min, ymax=max),position=position_dodge(.9),  fatten = 5)+
  scale_fill_manual(values = c("springgreen4","yellowgreen",
                               "#fb81bf","tan1","sienna3","brown"), guide="none") +
  scale_color_manual(values = c("#556B2F","#663399")) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = c(0.88, .9),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.title = element_text(size = 16, face="bold", color = "gray10"),
        legend.text = element_text(size = 16, color = "gray20")) +
  #annotate("text", label = 'atop(bold("B"))', parse= TRUE, size = 5, x = 0.6, y = 0.43) +
  ylab("Temporal change in species evenness") + xlab("Growth forms") 

within



# among life-forms -------------------------------------------------------


# load table ----------------------------------------------------------

data_biota <- read_csv("00_tables/life_form_2020_for_abudance.csv")
glimpse(data_biota)


# table with inverse simpsons calculated  -------------------------------------------------------------------


bt.invsim.mtz <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4, 24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = Month %>% stringr::str_replace("T", "")) %>% 
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
  rename(simp.inv.div = ".") %>% 
  hablar::rationalize() %>%
  #dplyr::select(c(1:6)) %>% 
  group_by(Site, Treatment, Plot, Time, Month) %>% 
  dplyr::summarize(Mean = mean(simp.inv.div)) %>% 
  ungroup(Site, Treatment, Plot, Time, Month) %>% 
  mutate(treatment = Treatment,
         time = Time) %>% 
  unite(treatment,treatment, time) %>% rename(PlotID = treatment) %>%
  dplyr::select(c(7,1:6)) %>% 
  mutate(Treatment = factor(Treatment, level = c("open", "closed")))


bt.invsim.mtz

# graphic -----------------------------------------------------------------



among <- ggplot(bt.invsim.mtz, aes(Time, Mean, color = Treatment, fill = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = FALSE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  labs(x = "Sampled period (months)", y = "Inverse Simpson index") +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        legend.position = c(0.15, .9),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96", "", "108"))

among



# life-form richness and shannon diversity among life-forms -----------------


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


# graph
lf_richness <- ggplot(bt.lf.rich, aes(as.factor(Time), n_lf, color = Treatment, fill = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = FALSE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  labs(x = "", y = "Growth form richness") +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96", "", "108"))

lf_richness



# life-form shannon diversity

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

#graph

lf_shannon <- ggplot(bt.shan.div.mtz, aes(as.factor(Time), shannon.div, color = Treatment, fill = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = FALSE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  labs(x = "", y = "Shannon diversity") +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96", "", "108"))

lf_shannon




# ggarrange ---------------------------------------------------------------
ggarrange(among, within, lf_richness, lf_shannon, ncol = 2, labels = c("A", "B", "C", "D"), nrow = 2, font.label = list(size = 18))

#ggsave("02_figures/life_form_species_diversity.jpeg", w = 30, h = 32, unit = "cm", dpi = 1000)


