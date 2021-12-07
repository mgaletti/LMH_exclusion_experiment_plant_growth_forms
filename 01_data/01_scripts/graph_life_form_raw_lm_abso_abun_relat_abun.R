

rm(list = ls())

# Load packages and open datatable -------------------------------------
if(!require("tidyverse")) install.packages("tidyverse", dependencies = TRUE)
if(!require("textclean")) install.packages("textclean", dependencies = TRUE)
if(!require("ggpubr")) install.packages("ggpubr", dependencies = TRUE)


data_biota <- read_csv("00_tables/life_form_2020_for_abudance.csv")


# Calculate absolute abundances -----------------------------------------


bt.lf.abn <- data_biota %>% 
  rename(life_form = `Life Form`) %>% 
  filter(!life_form == "indeterminate", !life_form == "arborescent fern") %>% 
  dplyr::select(-c(4,24:26,28)) %>% 
  gather(key = "Month", value = "value", 4:22) %>% 
  dplyr::mutate(Time = as.numeric(Month %>% stringr::str_replace("T", ""))) %>% 
  dplyr::mutate(Month = as.numeric(Month %>% stringr::str_replace("T", ""))) %>%
  dplyr::mutate(Month = as.numeric(Month)) %>%
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" cf. ", "..")) %>% 
  dplyr::mutate(Species = Species %>% 
                  stringr::str_replace(" ", ".")) %>% 
  dplyr::group_by(Site, Treatment, Plot, Month, Time, life_form) %>%
  summarise(Abundances = sum(value)) %>% 
  ungroup(Site, Treatment, Plot, Month, Time,life_form) %>%
  spread(life_form, Abundances) %>% 
  replace(is.na(.), 0) %>% 
  rename(bamboo = bamboo,
         herb = herb,
         liana = liana,
         palm = palm,
         shrub = shrub,
         tree = tree) %>% 
  mutate(Treatment = as.factor(Treatment)) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open", "closed"))) %>% 
  gather(key = "life_form", value = "abs_abund", 6:11) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "tree", "palm","liana","shrub","herb","bamboo")) %>%
  dplyr::mutate(life_form = recode(life_form, "tree" = "Tree", "palm" = "Palm","liana" = "Liana","shrub" = "Shrub","herb" = "Herb","bamboo" = "Bamboo"))
  
bt.lf.abn



# graph

ggplot(data = bt.lf.abn, aes(x = log(Month), y = log(abs_abund), shape = Treatment, colour = Treatment, fill = Treatment)) +
  geom_point(size = 2, alpha = 0.2) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), alpha = 0.5) + #
  theme_bw() +
  labs(x = "Sample period (log)",
       y = "Absolute abundance (log)") +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  facet_wrap(~ life_form, nrow = 3, ncol = 2) 



# Calculate relative abundances -----------------------------------------


bt.rel.abun <- data_biota %>% 
  dplyr::rename(life_form = `Life Form`) %>%
  textclean::drop_row("life_form", c("indeterminate", "fern")) %>%
  gather(key = "Month", value = "value", 5:23) %>% 
  dplyr::mutate(time = Month %>% stringr::str_replace("T", "")) %>% 
  dplyr::mutate(Month = Month %>% stringr::str_replace("T", "")) %>% 
  dplyr::mutate(Month = as.numeric(Month)) %>% 
  dplyr::mutate(Treatment = factor(Treatment)) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "open", "closed")) %>% 
  dplyr::mutate(Treatment = recode(Treatment, "open" = "Open", "closed" = "Closed")) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "tree", "palm","liana","shrub","herb","bamboo")) %>%
  dplyr::mutate(life_form = recode(life_form, "tree" = "Tree", "palm" = "Palm","liana" = "Liana","shrub" = "Shrub","herb" = "Herb","bamboo" = "Bamboo")) %>% 
  dplyr::group_by(Site, Plot, Treatment, Month, time, life_form) %>% 
  dplyr::summarise(abs.abun = sum(value)) %>% 
  dplyr::ungroup(Site, Plot, Treatment, Month, time, life_form) %>% 
  dplyr::group_by(Site,Plot,Treatment,Month,time) %>%
  dplyr::mutate(rel.abun = abs.abun / sum(abs.abun)) %>% 
  dplyr::ungroup(Site,Plot,Treatment,Month,time) %>% 
  dplyr::mutate(link_fun_bin = log(rel.abun/(1-rel.abun)))
bt.rel.abun 

# graph

abd.rel.link <- ggplot(data = bt.rel.abun, aes(x = log(Month+1), y = link_fun_bin, 
                               shape = Treatment, colour = Treatment, fill = Treatment)) +
  geom_point(size = 2, alpha = 0.2) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Sample period (log)",
       y = "Logit relative abundance") +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  facet_wrap(~ life_form, nrow = 3, ncol = 2) 
abd.rel.link

abd.rel.log.time <- ggplot(data = bt.rel.abun, aes(x = Month, y = rel.abun, 
                               shape = Treatment, colour = Treatment, fill = Treatment)) +
  geom_point(size = 2, alpha = 0.2) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "Sample period (months)",
       y = "Relative abundance (%)") +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  facet_wrap(~ life_form, nrow = 3, ncol = 2) 

abd.rel.log.time

ggarrange(abd.rel.log.time, abd.rel.link,
          legend = "none",
          labels = c("A","B"),
          ncol = 2, nrow = 1)



#ggsave("02_figures/abs_rel_abun_raw_loess_perc_loglogit.jpeg", width = 35, height = 25, units = "cm", dpi = 500)

