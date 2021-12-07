
# Load packages and set directory -------------------------------------

rm(list = ls())

if(!require("tidyverse"))install.packages("tidyverse", dependencies = TRUE)
if(!require("textclean"))install.packages("textclean", dependencies = TRUE)
if(!require("ggpubr"))install.packages("ggpubr", dependencies = TRUE)


# load data table ----------------------------------------------------------

data_biota <- read_csv("00_tables/life_form_2020.csv")
glimpse(data_biota)

data.biota.abun <- data_biota %>%
  dplyr::select(-c(24:28)) %>%  
  gather(key = "Month", value = "value", 5:23) %>% 
  rename(life_form = `Life Form`) %>% 
  mutate(time = Month %>% stringr::str_replace("T", "")) %>% 
  textclean::drop_row("life_form", c("indeterminate", "fern")) %>%
  dplyr::group_by(Site, Plot, Treatment, Month, time, life_form) %>% 
  summarise(abundances = sum(value)) %>% 
  ungroup(Site, Plot, Treatment, Month, time, life_form) %>% 
  mutate(Treatment = fct_relevel(Treatment, c("open", "closed")))

data.biota.abun




# A-Trees -------------------------------------------------------------------


arvores <- data.biota.abun %>% 
  filter(life_form == "tree")
arvores

tree <- ggplot(arvores, aes(time, abundances, color = Treatment, fill = Treatment, shape = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  #stat_summary(fun.data = mean_cl_boot,
  #             geom = "errorbar",
  #             width = 0.2,
  #             aes(group = Treatment),
  #             color = "black",
  #             fun.args = list(conf.int = .95, B = 2000)) +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12,
               show.legend = FALSE
               ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  labs(x = "Sampled period (months)", y = "Tree abundance") +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = c(0.15, .78),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18)) +
  annotate("text", label = 'atop(bold("A"))', parse= TRUE, size = 7, x = 1, y = 19.1) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96","","108")) +
  expand_limits(y=20)
  
tree


# B-Palms -------------------------------------------------------------------


palmeiras <- data.biota.abun %>% 
  filter(`life_form` == "palm")
palmeiras

palm <- ggplot(palmeiras, aes(time, abundances, color = Treatment, fill = Treatment, shape = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12,
               show.legend = FALSE
  ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  labs(x = "", y = "Palm abundance") +
  annotate("text", label = 'atop(bold("B"))', parse= TRUE, size = 7, x = 18.6, y = 11.7) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96","","108")) +
  expand_limits(y=12.5)
palm




# C-Lianas -------------------------------------------------------------------


lianas <- data.biota.abun %>% 
  filter(`life_form` == "liana")
lianas

liana <- ggplot(lianas, aes(time, abundances, color = Treatment, fill = Treatment, shape = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12,
               show.legend = FALSE
  ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  labs(x = "", y = "Liana abundance") +
  annotate("text", label = 'atop(bold("C"))', parse= TRUE, size = 7, x = 1, y = 7.5) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96","","108")) +
  expand_limits(y=8)
liana


# D-Shrubs -------------------------------------------------------------------

arbusto <- data.biota.abun %>% 
  filter(`life_form` == "shrub")
arbusto

shrub <- ggplot(arbusto, aes(time, abundances, color = Treatment, fill = Treatment, shape = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12,
               show.legend = FALSE
  ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  labs(x = "", y = "Shrub abundance") +
  annotate("text", label = 'atop(bold("D"))', parse= TRUE, size = 7, x = 18.6, y = 5.7) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96","","108")) +
  expand_limits(y=6)
shrub


# E-Herbs -------------------------------------------------------------------

ervas <- data.biota.abun %>% 
  filter(`life_form` == "herb")
ervas

herb <- ggplot(ervas, aes(time, abundances, color = Treatment, fill = Treatment, shape = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12,
               show.legend = FALSE
  ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  labs(x = "", y = "Herb abundance") +
  annotate("text", label = 'atop(bold("E"))', parse= TRUE, size = 7, x = 1, y = 12.5) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96","","108")) +
expand_limits(y=13)
herb



# F-bamboos -------------------------------------------------------------------

bamboos <- data.biota.abun %>% 
  filter(`life_form` == "bamboo")
bamboos

bamboo <- ggplot(bamboos, aes(time, abundances, color = Treatment, fill = Treatment, shape = Treatment)) +
  scale_color_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) +
  scale_shape_discrete(name = "Treatment", labels = c("Open","Closed")) +
  scale_fill_manual(values = c("#556B2F","#663399"), name = "Treatment", labels = c("Open","Closed")) + 
  theme_classic() +
  stat_summary(fun.data = mean_se,
               geom = "ribbon",
               aes(group = Treatment, fill = Treatment),
               color = "0.12",
               alpha = 0.12,
               show.legend = FALSE
  ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3,
               aes(shape = Treatment),
               show.legend = TRUE) +
  stat_summary(fun = mean, 
               geom = "line",
               aes(group = Treatment),
               size = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  labs(x = "", y = "Bamboo abundance") +
  annotate("text", label = 'atop(bold("F"))', parse= TRUE, size = 7, x = 18.6, y = 6.7) +
  scale_x_discrete(labels = c("00","","12","","24","","36","","48","","60","","72","","84","","96","","108"))+
expand_limits(y=7)
bamboo


# ggarrange ---------------------------------------------------------------


ggarrange(tree, palm, liana, shrub, herb, bamboo, ncol = 2, nrow = 3)

