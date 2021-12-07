

rm(list = ls())



# packages -----------------------------------------------------------------


library(tidyverse)



# data --------------------------------------------------------------------
r.a <- read_csv2("00_tables/relative_abundance.csv")
r.a <- r.a %>% 
  mutate(Estimates = as.numeric(Estimates)) %>% 
  mutate(se = as.numeric(se)) %>% 
  mutate(min = Estimates - se) %>% 
  mutate(max = Estimates + se) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "Open", "Closed")) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "Tree", "Palm","Liana","Shrub","Herb","Bamboo"))
  
r.a

s.d <- read_csv2("00_tables/species_diversity.csv")
s.d <- s.d %>% 
  mutate(Estimates = as.numeric(Estimates)) %>% 
  mutate(se = as.numeric(se)) %>% 
  mutate(min = Estimates - se) %>% 
  mutate(max = Estimates + se) %>% 
  dplyr::mutate(Treatment = fct_relevel(Treatment, "Open", "Closed")) %>% 
  dplyr::mutate(life_form = fct_relevel(life_form, "Tree", "Palm","Liana","Shrub","Herb","Bamboo"))
s.d


# relative abundance graph ------------------------------------------------

ggplot(r.a, aes(x=life_form, y=Estimates, color = Treatment, fill = life_form)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9), 
                show.legend = FALSE)  + 
  geom_rect(aes(xmin = 0.40, xmax = 1.5, ymin = -Inf, ymax = Inf),fill = "springgreen4", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf),fill = "yellowgreen", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf),fill = "#fb81bf", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf),fill = "tan1", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf),fill = "sienna3", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_rect(aes(xmin = 5.5, xmax = 6.6, ymin = -Inf, ymax = Inf),fill = "brown", alpha = 0.02, colour = "transparent", show.legend = FALSE) +
  geom_hline(yintercept = 0, colour="grey60", linetype = "longdash" ) +
  geom_pointrange(aes(ymin=min, ymax=max),position=position_dodge(.9),  fatten = 5) +
  scale_fill_manual(values = c("springgreen4","yellowgreen",
                               "#fb81bf","tan1","sienna3","brown"), guide="none") +
  scale_color_manual(values = c("#556B2F","#663399")) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.position = c(0.91, .90), #x e y
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.title = element_text(size=14, face="bold", color = "gray10"),
        legend.text = element_text(size=14, color = "gray20")) +
  annotate("text", label = 'atop(bold(""))', parse= TRUE, size = 5, x = 0.6, y = 0.242) +
  ylab("Temporal change in relative abundance") + xlab("Growth forms")



#ggsave("02_figures/relat_abund_error_bar.jpeg", unit = "cm", dpi = 1000, height = 12, width = 17)

  # species diversity graph -------------------------------------------------

ggplot(s.d, aes(x=life_form, y=Estimates, color = Treatment, fill = life_form)) + 
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
                               "#fb81bf","tan1","sienna3","brown"), guide=FALSE) +
  scale_color_manual(values = c("#556B2F","#663399")) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.position = c(0.91, .88),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.title = element_text(size=12, face="bold", color = "gray10"),
        legend.text = element_text(size=12, color = "gray20")) +
  annotate("text", label = 'atop(bold("B"))', parse= TRUE, size = 5, x = 0.6, y = 0.43) +
  ylab("Temporal change in species evenness") + xlab("Growth forms") 


#ggsave("02_figures/evenness_error_bar.jpeg", unit = "cm", dpi = 1000, height = 12, width = 17)