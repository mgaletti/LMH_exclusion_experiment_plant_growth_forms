
# load packages -----------------------------------------------------------


library(tidyverse)
library(circlize)
library(chorddiag)
library(patchwork)
library(hrbrthemes)
library(migest)

biota <- read_csv("00_tables/biota_interaction_life_forms_summed_estimates.csv")


# Matrix --------------------------------------------------------------
# open --------------------------------------------------------------------


tiff("02_figures/BIOTA_cord_graphic_open.tiff", width = 8, height = 8, units = 'cm', res = 1000)

par(mfrow = c(1, 1))

biota.open <- biota %>% 
  filter(treatment == "open") %>% 
  rename(from = "life_form1",
         to = "life_form2",
         value = "estimates") %>% 
  dplyr::mutate(from = from %>% stringr::str_to_title()) %>% 
  dplyr::mutate(to = to %>% stringr::str_to_title()) %>% 
  select(2:4) %>% 
  spread(to, value = "value") %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("from")

biota.open

biota.openm = biota.open
for(cn in intersect(rownames(biota.openm), colnames(biota.openm))) {
  biota.openm[cn, cn] = 0
}
biota.openm

biota.openm <- as.matrix(biota.openm)

circos.clear()

circos.par(gap.after = c("Tree" = 7, "Palm" = 7, "Liana" = 7, "Shrub" = 7, "Herb" = 7,
                         "Bamboo" = 7),
           start.degree = 30)

biota_open_lty_mat = matrix(1, nrow = nrow(biota.openm), ncol = ncol(biota.openm))
biota_open_lty_mat[biota.openm < 0 & biota.openm > -0.050] = 3
biota_open_lty_mat[biota.openm > 0 & biota.openm < 0.050] = 3
biota_open_lty_mat


biota_open_border_mat = matrix(NA, nrow = nrow(biota.openm), ncol = ncol(biota.openm))
biota_open_border_mat[biota.openm < 0 ] = "red"
biota_open_border_mat[biota.openm >= 0 ] = "black"
biota_open_border_mat


col_mat = c(Tree = "#0ac85f", Palm = "#0ac85f",
            Liana = "#0ac85f", Shrub = "#0ac85f",
            Herb = "#0ac85f", Bamboo = "#0ac85f")

col_mat[biota.openm < 0.00] = "#DDA0DD"
col_mat[biota.openm > 0.00] = "#0ac85f"

chordDiagram(biota.openm, 
             link.border = biota_open_border_mat,
             link.lty = biota_open_lty_mat,
             order = c("Liana", "Shrub", "Herb", "Bamboo", "Tree", "Palm"),
             grid.col = c(Tree = "#0ac85f", Palm = "#0ac85f",
                          Liana = "#0ac85f", Shrub = "#0ac85f",
                          Herb = "#0ac85f", Bamboo = "#0ac85f"), 
             col = col_mat,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.01, 0.03),
             diffHeight = mm_h(2),
             grid.border = TRUE,
             transparency = 0.2)

circos.clear()

dev.off()



# closed ------------------------------------------------------------------


tiff("02_figures//BIOTA_cord_graphic_closed.tiff", width = 8, height = 8, units = 'cm', res = 1000)

par(mfrow = c(1, 1))

biota.closed <- biota %>% 
  filter(treatment == "closed") %>% 
  rename(from = "life_form1",
         to = "life_form2",
         value = "estimates") %>% 
  dplyr::mutate(from = from %>% stringr::str_to_title()) %>% 
  dplyr::mutate(to = to %>% stringr::str_to_title()) %>% 
  select(2:4) %>% 
  spread(to, value = "value") %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("from")

biota.closed

biota.closedm = biota.closed
for(cn in intersect(rownames(biota.closedm), colnames(biota.closedm))) {
  biota.closedm[cn, cn] = 0
}
biota.closedm

biota.closedm <- as.matrix(biota.closedm)

circos.clear()

circos.par(gap.after = c("Tree" = 7, "Palm" = 7, "Liana" = 7, "Shrub" = 7, "Herb" = 7,
                         "Bamboos" = 7),
           start.degree = 60)


biota_lty_mat = matrix(NA, nrow = nrow(biota.closedm), ncol = ncol(biota.closedm))
biota_lty_mat[biota.closedm < 0 & biota.closedm > -0.050] = 3
biota_lty_mat[biota.closedm > 0 & biota.closedm < 0.050] = 3
biota_lty_mat


biota_border_mat = matrix(NA, nrow = nrow(biota.closedm), ncol = ncol(biota.closedm))
biota_border_mat[biota.closedm < 0 ] = "black"
biota_border_mat[biota.closedm >= 0 ] = "black"
biota_border_mat


col_mat = c(Tree = "#0ac85f", Palm = "#0ac85f",
            Liana = "#0ac85f", Shrub = "#0ac85f",
            Herb = "#0ac85f", Bamboo = "#0ac85f")

col_mat[biota.closedm < 0.00] = "#DDA0DD"
col_mat[biota.closedm > 0.00] = "#0ac85f"

chordDiagram(biota.closedm, 
             link.border = biota_border_mat,
             order = c("Liana", "Shrub", "Herb", "Bamboo", "Tree", "Palm"),
             grid.col = c(Tree = "#0ac85f", Palm = "#0ac85f",
                          Liana = "#0ac85f", Shrub = "#0ac85f",
                          Herb = "#0ac85f", Bamboo = "#0ac85f"), 
             col = col_mat,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.01, 0.03),
             diffHeight = mm_h(2),
             grid.border = TRUE,
             transparency = 0.2)


circos.clear()

dev.off()

