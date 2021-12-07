
# Load packages and set directory -------------------------------------

rm(list = ls())

if(!require("tidyverse"))install.packages("tidyverse", dependencies = TRUE)
if(!require("textclean"))install.packages("textclean", dependencies = TRUE)
if(!require("ggpubr"))install.packages("ggpubr", dependencies = TRUE)



# load data table ----------------------------------------------------------

#Precipitation
data_prec <- read_csv2("00_tables/base_santos_prec.csv")
data_temp <- read_csv2("00_tables/base_santos_temp.csv")


data_temp_mean <- data_temp %>% 
  group_by(ano) %>% 
  summarise(mean_temp = mean(temp_c)) 
data_temp_mean

data_prec_temp <- left_join(data_prec, data_temp_mean)
data_prec_temp

clime <- ggplot(data_prec_temp, aes(x = ano)) +
  geom_bar(aes(y = mm_total, color = "Precipitation", fill = "Precipitation"), stat = "identity", width = 0.5) +
  geom_line(aes(y = mean_temp*50, color = "Temperature", fill = "Temperature"), size=1.5) + 
  scale_colour_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  labs(x = "Year") +
  guides(color = FALSE, size = FALSE) +
  theme(legend.position = c(0.925, 0.87),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title=element_blank(),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10))) +
  scale_y_continuous("Total precipitation (mm)", 
                     sec.axis = sec_axis(~. /50, name = "Mean temperature (?C)") # Reverse transformation to match data
  )
clime
