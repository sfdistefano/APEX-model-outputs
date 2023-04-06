library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data/")

## Reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite.csv")

## APEX output
apexsad <- read.delim("D:/APEX model/APEX1905_New/APEX1905_New/CONUNN_AGM.sad",
                     sep = "", dec = ".", skip = 8) %>%
  select(Y,M,D,ID,X0.5CM:X75.85CM) %>%
  unique()

# summarizing for soil water data (volumetric water content)
apex_sw <- apexsad %>%
  mutate(date = paste(Y, M, D, sep="-")) %>%
  mutate(date = ymd(date)) %>%
  pivot_longer(cols = X0.5CM:X75.85CM, names_to = "depth.range",
               values_to = "vwc.daily") %>%
  mutate(depth.range = gsub("X","", depth.range),
         vwc.daily = vwc.daily*100) %>%
  filter(date >= "2018-01-01")

## Creating reference dataframe for soil depth
depth.range <- unique(apex_sw$depth.range)
depth.cm <- c(2.5, 10, 20, 30,40,50,80)

depth.ref <- data.frame(depth.range, depth.cm)

# adding in depths and pasture names
apex_sw <- merge(apex_sw, depth.ref, by = "depth.range") %>%
  merge(pastID_ecosite, by.x = "ID", by.y = "PastureID")

## CPER soil moisture field data (2018-2022)
cper_sw <- read.csv("19N-15E_soilmoisture_2018-2022.csv") %>%
  mutate(date = mdy(date),
         Y = year(date))

## Plotting volumetric water content
plot_sw <- function(past_name, depth){
  
  apex <- apex_sw %>% filter(Pasture == past_name & depth.cm == depth)
  cper <- cper_sw %>% filter(site == past_name & depth.cm == depth)
  
  plot <- ggplot() +
    geom_line(data = apex, aes(x = date, y = vwc.daily, color = "APEX")) +
    geom_line(data = cper, aes(x = date, y = vwc.daily, color = "CPER"))  +
    scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
    ylab("Volumetric Water Content (%)") +
    xlab("Time") +
    labs(color = "Data Source") +
    ggtitle(paste("Pasture",past_name, ", soil depth =", depth)) +
    theme_bw()
  
  return(plot)
}

plot_sw(past_name = "15E", depth = 50)
