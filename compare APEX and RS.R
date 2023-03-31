library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data")

## APEX output
apexsad = read.delim("D:/APEX model/APEX1905_New_20_constant/CONUNN_AGM_past15E.sad",
                     sep = "", dec = ".", skip = 8) 

# VOR biomass data
apex.bm_vor <- apexsad %>% 
  mutate(date = paste(Y, M, D, sep="-")) %>% # adding date column
  mutate(date = ymd(date)) %>% # transforming into date format
  mutate(vor.bm = (STL + STD)*1000) # biomass left after grazing

# subset for pastures
past15E_bm <- apex.bm_vor %>% filter(ID == 1) %>%
  group_by(Y, ID, date) %>% 
  summarize(vor.bm = sum(vor.bm)) %>%
  filter(date >= "2014-01-01")

past19N_bm <- apex.bm_vor %>% filter(ID == 6) %>%
  group_by(Y, ID, date) %>% 
  summarize(vor.bm = sum(vor.bm))

## RS biomass data
rs.biomass <- read.csv("cper_biomass_means_2014_2022.csv") %>%
  mutate(date = gsub(" 0:00", "", date)) %>% # removing unnecessary info
  mutate(date = mdy(date))

past15E_bm.rs <- rs.biomass %>% filter(Pasture == "15E") %>%
  mutate(Biomass_kg_ha = Biomass_lbs_ac*1.121) # converting to kg per ha

## Field biomass data
field.biomass <- read.csv("cper_vor_june_oct_2013_2021.csv")

# summarizing to plot
field.biomass_plot <- field.biomass %>%
  group_by(Year, Season, Pasture, Plot) %>%
  summarize(bm_plot = mean(bm))

# summarizing to pasture
field.biomass_past <- field.biomass_plot %>%
  group_by(Year, Season, Pasture) %>%
  summarize(bm_past = mean(bm_plot),
            bm_se = sd(bm_plot)/sqrt(length(bm_plot))) %>%
  mutate(month = ifelse(Season == "Spring", 6, 10)) %>%
  mutate(date = paste(Year, month, 15, sep="-")) %>%
  mutate(date = ymd(date))

past15E_bm.field <- field.biomass_past %>% filter(Pasture == "15E")

## subsetting for grazing season (May - Oct)
graze.season_rs <- past15E_bm.rs %>% 
  filter(month(ymd(date)) %in% c(5:10)) %>%
  mutate(day = yday(date))

graze.season_apex <- past15E_bm %>% 
  filter(month(ymd(date)) %in% c(5:10)) %>%
  mutate(day = yday(date)) %>%
  rename(Year = Y)

graze.season_ndvi <- past15E_ndvi %>% 
  filter(month(ymd(date)) %in% c(5:10)) %>%
  mutate(day = yday(date))

graze.season_field <- past15E_bm.field %>% 
  filter(month(ymd(date)) %in% c(5:10)) %>%
  mutate(day = yday(date))

## plotting comparison results
ggplot() +
  geom_line(data = graze.season_apex, aes(x = day, y = vor.bm, color = "APEX")) +
  geom_line(data = graze.season_rs, aes(x = day, y = Biomass_kg_ha, color = "Remote Sensing")) +
  geom_line(data = graze.season_ndvi, aes(x = day, y = NDVI*2000, color = "NDVI")) +
  facet_wrap(~Year) +
  geom_point(data = graze.season_field, aes(x = day, y = bm_past)) +
  geom_errorbar(data = graze.season_field,
                aes(x = day, 
                    ymin = pmax(0,(bm_past - bm_se)), 
                    ymax = (bm_past + bm_se)),
                width = 20, color = "black") +
  theme_bw() +
  ylab("VOR biomass (kg per ha)") +
  xlab("Day of the Year") +
  labs(color = "Data Source") +
  ggtitle("Pasture 15E")
