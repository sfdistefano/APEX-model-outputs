library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data")

## Reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite.csv")

## APEX output
apexsad <- read.delim("D:/APEX model/APEX1905_New/APEX1905_New/CONUNN_AGM.sad",
                     sep = "", dec = ".", skip = 8) 

# calculting equivalent of VOR biomass data
apex.bm_vor <- apexsad %>%
  merge(pastID_ecosite, by.x = "ID", by.y = "PastureID") %>%
  mutate(date = paste(Y, M, D, sep="-")) %>% # adding date column
  mutate(date = ymd(date)) %>% # transforming into date format
  mutate(vor.bm = (STL + STD)*1000) # biomass left after grazing

## RS biomass data, derived by Sean Kearney
rs.biomass <- read.csv("cper_biomass_means_2014_2022.csv") %>%
  mutate(date = gsub(" 0:00", "", date)) %>% # removing unnecessary info
  mutate(date = mdy(date),
         Biomass_kg_ha = Biomass_lbs_ac*1.121)

## Field VOR biomass data (biomass left after grazing)
field.biomass_jun <- read.csv("CPER VOR/CARM_VOR_JUN_cln_attr_ALL_HiLo2022-10-17.csv") %>%
  select(Year, Season, Pasture, Ecosite, Plot, Transect, Distance, HiLo_vor_kgPerha)

field.biomass_oct <- read.csv("CPER VOR/CARM_VOR_OCT_cln_attr_ALL_HiLo2022-10-19.csv") %>%
  select(Year, Season, Pasture, Ecosite, Plot, Transect, Distance, HiLo_vor_kgPerha)


field.biomass <- rbind(field.biomass_jun, field.biomass_oct) %>%
  na.omit(HiLo_vor_kgPerha) # removing distance points with NA

# summarizing to plot then pasture
field.biomass_past <- field.biomass %>%
  group_by(Year, Season, Pasture, Ecosite, Plot, Transect) %>%
  summarize(bm_transect = mean(HiLo_vor_kgPerha)) %>% # biomass of each transect
  group_by(Year, Season, Pasture, Ecosite, Plot) %>% 
  summarize(bm_plot = mean(bm_transect)) %>% # biomass of each plot
  group_by(Year, Season, Pasture) %>% 
  summarize(bm_past = mean(bm_plot), # biomass of each pasture
            bm_se = sd(bm_plot)/sqrt(length(bm_plot))) # SE across plots

# adding date of collection
field.biomass_past <- field.biomass_past %>%
  mutate(month = ifelse(Season == "Spring", 6, 10)) %>% # add date based on season value
  mutate(date = paste(Year, month, 15, sep="-")) %>%
  mutate(date = ymd(date))

## NDVI
ndvi_cper <- read.csv("cper_ndvi_means_2014_2022.csv") %>%
  mutate(date = gsub(" 0:00", "", date)) %>%
  mutate(date = mdy(date))

## Sub-setting data for pasture
# function
subset_past <- function(data, data.type, past_name) {
  if(data.type == "APEX") {
    past_bm <- data %>% filter(Pasture == past_name & 
                                 date >= "2014-01-01" &
                                 month(ymd(date)) %in% c(5:10)) %>%
      group_by(Y, date, Pasture) %>%
      summarize(vor.bm = sum(vor.bm)) %>%
      mutate(day = yday(date)) %>%
      rename(Year = Y)
    
    return(past_bm)
  }
  
  if(data.type == "RS" | data.type == "NDVI" | data.type == "field") {
    other_bm <- data %>% 
      filter(Pasture == past_name) %>% 
      filter(month(ymd(date)) %in% c(5:10) &
               date >= "2014-01-01") %>%
      mutate(day = yday(date))
    
    return(other_bm)
  }
}

past15E_apex <- subset_past(data = apex.bm_vor, data.type = "APEX", past_name = "15E")
past15E_rs <- subset_past(data = rs.biomass, data.type = "RS", past_name = "15E")
past15E_ndvi <- subset_past(data = ndvi_cper, data.type = "NDVI", past_name = "15E")
past15E_field <- subset_past(data = field.biomass_past, data.type = "field", past_name = "15E")

past19N_apex <- subset_past(data = apex.bm_vor, data.type = "APEX", past_name = "19N")
past19N_rs <- subset_past(data = rs.biomass, data.type = "RS", past_name = "19N")
past19N_ndvi <- subset_past(data = ndvi_cper, data.type = "NDVI", past_name = "19N")
past19N_field <- subset_past(data = field.biomass_past, data.type = "field", past_name = "19N")

## plotting comparison results
plot_bm.vor <- function(apex, rs, ndvi, field, past_name) {
  plot <- ggplot() +
    geom_line(data = apex, aes(x = date, y = vor.bm, color = "APEX")) +
    geom_line(data = rs, aes(x = date, y = Biomass_kg_ha, color = "Remote Sensing")) +
    geom_line(data = ndvi, aes(x = date, y = NDVI*2000, color = "NDVI")) +
    geom_point(data = field, aes(x = date, y = bm_past)) +
    geom_errorbar(data = field,
                  aes(x = date, 
                      ymin = pmax(0,(bm_past - bm_se)), 
                      ymax = (bm_past + bm_se)),
                  width = 20, color = "black") +
    facet_wrap(.~Year, scales = "free_x") +
    scale_x_date(date_breaks = "1 month",date_labels = "%m") +
    theme_bw() +
    ylab("VOR biomass (kg per ha)") +
    xlab("Day of the Year") +
    labs(color = "Data Source") +
    ggtitle(paste("Pasture", past_name))
  
  return(plot)
}

plot_bm.vor(apex = past15E_apex, rs = past15E_rs, 
            ndvi = past15E_ndvi, field = past15E_field,
            past_name = "15E")

plot_bm.vor(apex = past19N_apex, rs = past19N_rs, 
            ndvi = past19N_ndvi, field = past19N_field,
            past_name = "19N")
