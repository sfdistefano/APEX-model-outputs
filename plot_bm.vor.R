library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data")

## Reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite_92subareas.csv")

#### APEX output ###############################################################
apexsad <- read.delim("D:/APEX model/APEX1905_New/APEX1905_div plot_precip past_OPC diff - COPY/CONUNN_AGM.sad",
                      sep = "", dec = ".", skip = 8) %>%
  select(ID, Y, M, D, CPNM, STL, STD, PRCP)

## calculating equivalent of VOR biomass data
apex.bm_vor <- apexsad %>%
  mutate(date = paste(Y, M, D, sep="-")) %>% # adding date column
  mutate(date = ymd(date)) %>% # transforming into date format
  rename(Year = Y) %>%
  filter(date >= "2014-01-01" & # filter for start of CARM
           month(ymd(date)) %in% c(5:10)) %>% # filter for growing season
  merge(pastID_ecosite, by.x = "ID", by.y = "PastureID") %>%
  mutate(vor.bm = (STL + STD)*1000,
         vor.bm_wt = vor.bm*Proportion)

# VOR biomass at pasture level
apex.bm_vor_pasture <- apex.bm_vor %>%
  group_by(Year, date, Pasture, graze.trt) %>%
  summarize(vor.bm_past = sum(vor.bm_wt)) %>% # total standing biomass of a pasture
  mutate(day = yday(date))
# VOR biomass within grazing treatment
apex.bm_vor_trt <- apex.bm_vor_pasture %>%
  group_by(Year, date, graze.trt) %>%
  summarize(vor.bm_trt = mean(vor.bm_past))

#### RS biomass data, derived by Sean Kearney ####
rs.biomass <- read.csv("cper_biomass_means_2014_2022.csv") %>%
  mutate(date = gsub(" 0:00", "", date)) %>% # removing unnecessary info
  mutate(date = mdy(date),
         Biomass_kg_ha = Biomass_lbs_ac*1.121) %>% 
  filter(date >= "2014-01-01" & # filter for start of CARM
           month(ymd(date)) %in% c(5:10))

# RS biomass within grazing treatment
rs.biomass_trt <- rs.biomass %>%
  merge(pastID_ecosite, by = "Pasture") %>%
  group_by(Year, date, graze.trt) %>%
  summarize(Biomass_kg_ha = mean(Biomass_kg_ha))

#### Field biomass data ########################################################
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
            bm_se = sd(bm_plot)/sqrt(length(bm_plot))) %>%
  mutate(month = ifelse(Season == "Spring", 6, 10)) %>% # add date based on season value
  mutate(date = paste(Year, month, 15, sep="-")) %>%
  mutate(date = ymd(date)) %>% # transforming to date format
  filter(date >= "2014-01-01")

# summarizing by grazing treatment (TRM vs CARM)
field.biomass_trt <- field.biomass_past %>%
  merge(pastID_ecosite, by = "Pasture") %>%
  group_by(Year, date, graze.trt) %>%
  summarize(bm_trt = mean(bm_past), 
            bm_se = sd(bm_past)/sqrt(length(bm_past)))

#### NDVI ######################################################################
ndvi_cper <- read.csv("cper_ndvi_means_2014_2022.csv") %>%
  mutate(date = gsub(" 0:00", "", date)) %>%
  mutate(date = mdy(date)) %>% # transforming to date format
  filter(month(ymd(date)) %in% c(5:10) &
           date >= "2014-01-01") %>% # filter for start of CARM
  mutate(day = yday(date)) # adding day column

# NDVI within grazing treatment
ndvi_trt <- ndvi_cper %>%
  merge(pastID_ecosite, by = "Pasture") %>%
  group_by(Year, date, graze.trt) %>%
  summarize(NDVI = mean(NDVI))

##### PLOTTING RESULTS #########################################################
#### Plotting comparison results by PASTURE
plot_bm.vor_past <- function(past_name) {
  
  apex <- apex.bm_vor_pasture %>% filter(Pasture == past_name)
  rs <- rs.biomass %>% filter(Pasture == past_name)
  ndvi <- ndvi_cper %>% filter(Pasture == past_name)
  field <- field.biomass_past %>% filter(Pasture == past_name)
  
  coeff <- max(apex$vor.bm_past)/max(ndvi$NDVI)
  
  plot <- ggplot() +
    geom_line(data = apex, aes(x = date, y = vor.bm_past, color = "APEX")) +
    geom_line(data = rs, aes(x = date, y = Biomass_kg_ha, color = "Remote Sensing")) +
    geom_line(data = ndvi, aes(x = date, y = NDVI*coeff, color = "NDVI")) +
    geom_point(data = field, aes(x = date, y = bm_past)) +
    geom_errorbar(data = field,
                  aes(x = date, 
                      ymin = pmax(0,(bm_past - bm_se)), 
                      ymax = (bm_past + bm_se)),
                  width = 20, color = "black") +
    facet_wrap(.~Year, scales = "free_x") +
    scale_x_date(date_breaks = "1 month",date_labels = "%m") +
    theme_bw() +
    scale_y_continuous("VOR biomass (kg per ha)",
                       sec.axis = sec_axis(~./coeff, name = "NDVI")) +
    xlab("Month of the Year") +
    labs(color = "Data Source") +
    ggtitle(paste("Pasture", past_name))
  
  return(plot)
}

# plot_bm.vor_past(past_name = "15E")
plot_bm.vor_past(past_name = "15E")

#### Plotting results by grazing TREATMENT (TRM vs AGM)
plot_bm.vor_trt <- function(trt_name) {
  
  apex <- apex.bm_vor_trt %>% filter(graze.trt == trt_name) 
  rs <- rs.biomass_trt %>% filter(graze.trt == trt_name)
  ndvi <- ndvi_trt %>% filter(graze.trt == trt_name) 
  field <- field.biomass_trt %>% filter(graze.trt == trt_name)
  
  coeff <- max(apex$vor.bm_trt)/max(ndvi$NDVI)
  
  plot <- ggplot() +
    geom_line(data = apex, aes(x = date, y = vor.bm_trt, color = "APEX")) +
    geom_line(data = rs, aes(x = date, y = Biomass_kg_ha, color = "Remote Sensing")) +
    geom_line(data = ndvi, aes(x = date, y = NDVI*coeff, color = "NDVI")) +
    geom_point(data = field, aes(x = date, y = bm_trt)) +
    geom_errorbar(data = field,
                  aes(x = date, 
                      ymin = pmax(0,(bm_trt - bm_se)), 
                      ymax = (bm_trt + bm_se)),
                  width = 20, color = "black") +
    facet_wrap(.~Year, scales = "free_x") +
    scale_x_date(date_breaks = "1 month",date_labels = "%m") +
    theme_bw() +
    scale_y_continuous("VOR biomass (kg per ha)",
                       sec.axis = sec_axis(~./coeff, name = "NDVI")) +
    xlab("Month of the Year") +
    labs(color = "Data Source") +
    ggtitle(paste("Treatment", trt_name))
  
  return(plot)
}

plot_bm.vor_trt(trt_name = "TRM")
plot_bm.vor_trt(trt_name = "CARM") # doesn't work b/c only have NDVI for 15E and 19N
