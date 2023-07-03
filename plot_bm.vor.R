library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data")

##### APEX DATA ################################################################

## Data pathways
# Continuously grazed pastures (TRM)
# TRM <- "D:/APEX model/APEX1905_New/APEX1905_div plot_precip past_OPC diff - TRM/CONUNN_AGM.sad"

# Rotationally grazed pastures (CARM)
# CARM <- "D:/APEX model/APEX1905_New/APEX1905_div plot_precip past_OPC diff - CARM/CONUNN_AGM.sad"

# All 92 plots
all92 <- "D:/APEX model/APEX1905_New/APEX1905_div plot_precip past_OPC diff - all 92/CONUNN_AGM.sad"

## Import data
# TRM
# apexsad_TRM <- read.delim(TRM, sep = "", dec = ".", skip = 8) %>%
#   select(ID, Y, M, D, CPNM, STL, STD, PRCP)

# CARM
# apexsad_CARM <- read.delim(CARM, sep = "", dec = ".", skip = 8) %>%
#   select(ID, Y, M, D, CPNM, STL, STD, PRCP)

# ALL 92
apexsad <- read.delim(all92, sep = "", dec = ".", skip = 8) %>%
    select(ID, Y, M, D, CPNM, STL, STD, PRCP) %>%
  mutate(date = ymd(paste(Y, M, D, sep="-")),
         Year = year(date)) %>% # add date column
  filter(date >= "2014-01-01" & # filter for start of experiment
           month(ymd(date)) %in% c(5:10)) # filter for growing season

# Reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite_92subareas.csv")

## Data preparation
# apexsad <- rbind(apexsad_TRM, apexsad_CARM) %>%
#   mutate(date = ymd(paste(Y, M, D, sep="-")),
#          Year = year(date)) %>% # add date column
#   filter(date >= "2014-01-01" & # filter for start of experiment
#            month(ymd(date)) %in% c(5:10)) # filter for growing season

apex.bm_vor <- apexsad %>% 
  merge(pastID_ecosite, by.x = "ID", by.y = "PastureID") %>%
  mutate(vor.bm = (STL + STD)*1000,
         vor.bm_wt = vor.bm*Proportion)

# VOR biomass at pasture level
apex.bm_vor_pasture <- apex.bm_vor %>%
  group_by(Year, date, Pasture, Ecosite, graze.trt, CARM_name) %>%
  summarize(vor.bm_past = sum(vor.bm_wt))

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

# RS biomass by ecological site
pastID_ecosite_sub <- pastID_ecosite %>% 
  select(PastureID, Pasture, Ecosite, CARM_name) 

rs.biomass_ecosite <- rs.biomass %>%
  merge(pastID_ecosite_sub, by = "Pasture") %>%
  group_by(Year, date, Ecosite) %>%
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

# summarizing by ecological site
field.biomass_ecosite <- field.biomass %>%
  group_by(Year, Season, Pasture, Ecosite, Plot, Transect) %>%
  summarize(bm_transect = mean(HiLo_vor_kgPerha)) %>% # biomass of each transect
  group_by(Year, Season, Pasture, Ecosite, Plot) %>% 
  summarize(bm_plot = mean(bm_transect)) %>%
  group_by(Year, Season, Ecosite) %>% 
  summarize(bm_ecosite = mean(bm_plot), 
            bm_se = sd(bm_plot)/sqrt(length(bm_plot))) %>%
  mutate(month = ifelse(Season == "Spring", 6, 10)) %>% # add date based on season value
  mutate(date = paste(Year, month, 15, sep="-")) %>%
  mutate(date = ymd(date)) %>% # transforming to date format
  filter(date >= "2014-01-01")

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

##### PLOTTING DATA ############################################################

### Pasture, TOTAL biomass
# this code only works for pastures w/ one ecological site
plot_bm.vor_ecosite <- function(past_name) {
  
  # model output for a single pasture
  apex <- apex.bm_vor_pasture %>% filter(Pasture == past_name)
  
  ## pulling out filter criteria
  past_ecosite <- apex$Ecosite %>% unique() # ecological site name
  carm_name <- apex$CARM_name %>% unique() # CARM name, if applicable
  
  ## comparison data sets
  rs <- rs.biomass_ecosite %>% filter(Ecosite == past_ecosite)
  ndvi <- ndvi_cper %>% filter(Pasture == past_name | Pasture == carm_name)
  field <- field.biomass_ecosite %>% filter(Ecosite == past_ecosite)
  
  ## scaling for secondary y axis
  coeff <- max(apex$vor.bm_past)/max(ndvi$NDVI)
  
  plot <- ggplot() +
    geom_line(data = apex, aes(x = date, y = vor.bm_past, color = "APEX"),
              linewidth = 1) +
    geom_line(data = rs, aes(x = date, y = Biomass_kg_ha, color = "Remote Sensing"),
              linewidth = 1) +
    geom_line(data = ndvi, aes(x = date, y = NDVI*coeff, color = "NDVI"),
              linewidth = 1) +
    geom_point(data = field, aes(x = date, y = bm_ecosite)) +
    geom_errorbar(data = field,
                  aes(x = date, 
                      ymin = pmax(0,(bm_ecosite - bm_se)), 
                      ymax = (bm_ecosite + bm_se)),
                  width = 20, color = "black") +
    facet_wrap(.~Year, scales = "free_x") +
    scale_x_date(date_breaks = "1 month",date_labels = "%m") +
    theme_bw() +
    scale_y_continuous("VOR biomass (kg per ha)",
                       sec.axis = sec_axis(~./coeff, name = "NDVI")) +
    ggtitle(paste("Pasture", past_name, ",", 
                  past_ecosite, "Ecosite, CARM name:", carm_name, sep = " "))
  
  return(plot)
}

plot_bm.vor_ecosite(past_name = "7NW")
