library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data")


## Reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite_92subareas.csv")

#### APEX output ###############################################################
# All 92 plots
all92 <- "D:/APEX model/APEX1905_Sean/APEX1905_div plot_precip past_OPC diff - all 92/CONUNN_AGM.sad"

apexsad <- read.delim(all92, sep = "", dec = ".", skip = 8) %>%
  select(ID, Y, M, D, CPNM, STL, STD, PRCP, A_DDM) %>%
  mutate(date = ymd(paste(Y, M, D, sep="-")),
         Year = year(date)) %>% # add date column
  filter(date >= "2014-01-01" & # filter for start of experiment
           month(ymd(date)) %in% c(5:10)) # filter for growing season

## calculating equivalent of VOR biomass data
apex.bm_vor <- apexsad %>%
  filter(date >= "2014-01-01" & # filter for start of CARM
           month(ymd(date)) %in% c(5:10)) %>% # filter for growing season
  merge(pastID_ecosite, by.x = "ID", by.y = "PastureID") %>%
  mutate(vor.bm = (STL + STD)*1000,
         vor.bm_wt = vor.bm*Proportion) %>%
  filter(!(CPNM == "VUOC"))

## VOR biomass at pasture level
apex.bm_vor_pasture <- apex.bm_vor %>%
  group_by(Year, date, Pasture, Ecosite, graze.trt, CARM_name) %>%
  summarize(vor.bm_past = sum(vor.bm_wt)) %>% # total standing biomass of a pasture
  mutate(day = yday(date)) 

# divided by ecosites w/in pasture
apex.bm_vor_pasture.es <- apex.bm_vor %>%
  group_by(Year, date, Pasture, Ecosite, graze.trt, CARM_name, Proportion) %>%
  summarize(vor.bm_past.es = sum(vor.bm_wt)) %>% # total standing biomass of a pasture
  mutate(day = yday(date))

apex.bm_vor_pasture <- apex.bm_vor_pasture.es %>%
  group_by(Year, date, day, Pasture, graze.trt, CARM_name, Proportion) %>%
  summarize(vor.bm_past = sum(vor.bm_past.es))

# VOR biomass within grazing treatment
apex.bm_vor_trt <- apex.bm_vor_pasture %>%
  group_by(Year, date, graze.trt) %>%
  summarize(vor.bm_trt = mean(vor.bm_past))

#### RS biomass data, derived by Sean Kearney ####
rs.biomass <- read.csv("Sean K Biomass/cper_biomass_means_2014_2022.csv") %>%
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
# 92 subareas
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

# summarizing to plot
field.biomass_plot <- field.biomass %>%
  group_by(Year, Season, Pasture, Ecosite, Plot, Transect) %>%
  summarize(bm_transect = mean(HiLo_vor_kgPerha)) %>% # biomass of each transect
  group_by(Year, Season, Pasture, Ecosite, Plot) %>% 
  summarize(bm_plot = mean(bm_transect)) %>%
  mutate(month = ifelse(Season == "Spring", 6, 10)) %>% # add date based on season value
  mutate(date = paste(Year, month, 15, sep="-")) %>%
  mutate(date = ymd(date)) %>% # transforming to date format
  filter(date >= "2014-01-01")  %>%
  filter(!(Season == "Fall Pre-Burn"))

# summarizing by pasture
field.biomass_past <- field.biomass_plot %>% # biomass of each plot
  group_by(Year, Season, Pasture, date) %>% 
  summarize(bm_past = mean(bm_plot), # biomass of each pasture
            bm_se = sd(bm_plot)/sqrt(length(bm_plot)))

# summarizing by grazing treatment (TRM vs CARM)
field.biomass_trt <- field.biomass_past %>%
  merge(pastID_ecosite, by = "Pasture") %>%
  group_by(Year, date, graze.trt) %>%
  summarize(bm_trt = mean(bm_past), 
            bm_se = sd(bm_past)/sqrt(length(bm_past)))

# summarizing by ecological site
field.biomass_ecosite <-  field.biomass_plot %>%
  group_by(Year, Season, date, Ecosite) %>% 
  summarize(bm_ecosite = mean(bm_plot), 
            bm_se = sd(bm_plot)/sqrt(length(bm_plot)))

#### NDVI ######################################################################
ndvi_cper <- read.csv("Sean K Biomass/cper_ndvi_means_2014_2022.csv") %>%
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

###### PLOTTING RESULTS ########################################################
##### Plotting comparison results by PASTURE ###################################
plot_bm.vor_past <- function(past_name) {
  
  apex <- apex.bm_vor_pasture %>% filter(Pasture == past_name)
  apex.names <- apex.bm_vor_pasture.es %>% filter(Pasture == past_name)
  
  ## pulling out filter criteria
  past_ecosite <- apex.names$Ecosite %>% unique() # ecological site name
  carm_name <- apex$CARM_name %>% unique() # CARM name, if applicable
  prop.num <- 1/(past_ecosite %>% length())
  
  ## comparison data sets
  ndvi <- ndvi_cper %>% filter(Pasture == past_name | Pasture == carm_name)
  
  rs <- rs.biomass_ecosite %>% 
    filter(Ecosite %in% past_ecosite) %>%
    mutate(Biomass_wt = Biomass_kg_ha*prop.num) %>%
    group_by(Year, date) %>%
    summarize(Biomass_kg_ha = sum(Biomass_wt))
  
  field <- field.biomass_ecosite %>% 
    filter(Ecosite %in% past_ecosite) %>%
    mutate(bm_wt = bm_ecosite*prop.num) %>%
    group_by(Year, Season, date) %>%
    summarize(bm_ecosite = sum(bm_wt))
  
  # SE across plots of the same SE
  field_se <- field.biomass_plot %>%
    filter(Ecosite %in% past_ecosite) %>%
    group_by(Year, Season, date) %>%
    summarize(bm_se = sd(bm_plot)/sqrt(length(bm_plot)))
  
  # coefficient for scaling secondary axis (NDVI)
  coeff <- max(apex$vor.bm_past)/max(ndvi$NDVI)
  
  # Calculate RMSE
  vor_rmse <- rmse(actual = rs$Biomass_kg_ha, predicted = apex$vor.bm_past) %>%
    round(digits = 2)
  
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
                      ymin = pmax(0,(bm_ecosite - field_se$bm_se)),
                      ymax = (bm_ecosite + field_se$bm_se)),
                  width = 20, color = "black") +
    facet_wrap(.~Year, scales = "free_x") +
    scale_x_date(date_breaks = "1 month",date_labels = "%m") +
    theme_bw() +
    xlab("Month of the Year") +
    labs(color = "Data Source") +
    scale_y_continuous("VOR biomass (kg per ha)", limits = c(0, NA),
                       sec.axis = sec_axis(~./coeff, name = "NDVI")) +
    ggtitle(paste("Pasture", past_name, ",",
                  past_ecosite, "Ecosite, CARM name:", carm_name, sep = " ")) +
    theme(text = element_text(size = 15, family = 'serif')) +
    labs(caption = paste("RMSE =", vor_rmse))
  
  return(plot)
}

plot_bm.vor_past(past_name = "19N")

##### Plotting results by grazing TREATMENT (TRM vs AGM) #######################
plot_bm.vor_trt <- function(trt_name) {

  apex <- apex.bm_vor_trt %>% filter(graze.trt == trt_name)
  rs <- rs.biomass_trt %>% filter(graze.trt == trt_name)
  ndvi <- ndvi_trt %>% filter(graze.trt == trt_name)
  field <- field.biomass_trt %>% filter(graze.trt == trt_name)
  
  # coefficient for scaling secondary axis (NDVI)
  coeff <- max(apex$vor.bm_trt)/max(ndvi$NDVI)
  
  # Calculate RMSE
  vor_rmse <- rmse(actual = rs$Biomass_kg_ha, predicted = apex$vor.bm_trt) %>%
    round(digits = 2)

  plot <- ggplot() +
    geom_line(data = apex, aes(x = date, y = vor.bm_trt, color = "APEX"),
              linewidth = 1) +
    geom_line(data = rs, aes(x = date, y = Biomass_kg_ha, color = "Remote Sensing"),
              linewidth = 1) +
    geom_line(data = ndvi, aes(x = date, y = NDVI*coeff, color = "NDVI"),
              linewidth = 1) +
    geom_point(data = field, aes(x = date, y = bm_trt)) +
    geom_errorbar(data = field,
                  aes(x = date,
                      ymin = pmax(0,(bm_trt - bm_se)),
                      ymax = (bm_trt + bm_se)),
                  width = 20, color = "black") +
    facet_wrap(.~Year, scales = "free_x") +
    scale_x_date(date_breaks = "1 month",date_labels = "%m") +
    theme_bw() +
    scale_y_continuous("VOR biomass (kg per ha)", limits = c(0, NA),
                       sec.axis = sec_axis(~./coeff, name = "NDVI")) +
    xlab("Month of the Year") +
    labs(color = "Data Source") +
    ggtitle(paste("Grazing Treatment", trt_name)) +
    theme(text = element_text(size = 15, family = 'serif')) +
    labs(caption = paste("RMSE =", vor_rmse))

  return(plot)
}
 
plot_bm.vor_trt(trt_name = "TRM")
# plot_bm.vor_trt(trt_name = "CARM") # doesn't work b/c only have NDVI for 15E and 19N
