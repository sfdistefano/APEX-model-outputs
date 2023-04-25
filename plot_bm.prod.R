library(tidyverse)
library(lubridate)
library(zoo)

setwd("D:/APEX data and scripts/Data/")

## reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite.csv")

## APEX output
apexsad <- read.delim("D:/APEX model/APEX1905_New/APEX1905_New/CONUNN_AGM.sad",
                     sep = "", dec = ".", skip = 8) %>%
  mutate(date = paste(Y,M,D, sep = "-")) %>%
  mutate(date = ymd(date))
# Filtering for growing season (May - Oct)
apex_growing.season <- apexsad %>% 
  filter(M >= 5 & M <= 10) %>% # June - October
  filter(Y >= 2013) # starting year of CARM

## Calculating cumulative sum --> biomass production
# cumulative sum for EACH CPNM
apex_cumsum_CPNM <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, A_DDM) %>%
  group_by(ID, Y, CPNM) %>%
  dplyr::mutate(cum_DDM = cumsum(A_DDM)*1000) # conversion from metric ton/ha to kg/ha

# cumulative sum of ALL CPNM (i.e., total biomass production)
apex_cumsum_TOTAL <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, A_DDM) %>%
  group_by(date, ID, Y) %>%
  summarize(A_DDM = sum(A_DDM)) %>%
  group_by(ID, Y) %>%
  dplyr::mutate(cum_DDM = cumsum(A_DDM)*1000)

## STD and STL
apex_stdl <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, STLD.INTAKEkg.ha) %>%
  group_by(date, Y, ID) %>%
  summarize(STDL = sum(STLD.INTAKEkg.ha))

apex_stdl_CPNM <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, STLD.INTAKEkg.ha) %>%
  group_by(date, Y, ID, CPNM) %>%
  summarize(STDL = sum(STLD.INTAKEkg.ha))

## original cage biomass used (day*pasture*plot*functional group)
biomass.plot <- read.csv("CPER Biomass/CARM_Biomass_cln_attr2023-01-19_DJA_pivots.csv") %>% 
  group_by(YearSampled, Treatment, Pasture, Ecosite, Plot, FGCode) %>%
  summarize(kgPerHa = mean(kgPerHa))%>% # summarize by plot of each pasture
  filter(!(Plot == 5 & Pasture == "18S")) %>% # prescribed burn plots
  filter(!(Plot == 6 & Pasture == "18S")) %>%
  filter(!(Plot == 5 & Pasture == "19N")) %>%
  filter(!(Plot == 6 & Pasture == "19N")) %>%
  rename(Year = YearSampled)

# preparing dataframe for plotting SE across plots
biomass.plot_info <- merge(biomass.plot, pastID_ecosite,
                           by = c("Pasture", "Ecosite"), all.x = TRUE) %>%
  mutate(Month = 8, Day = 12,
         date = ymd(paste(Year, Month, Day, sep = "-"))
  ) %>% # adding date info, when collected in field
  rename(CPNM = FGCode, Y = Year)

biomass.plot_info$CPNM2 <- biomass.plot_info$CPNM %>% 
  gsub(pattern = "C3PG", replacement = "CSPG") %>%
  gsub(pattern = "FORB", replacement = "FRB3") %>%
  gsub(pattern = "SS", replacement = "SSHB")

# biomass data is summarized to day*pasture*functional group
biomass.pasture <-  biomass.plot %>% 
  mutate(Month = 8, Day = 12,
         date = ymd(paste(Year, Month, Day, sep = "-"))
         ) %>% # adding date info, when collected in field
  group_by(date, Year, Pasture, Ecosite, Treatment, FGCode) %>%
  summarize(Biomass = mean(kgPerHa)) %>% # summarized by pasture
  rename(CPNM = FGCode, Y = Year) %>%
  filter(!(Pasture %in% c("1W", "28N", "32W"))) %>% # pastures not included in CARM or TRM
  filter(!(CPNM == "SD")) # removing standing dead (growth from previous season)

## adding pasture ID to match with APEX output files
biomass.pasture <- merge(biomass.pasture, pastID_ecosite, 
                         by = c("Pasture", "Ecosite"), all.x = TRUE) 

## fixing functional groups names so they match APEX names
biomass.pasture$CPNM2 <- biomass.pasture$CPNM %>% 
  gsub(pattern = "C3PG", replacement = "CSPG") %>%
  gsub(pattern = "FORB", replacement = "FRB3") %>%
  gsub(pattern = "SS", replacement = "SSHB")

## function for plotting biomass
# NOTE: SE across plots vs pastures aren't much different
plot_bm.prod <- function(past.output, past.name, ecosite, veg.output, func.group) {
  
  if(veg.output == "CPNM" & past.output == "Pasture") {
    
    cumsum <- merge(x = apex_cumsum_CPNM, y = pastID_ecosite, 
                    by.x = "ID", by.y = "PastureID") %>% 
      filter(Pasture == past.name & CPNM == func.group) %>%
      mutate(cum_DDM_wt = cum_DDM*Proportion) %>% # weighting values by proportion of pasture
      group_by(Y, date, Pasture, CPNM) %>%
      summarize(cum_DDM = sum(cum_DDM_wt))
    
    stdl <- merge(x = apex_stdl_CPNM, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>% 
      filter(Pasture == past.name & CPNM == func.group) %>%
      mutate(STDL_wt = STDL*Proportion) %>%
      group_by(Y, date, Pasture, CPNM) %>%
      summarize(STDL = sum(STDL_wt))
    
    # mean calculated at pasture level
    field01 <- biomass.pasture %>% 
      filter(Pasture == past.name & CPNM2 == func.group) %>%
      mutate(Biomass_wt = Biomass*Proportion) %>%
      group_by(Y, date, Pasture, CPNM2) %>%
      summarize(Biomass = sum(Biomass_wt)) # sum of 1 func.group across ecosites w/in pasture
    
    # SE calculated from plots; some pastures may have 1 plot, SE = 0
    field02 <- biomass.plot_info %>%
      filter(Pasture == past.name & CPNM2 == func.group) %>%
      mutate(Biomass_wt = kgPerHa*Proportion) %>%
      group_by(date, Y, Pasture, CPNM2) %>%
      summarize(se = sd(Biomass_wt)/sqrt(length(Biomass_wt)))
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass")) +
      geom_point(data = field01, aes(x = date, y = Biomass)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass - se, 
                        ymax = field01$Biomass + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass Production (kg/ha)") +
      xlab("Month of the Year (Jun - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Pasture", past.name, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "CPNM" & past.output == "Ecosite") {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      group_by(date, Y, Ecosite, CPNM) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM)) %>%
      filter(Ecosite == ecosite, CPNM == func.group)
    
    stdl <- merge(x = apex_stdl_CPNM, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      group_by(date, Y, Ecosite, CPNM) %>%
      summarize(STDL = mean(STDL)) %>%
      filter(Ecosite == ecosite & CPNM == func.group)
    
    field <- biomass.pasture %>%
      group_by(date, Y, Ecosite, CPNM2) %>%
      summarize(Biomass_mean = mean(Biomass), # mean across pastures
                se = sd(Biomass)/sqrt(length(Biomass))) %>% # standard error across pastures
      filter(Ecosite == ecosite, CPNM2 == func.group)
    
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production"))  +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass")) +
      geom_point(data = field, aes(x = date, y = Biomass_mean)) +
      geom_errorbar(data = field,
                    aes(x = date, ymin = Biomass_mean - se, ymax = Biomass_mean + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass Production (kg/ha)") +
      xlab("Month of the Year (Jun - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Ecosite", ecosite, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "Total" & past.output == "Pasture"){
    
    cumsum <- merge(x = apex_cumsum_TOTAL, y = pastID_ecosite, 
                    by.x = "ID", by.y = "PastureID") %>% 
      filter(Pasture == past.name) %>%
      mutate(cum_DDM_wt = cum_DDM*Proportion) %>%
      group_by(date, Y, Pasture) %>%
      summarize(cum_DDM = sum(cum_DDM_wt))
    
    stdl <- merge(x = apex_stdl, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>% 
      filter(Pasture == past.name) %>%
      mutate(stdl_wt = STDL*Proportion) %>%
      group_by(date, Y, Pasture) %>%
      summarize(STDL = sum(stdl_wt))
    
    # mean calculated at pasture level
    field01 <- biomass.pasture %>%
      filter(Pasture == past.name) %>%
      mutate(Biomass_wt = Biomass*Proportion) %>%
      group_by(date, Y, Pasture) %>%
      summarize(Biomass_total = sum(Biomass_wt)) # sum per pasture
    
    # SE calculated from plots; some pastures may have 1 plot, SE = 0            
    field02 <- biomass.plot_info %>%
      filter(Pasture == past.name) %>%
      mutate(Biomass_wt = kgPerHa*Proportion) %>%
      group_by(date, Y, Pasture, Ecosite, Plot) %>%
      summarize(Biomass_wt = sum(Biomass_wt)) %>% # total biomass for each plot
      group_by(date, Y, Pasture) %>%
      summarize(se = sd(Biomass_wt)/sqrt(length(Biomass_wt))) # SE across all plots of a pasture
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass")) +
      geom_point(data = field01, aes(x = date, y = Biomass_total)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_total - se, 
                        ymax = field01$Biomass_total + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (Jun - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Pasture", past.name, "TOTAL biomass"))
    
    return(plot)
  }
  
  if(veg.output == "Total" & past.output == "Ecosite") {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      group_by(date, Y, Ecosite, Pasture) %>%
      summarize(A_DDM = sum(A_DDM),
                cum_DDM = sum(cum_DDM)) %>%
      group_by(date, Y, Ecosite) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM)) %>%
      filter(Ecosite == ecosite)
    
    stdl <- merge(x = apex_stdl, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      group_by(date, Y, Ecosite) %>%
      summarize(STDL = mean(STDL)) %>%
      filter(Ecosite == ecosite)
    
    field <- biomass.pasture %>%
      group_by(date, Y, Ecosite, Pasture) %>%
      summarize(Biomass_total = sum(Biomass), 
                se = sd(Biomass)/sqrt(length(Biomass))) %>%
      group_by(date, Y, Ecosite) %>%
      summarize(Biomass_mean = mean(Biomass_total), # mean across pastures
                se = sd(Biomass_total)/sqrt(length(Biomass_total))) %>% # SE across pastures
      filter(Ecosite == ecosite)
    
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass")) +
      geom_point(data = field, aes(x = date, y = Biomass_mean)) +
      geom_errorbar(data = field,
                    aes(x = date, ymin = Biomass_mean - se,
                        ymax = Biomass_mean + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass Production (kg/ha)") +
      xlab("Month of the Year (Jun - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Ecosite", ecosite, "TOTAL biomass"))
    
    return(plot)
    
  }
}

plot_bm.prod(past.output = "Pasture", past.name = "15E", 
             veg.output = "CPNM", func.group = "BOBU")

plot_bm.prod(past.output = "Pasture", past.name = "7NW", 
             veg.output = "Total")


plot_bm.prod(past.output = "Ecosite", ecosite = "Sandy", 
             veg.output = "CPNM", func.group = "BOBU")

plot_bm.prod(past.output = "Ecosite", ecosite = "Loamy",
             veg.output = "Total")
