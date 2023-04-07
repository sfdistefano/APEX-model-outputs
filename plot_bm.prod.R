library(tidyverse)
library(lubridate)
library(zoo)

setwd("D:/APEX data and scripts/Data/")

## APEX output
apexsad = read.delim("D:/APEX model/APEX1905_New/APEX1905_New/CONUNN_AGM.sad",
                     sep = "", dec = ".", skip = 8) %>%
  mutate(date = paste(Y,M,D, sep = "-")) %>%
  mutate(date = ymd(date))
## Filtering for growing season (Jun - Oct)
apex_growing.season <- apexsad %>% 
  filter(M >= 6 & M <= 10) %>% # June - October
  filter(Y >= 2013) # starting year of CARM

## Calculating cumulative sum --> biomass production
# cumulative sum for EACH CPNM
apex_cumsum_CPNM <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, A_DDM) %>%
  group_by(ID, Y, CPNM) %>%
  mutate(cum_DDM = cumsum(A_DDM)*1000)

# cumulative sum of ALL CPNM (i.e., total biomass production)
apex_cumsum_TOTAL <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, A_DDM) %>%
  group_by(ID, Y) %>%
  mutate(cum_DDM = cumsum(A_DDM)*1000,
         date = ymd(date))

## reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite.csv")

## original cage biomass used (day*pasture*plot*functional group)
biomass.plot <- read.csv("CARM_Biomass_cln_attr2021-12-01_wide.csv") %>% 
  select(-Forage, -HERB) %>%
  pivot_longer(cols = BOBU:SD, names_to = "FGCode", values_to = "kgPerHa") %>%
  group_by(Year, Treatment, Pasture, Ecosite, Plot, FGCode) %>%
  summarize(kgPerHa = mean(kgPerHa))%>% # summarize by plot of each pasture
  filter(!(Plot == 5 & Pasture == "18S")) %>% # prescribed burn plots
  filter(!(Plot == 6 & Pasture == "18S")) %>%
  filter(!(Plot == 5 & Pasture == "19N")) %>%
  filter(!(Plot == 6 & Pasture == "19N")) 

# biomass data is summarized to day*pasture*functional group
biomass.pasture <-  biomass.plot %>% 
  mutate(Month = 8, Day = 5,
         date = ymd(paste(Year, Month, Day, sep = "-"))
  ) %>% # adding date info, when collected in field
  group_by(date, Year, Pasture, Ecosite, Treatment, FGCode) %>%
  summarize(Biomass = mean(kgPerHa), 
            se = sd(kgPerHa)/sqrt(length(kgPerHa))) %>% # summarized by pasture
  rename(CPNM = FGCode, Y = Year) %>%
  filter(!(Pasture %in% c("1W", "28N", "32W"))) %>% # pastures not included in CARM or TRM
  filter(!(CPNM == "SD")) # removing standing dead

## adding pasture ID to match with APEX output files
biomass.pasture <- merge(biomass.pasture, pastID_ecosite, 
                         by = c("Pasture", "Ecosite"), all.x = TRUE) 

## fixing functional groups names so they match APEX names
biomass.pasture$CPNM2 <- biomass.pasture$CPNM %>% 
  gsub(pattern = "C3PG", replacement = "CSPG") %>%
  gsub(pattern = "FORB", replacement = "FRB3") %>%
  gsub(pattern = "SS", replacement = "SSHB")

## function for plotting biomass
plot_bm.prod <- function(past.output, past.name, ecosite, veg.output, func.group) {
  
  if(veg.output == "CPNM" & past.output == "Pasture") {
    
    cumsum <- merge(x = apex_cumsum_CPNM, y = pastID_ecosite, 
                    by.x = "ID", by.y = "PastureID") %>% 
      filter(Pasture == past.name & CPNM == func.group)
    
    field <- biomass.pasture %>% filter(Pasture == past.name & CPNM2 == func.group)
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM), color = "red") +
      geom_point(data = field, aes(x = date, y = Biomass)) +
      geom_errorbar(data = field,
                    aes(x = date, ymin = Biomass - se, ymax = Biomass + se),
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
    
    field <- biomass.pasture %>%
      group_by(date, Y, Ecosite, CPNM2) %>%
      summarize(Biomass_mean = mean(Biomass), se = sd(Biomass)/sqrt(length(Biomass))) %>%
      filter(Ecosite == ecosite, CPNM2 == func.group)
    
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM), color = "red") +
      geom_point(data = field, aes(x = date, y = Biomass_mean)) +
      geom_errorbar(data = field,
                    aes(x = date, ymin = Biomass_mean - se, ymax = Biomass_mean + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass Production (kg/ha)") +
      xlab("Month of the Year (Jun - Oct)") +
      ggtitle(paste("Ecosite", ecosite, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "Total" & past.output == "Pasture"){
    
    cumsum <- merge(x = apex_cumsum_TOTAL, y = pastID_ecosite, 
                    by.x = "ID", by.y = "PastureID") %>% 
      filter(Pasture == past.name)
    
    field <- biomass.pasture %>%
      group_by(date, Y, Pasture) %>%
      summarize(Biomass_total = sum(Biomass), se = sd(Biomass)/sqrt(length(Biomass))) %>%
      filter(Pasture == past.name)
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM), color = "red") +
      geom_point(data = field, aes(x = date, y = Biomass_total)) +
      geom_errorbar(data = field,
                    aes(x = date, ymin = Biomass_total - se, ymax = Biomass_total + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass Production (kg/ha)") +
      xlab("Month of the Year (Jun - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Pasture", past.name, "TOTAL biomass"))
    
    return(plot)
  }
  
  if(veg.output == "Total" & past.output == "Ecosite") {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      group_by(date, Y, Ecosite) %>%
      summarize(A_DDM = sum(A_DDM),
                cum_DDM = sum(cum_DDM)) %>%
      filter(Ecosite == ecosite)
    
    field <- biomass.pasture %>%
      group_by(date, Y, Ecosite, Pasture) %>%
      summarize(Biomass_total = sum(Biomass), 
                se = sd(Biomass)/sqrt(length(Biomass))) %>%
      group_by(date, Y, Ecosite) %>%
      summarize(Biomass_mean = mean(Biomass_total),
                se = sd(Biomass_total)/sqrt(length(Biomass_total))) %>%
      filter(Ecosite == ecosite)
    
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM), color = "red") +
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
      ggtitle(paste("Ecosite", ecosite, "TOTAL biomass"))
    
    return(plot)
    
  }
}

plot_bm.prod(past.output = "Pasture", past.name = "19N", 
             veg.output = "CPNM", func.group = "BOBU")

plot_bm.prod(past.output = "Pasture", past.name = "19N", 
             veg.output = "Total")


plot_bm.prod(past.output = "Ecosite", ecosite = "Sandy", 
             veg.output = "CPNM", func.group = "BOBU")

plot_bm.prod(past.output = "Ecosite", ecosite = "Sandy",
             veg.output = "Total",)
