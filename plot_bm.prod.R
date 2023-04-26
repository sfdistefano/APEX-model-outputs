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
  select(date, Y, ID, CPNM, STD, STL, GZSDkg.ha, GZSLkg.ha) %>%
  mutate(STD.kgha = STD*1000, STL.kgha = STL*1000) %>%
  mutate(STLD = STD.kgha + STL.kgha + GZSDkg.ha + GZSLkg.ha) %>%
  group_by(date, Y, ID) %>%
  summarize(STDL = sum(STLD))

apex_stdl_CPNM <- apex_growing.season %>% 
  select(date, Y, ID, CPNM, STD, STL, GZSDkg.ha, GZSLkg.ha) %>%
  mutate(STD.kgha = STD*1000, STL.kgha = STL*1000) %>%
  mutate(STLD = STD.kgha + STL.kgha + GZSDkg.ha + GZSLkg.ha) %>%
  group_by(date, Y, ID, CPNM) %>%
  summarize(STDL = sum(STLD))

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
plot_bm.prod <- function(veg.output,
                         pasture = F, 
                         ecosite = F, treatment = F, es_treatment = F,
                         past.name, func.group, es, trt) {
  
  if(veg.output == "CPNM" & pasture == T) {
    
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
      summarize(Biomass_fg = sum(Biomass_wt)) # sum of 1 func.group across ecosites w/in pasture
    
    # SE calculated from plots; some pastures may have 1 plot, SE = 0
    field02 <- biomass.plot_info %>%
      filter(Pasture == past.name & CPNM2 == func.group) %>%
      mutate(Biomass_wt = kgPerHa*Proportion) %>%
      group_by(date, Y, Pasture, CPNM2) %>%
      summarize(se = sd(Biomass_wt)/sqrt(length(Biomass_wt)))
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_fg)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_fg - se, 
                        ymax = field01$Biomass_fg + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Pasture", past.name, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "CPNM" & ecosite == T) {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es, CPNM == func.group) %>%
      group_by(date, Y, Ecosite, CPNM) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM))
    
    
    stdl <- merge(x = apex_stdl_CPNM, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es & CPNM == func.group) %>%
      group_by(date, Y, Ecosite, CPNM) %>%
      summarize(STDL = mean(STDL))
    
    # mean calculated across pastures w/in ecosite
    field01 <- biomass.pasture %>% 
      filter(Ecosite == es & CPNM2 == func.group) %>%
      group_by(date, Y, Ecosite, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass))
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(Ecosite == es & CPNM2 == func.group) %>%
      group_by(date, Y, Ecosite, Plot, CPNM2) %>%
      summarize(Biomass_fg = mean(kgPerHa)) %>%
      group_by(date, Y, Ecosite, CPNM2) %>%
      summarize(se = sd(Biomass_fg, na.rm = T)/sqrt(length(Biomass_fg))) # SE across plots w/in ecosite
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)"))  +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_fg)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_fg - se, 
                        ymax = field01$Biomass_fg + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Ecosite", es, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "CPNM" & treatment == T) {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      filter(graze.trt == trt, CPNM == func.group) %>%
      group_by(date, Y, graze.trt, CPNM) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM))
    
    stdl <- merge(x = apex_stdl_CPNM, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      filter(graze.trt == trt & CPNM == func.group) %>%
      group_by(date, Y, graze.trt, CPNM) %>%
      summarize(STDL = mean(STDL)) %>%
      filter(graze.trt == trt & CPNM == func.group)
    
    # mean calculated across pastures w/in grazing treatment
    field01 <- biomass.pasture %>% 
      filter(graze.trt == trt & CPNM2 == func.group) %>%
      group_by(date, Y, graze.trt, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass)) # mean of each func.group w/in treatment
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(graze.trt == trt & CPNM2 == func.group) %>%
      group_by(date, Y, graze.trt, Plot, CPNM2) %>%
      summarize(Biomass_fg = mean(kgPerHa)) %>%
      group_by(date, Y, graze.trt, CPNM2) %>%
      summarize(se = sd(Biomass_fg, na.rm = T)/sqrt(length(Biomass_fg))) # SE across plots w/in treatment
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)"))  +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_fg)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_fg - se, 
                        ymax = field01$Biomass_fg + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Grazing", trt, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "CPNM" & es_treatment == T) {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es, graze.trt == trt, CPNM == func.group) %>%
      group_by(date, Y, Ecosite, graze.trt, CPNM) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM))
    
    stdl <- merge(x = apex_stdl_CPNM, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es, graze.trt == trt, CPNM == func.group) %>%
      group_by(date, Y, Ecosite, graze.trt, CPNM) %>%
      summarize(STDL = mean(STDL))
    
    # mean calculated across pastures w/in ecosite*grazing treatment
    field01 <- biomass.pasture %>% 
      filter(Ecosite == es & graze.trt == trt & CPNM2 == func.group) %>%
      group_by(date, Y, Ecosite, graze.trt, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass)) # mean of each func.group w/in treatment
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(Ecosite == es & graze.trt == trt & CPNM2 == func.group) %>%
      group_by(date, Y, Ecosite, graze.trt, Plot, CPNM2) %>%
      summarize(Biomass_fg = mean(kgPerHa)) %>%
      group_by(date, Y, Ecosite, graze.trt, CPNM2) %>%
      summarize(se = sd(Biomass_fg, na.rm = T)/sqrt(length(Biomass_fg))) # SE across plots w/in ecosite*treatment
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)"))  +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_fg)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_fg - se, 
                        ymax = field01$Biomass_fg + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Ecosite", es,"- Grazing", trt, "CPNM = ", func.group))
    
    return(plot)
    
  }
  
  if(veg.output == "Total" & pasture == T){
    
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
    
    # mean calculated across pastures 
    field01 <- biomass.pasture %>% 
      filter(Pasture == past.name) %>%
      group_by(date, Y, Pasture, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass)) %>%  # mean of each func.group w/in treatment
      group_by(date, Y, Pasture) %>%
      summarize(Biomass_tot = sum(Biomass_fg, na.rm = T)) # sum of all func.group
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(Pasture == past.name) %>%
      group_by(date, Y, Pasture, Plot, CPNM2) %>%
      summarize(Biomass_fg = mean(kgPerHa)) %>%
      group_by(date, Y, Pasture, Plot) %>%
      summarize(Biomass_tot = sum(Biomass_fg)) %>%
      group_by(date, Y, Pasture) %>%
      summarize(se = sd(Biomass_tot, na.rm = T)/sqrt(length(Biomass_tot))) # SE across plots w/in pasture
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_tot)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_tot - se, 
                        ymax = field01$Biomass_tot + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Pasture", past.name, "TOTAL biomass"))
    
    return(plot)
  }
  
  if(veg.output == "Total" & ecosite == T) {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es) %>%
      group_by(date, Y, Ecosite, Pasture) %>%
      summarize(A_DDM = sum(A_DDM),
                cum_DDM = sum(cum_DDM)) %>%
      group_by(date, Y, Ecosite) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM))
    
    stdl <- merge(x = apex_stdl, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es) %>%
      group_by(date, Y, Ecosite) %>%
      summarize(STDL = mean(STDL))
    
    # mean calculated across pastures w/in grazing treatment
    field01 <- biomass.pasture %>% 
      filter(Ecosite == es) %>%
      group_by(date, Y, Ecosite, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass)) %>%  # mean of each func.group w/in treatment
      group_by(date, Y, Ecosite) %>%
      summarize(Biomass_tot = sum(Biomass_fg, na.rm = T)) # sum of all func.group
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(Ecosite == es) %>%
      group_by(date, Y, Ecosite, Plot, CPNM2) %>%
      summarize(Biomass_fg = mean(kgPerHa)) %>%
      group_by(date, Y, Ecosite, Plot) %>%
      summarize(Biomass_tot = sum(Biomass_fg)) %>%
      group_by(date, Y, Ecosite) %>%
      summarize(se = sd(Biomass_tot, na.rm = T)/sqrt(length(Biomass_tot))) # SE across plots w/in treatment
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_tot)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_tot - se,
                        ymax = field01$Biomass_tot + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Ecosite", es, "TOTAL biomass"))
    
    return(plot)
    
  }
  
  if(veg.output == "Total" & treatment == T) {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      filter(graze.trt == trt) %>%
      group_by(date, Y, graze.trt, Pasture) %>%
      summarize(A_DDM = sum(A_DDM),
                cum_DDM = sum(cum_DDM)) %>%
      group_by(date, Y, graze.trt) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM))
    
    stdl <- merge(x = apex_stdl, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      filter(graze.trt == trt) %>%
      group_by(date, Y, graze.trt) %>%
      summarize(STDL = mean(STDL))
    
    # mean calculated across pastures w/in grazing treatment
    field01 <- biomass.pasture %>% 
      filter(graze.trt == trt) %>%
      group_by(date, Y, graze.trt, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass)) %>%  # mean of each func.group w/in treatment
      group_by(date, Y, graze.trt) %>%
      summarize(Biomass_tot = sum(Biomass_fg, na.rm = T)) # sum of all func.group
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(graze.trt == trt) %>%
      group_by(date, Y, graze.trt, Plot, CPNM2) %>%
      summarize(kgPerHa = mean(kgPerHa)) %>%
      group_by(date, Y, graze.trt, Plot) %>%
      summarize(Biomass_tot = sum(kgPerHa)) %>%
      group_by(date, Y, graze.trt) %>%
      summarize(se = sd(Biomass_tot, na.rm = T)/sqrt(length(Biomass_tot))) # SE across plots w/in treatment
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_tot)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_tot - se,
                        ymax = field01$Biomass_tot + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Grazing", trt, "TOTAL biomass"))
    
    return(plot)
    
  }
  
  if(veg.output == "Total" & es_treatment == T) {
    
    cumsum <- merge(apex_cumsum_CPNM, pastID_ecosite,
                    by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es & graze.trt == trt) %>%
      group_by(date, Y, Ecosite, graze.trt, Pasture) %>%
      summarize(A_DDM = sum(A_DDM), # summing all func.group of each pasture
                cum_DDM = sum(cum_DDM)) %>%
      group_by(date, Y, Ecosite, graze.trt) %>%
      summarize(A_DDM = mean(A_DDM),
                cum_DDM = mean(cum_DDM))
    
    
    stdl <- merge(x = apex_stdl, y = pastID_ecosite,
                  by.x = "ID", by.y = "PastureID") %>%
      filter(Ecosite == es & graze.trt == trt) %>%
      group_by(date, Y, Ecosite, graze.trt == trt) %>%
      summarize(STDL = mean(STDL))
    
    # mean calculated first at pasture level then grazing trt*ecosite
    field01 <- biomass.pasture %>% 
      filter(Ecosite == es & graze.trt == trt) %>%
      mutate(Biomass_wt = Biomass*Proportion) %>%
      group_by(date, Y, Ecosite, graze.trt, CPNM2) %>%
      summarize(Biomass_fg = mean(Biomass_wt)) %>%
      group_by(date, Y, Ecosite, graze.trt) %>%
      summarize(Biomass_tot = sum(Biomass_fg))
    
    # SE calculated from plots
    field02 <- biomass.plot_info %>%
      filter(Ecosite == es & graze.trt == trt) %>%
      mutate(Biomass_wt = kgPerHa*Proportion) %>%
      group_by(date, Y, Ecosite, graze.trt, Plot, CPNM2) %>%
      summarize(Biomass_fg = mean(kgPerHa)) %>%
      group_by(date, Y, Ecosite, graze.trt, Plot) %>%
      summarize(Biomass_tot = sum(Biomass_fg)) %>%
      group_by(date, Y, Ecosite, graze.trt) %>%
      summarize(se = sd(Biomass_tot)/sqrt(length(Biomass_tot)))
    
    plot <- ggplot() +
      geom_line(data = cumsum, aes(x = date, y = cum_DDM, color = "Biomass Production (above&below)")) +
      geom_line(data = stdl, aes(x = date, y = STDL, color = "Ungrazed Biomass (above)")) +
      geom_point(data = field01, aes(x = date, y = Biomass_tot)) +
      geom_errorbar(data = field02,
                    aes(x = date, ymin = field01$Biomass_tot - se,
                        ymax = field01$Biomass_tot + se),
                    width = 10) +
      facet_wrap(.~Y, scales = "free_x") +
      scale_x_date(date_breaks = "1 month",date_labels = "%m") +
      theme_bw() +
      ylab("Biomass (kg/ha)") +
      xlab("Month of the Year (May - Oct)") +
      labs(color = "Biomass") +
      ggtitle(paste("Ecosite", es,"- Grazing", trt, "TOTAL biomass"))
    
    return(plot)
    
  }
}


plot_bm.prod(veg.output = "CPNM", pasture = T,
             func.group = "BOBU", past.name = "7NW")
plot_bm.prod(veg.output = "CPNM", treatment = T,
             func.group = "BOBU", trt = "TRM")
plot_bm.prod(veg.output = "CPNM", ecosite = T,
             func.group = "BOBU", es = "Sandy")
plot_bm.prod(veg.output = "CPNM", es_treatment = T,
             func.group = "BOBU", trt = "TRM", es = "Sandy")


plot_bm.prod(veg.output = "Total", pasture = T, 
             past.name = "15E")
plot_bm.prod(veg.output = "Total", ecosite = T, 
             es = "Sandy")
plot_bm.prod(veg.output = "Total", treatment = T,
             trt = "TRM")
plot_bm.prod(veg.output = "Total", es_treatment = T,
             trt = "TRM", es = "Loamy")
