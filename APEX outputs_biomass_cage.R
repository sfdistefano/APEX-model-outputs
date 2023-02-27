library(lubridate)
library(tidyverse)

# setwd("C:/Users/sfper/Dropbox/USDA-ARS_FC/APEX1905_New/APEX1905_New/APEX1905_CPER_SCENARIO") # CHANGE THIS TO WHERE YOU HAVE SAVED THE BELOW CSV FILE WITH DATA
setwd("D:/APEX data and scripts/Data")

##### IMPORT DATA ##############################################################
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
  mutate(Month = 8, Day = 5) %>% # adding date info, when collected in field
  group_by(Month, Day, Year, Pasture, Ecosite, Treatment, FGCode) %>%
  summarize(Biomass = mean(kgPerHa), Stdev = sd(kgPerHa)) %>% # summarized by pasture
  rename(CPNM = FGCode) %>%
  filter(!(Pasture %in% c("1W", "28N", "32W"))) # pastures not included in CARM or TRM

## APEX output files
apexsad1 = read.delim("D:/APEX model/APEX1905_New_34_peak biomass/CONUNN_AGM.sad",
                      sep = "", dec = ".", skip = 8)
# additional runs of APEX with different parameters
apexsad2 = read.delim("D:/APEX model/APEX1905_New_34_peak biomass/CONUNN_AGM.sad",
                      sep = "", dec = ".", skip = 8)

apexsad3 = read.delim("D:/APEX model/APEX1905_New_34_peak biomass/CONUNN_AGM.sad",
                      sep = "", dec = ".", skip = 8)

#### DATA FORMATTING ###########################################################
# formatting date columns
biomass.pasture$Date2 <- paste(biomass.pasture$Year, 
                               biomass.pasture$Month, 
                               biomass.pasture$Day,sep="-") %>% 
  ymd() # changing to date format

apexsad1$Date2<-with(apexsad1,paste(Y,M,D,sep="-")) %>% ymd()
# apexsad2$Date2<-with(apexsad2,paste(Y,M,D,sep="-")) %>% ymd()
# apexsad3$Date2<-with(apexsad3,paste(Y,M,D,sep="-")) %>% ymd()

#### Add pasture and adjust functional group names #############################
# adding pasture ID to match with APEX output files
dat2 <- merge(biomass.pasture, pastID_ecosite, 
              by = c("Pasture", "Ecosite"), all.x = TRUE) 

# fixing functional groups names so they match APEX names
dat2$CPNM2 <- dat2$CPNM %>% 
  gsub(pattern = "C3PG", replacement = "CSPG") %>%
  gsub(pattern = "FORB", replacement = "FRB3") %>%
  gsub(pattern = "SS", replacement = "SSHB")

## Adding ecosite column and summarizing data (day*ecosite*func.group)
# summarizing to 3 significant figures, signif()
apex.dat1 <- merge(apexsad1, pastID_ecosite,
                   by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
  group_by(Date2,Ecosite,CPNM) %>%
  summarize(STL = signif(mean(STL), 3), # STL = standing live
            STD = signif(mean(STD), 3), # STD = standing dead
            GZSL = mean(GZSLkg.ha), # grazing intake of standing live
            GZSD = mean(GZSDkg.ha)) %>% # grazing intake of standing dead
  mutate(bm = STL + STD + GZSL + GZSD)

apex.dat2 <- merge(apexsad2, pastID_ecosite,
                   by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
  group_by(Date2,Ecosite,CPNM) %>%
  summarize(STL = signif(mean(STL), 3), # STL = standing live
            STD = signif(mean(STD), 3), # STD = standing dead
            GZSL = mean(GZSLkg.ha), # grazing intake of standing live
            GZSD = mean(GZSDkg.ha)) %>% # grazing intake of standing dead
  mutate(bm = STL + STD + GZSL + GZSD)

apex.dat3 <- merge(apexsad3, pastID_ecosite,
                   by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
  group_by(Date2,Ecosite,CPNM) %>%
  summarize(STL = signif(mean(STL), 3), # STL = standing live
            STD = signif(mean(STD), 3), # STD = standing dead
            GZSL = mean(GZSLkg.ha), # grazing intake of standing live
            GZSD = mean(GZSDkg.ha)) %>% # grazing intake of standing dead
  mutate(bm = STL + STD + GZSL + GZSD)


# cage biomass data (day*ecosite*func.group)
dat3 <- dat2 %>% group_by(Date2, Ecosite, CPNM2) %>%
  summarize(bm = mean(Biomass), se = sd(Biomass)/sqrt(length(Biomass)))


##### APEX vs. FIELD DATA ######################################################

#### PLOTTING TOTAL BIOMASS FOR ONE ECOLOGICAL SITE ############################
# STL = standing live
# STD = standing deadapex.dat1_total

apex.dat1_total <- apex.dat1 %>% group_by(Date2, Ecosite) %>%
  summarize(tot.bm = sum(bm))

# apex.dat2_total <- apex.dat2 %>% group_by(Date2, Ecosite) %>%
#   summarize(tot.bm = sum(bm))
# 
# apex.dat3_total <- apex.dat3 %>% group_by(Date2, Ecosite) %>%
#   summarize(tot.bm = sum(bm))

# cage biomass data from field (date*ecosite*sum of all func.group)
dat4 <- dat3 %>% group_by(Date2, Ecosite) %>%
  summarize(tot.bm = sum(bm), se = sd(bm)/sqrt(length(bm)))

## plotting total biomass (day*ecosite[sum of func.group])
plot.biomass01 <- function(dat, eco, run1, run2, run3, run4){
  title <- paste("TOTAL biomass for",eco, "(APEX vs field data)", sep = " ")

  dat <- dat %>% filter(Ecosite == eco)
  napexsad1 <- apex.dat1_total %>% filter(Ecosite == eco & Date2 > "2013-01-01")
  # napexsad2 <- apex.dat2_total %>% filter(Ecosite == eco & Date2 > "2013-01-01")
  # napexsad3 <- apex.dat3_total %>% filter(Ecosite == eco & Date2 > "2013-01-01")

  # needed for added lines for start of each year
  date_range <- which(napexsad1$Date2 %in% as.Date(
    c("2013-01-01",
      "2014-01-01", "2015-01-01", "2016-01-01", "2017-01-01",
      "2018-01-01", "2019-01-01", "2020-01-01", "2021-01-01")
  ))

  ggplot() +
    geom_line(data=napexsad1, aes(x=Date2,y=tot.bm, color = run1),
              linetype="solid", size = 1) +
    # geom_line(data=napexsad2, aes(x=Date2,y=tot.bm, color = run2),
    #           linetype="solid", size = 1) +
    # geom_line(data=napexsad3, aes(x=Date2,y=tot.bm, color = run3),
    #           linetype="solid", size = 1) +
    # geom_line(data=napexsad4, aes(x=Date2,y=tot.bm, color = run4),
    #           linetype="solid", size = 1) +
    geom_vline(xintercept = as.numeric(napexsad1$Date2[date_range]),
               linetype = "dashed") +
    geom_point(data = dat, aes(x=Date2, y=tot.bm), color = 'black', size = 2) +
    geom_errorbar(data=dat,
                  aes(x=Date2, ymin=tot.bm-se, ymax=tot.bm+se),
                  width=30, colour="black", size = 1) +
    ylab("Biomass (kg per ha)") +
    xlab("Time") +
    # ylim(0,3750) +
    ggtitle(title) +
    theme_classic() +
    theme(text = element_text(family = "serif", size = 15)) +
    scale_color_brewer(palette = "Set1") +
    scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
    labs(color = "Number of Subareas")

}

# plot.biomass01(dat = dat4, eco = "Loamy",
#                run1 = "20 subareas - constant",
#                run2 = "20 subareas - field",
#                run3 = "34 subareas - constant",
#                run4 = "34 subareas - field")

#### PLOTTING BIOMASS ACROSS ECOLOGICAL SITE FOR ONE FUNCTIONAL GROUP ##########
# function for plotting across all pastures of an ecological site
plot.biomass02 <- function(dat, func.group, eco, parm, leg.title,
                           default, run1, run2, run3)
  {
  title <- paste(parm, "Biomass for",eco, func.group, sep = " ")
  
  dat2 <- dat %>% filter(Ecosite == eco & CPNM2 == func.group)
  napexsad1 <- apex.dat1 %>% filter(Ecosite == eco & CPNM == func.group)
  napexsad2 <- apex.dat2 %>% filter(Ecosite == eco & CPNM == func.group)
  napexsad3 <- apex.dat3 %>% filter(Ecosite == eco & CPNM == func.group)
  napexsad4 <- apex.dat4 %>% filter(Ecosite == eco & CPNM == func.group)

  date_range <- which(napexsad1$Date2 %in% as.Date(
    c("2013-01-01",
      "2014-01-01", "2015-01-01", "2016-01-01", "2017-01-01",
      "2018-01-01", "2019-01-01", "2020-01-01", "2021-01-01")
  ))
  
  ggplot() +
    geom_point(data = dat2, aes(x=Date2, y=bm), color = 'black') +
    geom_errorbar(data=dat2,
                  aes(x=Date2, ymin=pmax(0,(bm-se)), ymax=(bm+se)),
                  width=30,colour="black") +
    geom_line(data=napexsad1, aes(x = Date2, y = bm, color = default),
              linetype="solid", size=1) +
    geom_line(data=napexsad2, aes(x = Date2, y = bm, color = run1),
              linetype="dashed", size=1) +
    geom_line(data=napexsad3, aes(x = Date2, y = bm, color = run2),
              linetype="dashed", size=1) +
    geom_line(data=napexsad4, aes(x = Date2, y = bm, color = run3),
              linetype="dashed", size=1) +
    geom_vline(xintercept = as.numeric(napexsad1$Date2[date_range]),
               linetype = "dashed", size = 1) +
    ylab("Biomass (kg per ha)") +
    xlab("Time") +
    # ylim(0,ymax) +
    ggtitle(title) +
    theme_classic() +
    theme(text = element_text(family = "serif", size = 15)) +
    scale_color_brewer(palette = "Set1") +
    scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
    labs(color = leg.title)
  
}

# Salt Flats doesn't work for 20 pastures bc not added in - need to to total biomass
# plot.biomass02(dat = dat3, func.group = "CSPG", eco = "Sandy",
#                leg.title = "Number of Subareas",
#                parm = "APEX vs. Cage Biomass:",
#                default = "20 subareas - constant",
#                run1 = "20 subareas - field",
#                run2 = "34 subareas - constant",
#                run3 = "34 subareas - field")

#### PLOTTING MULTIPLE FUNCTIONAL GROUPS FOR ONE APEX RUN ######################
# eco <- "Salt Flats"
# 
# napexsad1 <- apex.dat1 %>% filter(Ecosite == eco) %>% filter(CPNM == "BOBU")
# napexsad2 <- apex.dat1 %>% filter(Ecosite == eco) %>% filter(CPNM == "CSPG")
# napexsad3 <- apex.dat1 %>% filter(Ecosite == eco) %>% filter(CPNM == "WSPG")
# 
# ggplot() +
#   geom_line(data=napexsad1, aes(x=Date2,y=(STL+STD)*1000, color = "BOBU"),
#             linetype="solid", size=1) +
#   geom_line(data=napexsad2, aes(x=Date2,y=(STL+STD)*1000, color = "CSPG"),
#             linetype="solid", size=1) +
#   geom_line(data=napexsad3, aes(x=Date2,y=(STL+STD)*1000, color = "WSPG"),
#             linetype="solid", size=1) +
#   geom_vline(xintercept = as.Date(apexsad1$Date2[date_range]),
#              linetype = "dashed", size = 1) +
#   ylab("Biomass (kg per ha)") +
#   xlab("Time") +
#   ylim(0,800) +
#   theme_classic() +
#   theme(text = element_text(family = "serif", size = 15)) +
#   scale_color_brewer(palette = "Set1") +
#   scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
#   ggtitle("SALT FLATS: biomass by functional group")

#### PLOTTING BIOMASS FOR INDIVIDUAL PASTURE ###################################

## formatting data
dat2_sub <- dat2 %>% group_by(Date2, Pasture, PastureID, Ecosite) %>% 
  summarize(tot.bm = sum(Biomass)) %>%
  group_by(Pasture, Date2, Ecosite) %>%
  summarize(sub.bm = mean(tot.bm))

dat2_sub7 <- dat2_sub %>% 
  filter(Pasture == "7NW" & Ecosite == "Sandy") %>%
  mutate(w.bm = sub.bm/6)
dat2_sub27 <- dat2_sub %>% 
  filter(Pasture == "7NW" & Ecosite == "Loamy") %>%
  mutate(w.bm = sub.bm/2)
dat2_sub28 <- dat2_sub %>% 
  filter(Pasture == "7NW" & Ecosite == "Salt Flats") %>%
  mutate(w.bm = sub.bm/3)

dat2_past <- rbind(dat2_sub7, dat2_sub27, dat2_sub28) %>%
  group_by(Pasture, Date2) %>%
  summarize(past.bm = sum(w.bm), se = sd(w.bm)/sqrt(length(w.bm)))



# apex.dat1_past <- merge(apexsad1, pastID_ecosite,
#                         by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
#   mutate(bm = ((STL+STD)*1000) + GZSLkg.ha + GZSDkg.ha) %>%
#   group_by(Date2, Pasture) %>%
#   summarize(past.bm = sum(bm))


## 34 pastures
apex.dat2_sub <- merge(apexsad1, pastID_ecosite,
                       by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
  mutate(bm = ((STL+STD)*1000) + GZSLkg.ha + GZSDkg.ha) %>%
  group_by(Date2, Pasture, Ecosite) %>%
  summarize(sub.bm = sum(bm))

apex.dat2_sub7 <- apex.dat2_sub %>% 
  filter(Pasture == "7NW" & Ecosite == "Sandy") %>%
  mutate(w.bm = sub.bm/6)
apex.dat2_sub27 <- apex.dat2_sub %>% 
  filter(Pasture == "7NW" & Ecosite == "Loamy") %>%
  mutate(w.bm = sub.bm/2)
apex.dat2_sub28 <- apex.dat2_sub %>% 
  filter(Pasture == "7NW" & Ecosite == "Salt Flats") %>%
  mutate(w.bm = sub.bm/3)

apex.dat2_past <- rbind(apex.dat2_sub7, apex.dat2_sub27, apex.dat2_sub28) %>%
  group_by(Pasture, Date2) %>%
  summarize(past.bm = sum(w.bm))


# apex.dat3_sub <- merge(apexsad3, pastID_ecosite,
#                        by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
#   mutate(bm = ((STL+STD)*1000) + GZSLkg.ha + GZSDkg.ha) %>%
#   group_by(Date2, Pasture, Ecosite) %>%
#   summarize(sub.bm = sum(bm))
# 
# apex.dat3_sub7 <- apex.dat3_sub %>% 
#   filter(Pasture == "7NW" & Ecosite == "Sandy") %>%
#   mutate(w.bm = sub.bm/6)
# apex.dat3_sub27 <- apex.dat3_sub %>% 
#   filter(Pasture == "7NW" & Ecosite == "Loamy") %>%
#   mutate(w.bm = sub.bm/2)
# apex.dat3_sub28 <- apex.dat3_sub %>% 
#   filter(Pasture == "7NW" & Ecosite == "Salt Flats") %>%
#   mutate(w.bm = sub.bm/3)
# 
# apex.dat3_past <- rbind(apex.dat3_sub7, apex.dat3_sub27, apex.dat3_sub28) %>%
#   group_by(Pasture, Date2) %>%
#   summarize(past.bm = sum(w.bm))


## plotting individual pasture
date_range <- which(apex.dat2_past$Date2 %in% as.Date(
  c("2013-01-01",
    "2014-01-01", "2015-01-01", "2016-01-01", "2017-01-01",
    "2018-01-01", "2019-01-01", "2020-01-01", "2021-01-01")
))

plot.biomass03 <- function(dat, pastname, parm, leg.title,
                          default, run1, run2, run3
                          ){
  title <- paste("TOTAL biomass for Pasture", pastname, "(APEX vs field data)", sep = " ")
  
  dat2 <- dat %>% filter(Pasture == pastname)
  
  # napexsad1 <- apex.dat1_past %>% filter(Pasture == pastname & Date2 > "2013-01-01")
  napexsad2 <- apex.dat2_past %>% filter(Pasture == pastname & Date2 > "2013-01-01")
  # napexsad3 <- apex.dat3_past %>% filter(Pasture == pastname & Date2 > "2013-01-01")

  ggplot() +
    # geom_line(data=napexsad1, aes(x = Date2, y = past.bm, color = run1),
    #           linetype="solid", size=1) +
    geom_line(data=napexsad2, aes(x = Date2, y = past.bm, color = run2),
              size=1) +
    # geom_line(data=napexsad3, aes(x = Date2, y = past.bm, color = run3),
    #           size=1) +
    geom_point(data = dat2, aes(x=Date2, y=past.bm), color = 'black') +
    geom_errorbar(data=dat2,
                  aes(x=Date2, ymin=pmax(0,(past.bm-se)), ymax=(past.bm+se)),
                  width=30,colour="black", size = 1) +
    geom_vline(xintercept = as.numeric(apex.dat2_past$Date2[date_range]),
               linetype = "dashed", size = 1) +
    ylab("Total Biomass (kg per ha)") +
    xlab("Time") +
    # ylim(0,ymax) +
    theme_classic() +
    theme(text = element_text(family = "serif", size = 20)) +
    scale_color_brewer(palette = "Set1") +
    scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
    labs(color = leg.title)
  
}

plot.biomass03(dat = dat2_past, pastname = "7NW",
               leg.title = "# of Subareas-Basal Area",
               parm = "APEX vs. Cage Biomass:",
               run1 = "20-constant",
               run2 = "34-constant",
               run3 = "34-peak biomass")
