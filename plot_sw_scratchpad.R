library(lubridate)
library(tidyverse)

apex_sw_subset <- apex_sw %>% filter(depth.cm <= 50, M %in% c(5:9)) %>% # filtering for growing season
  mutate(precip_event = ifelse(PRCP >= 2, 1, 0)) # ppt events >= 2mm

# dataframe of all precipitation events
precip.events <- apex_sw_subset %>% filter(precip_event == 1) %>%
  group_by(Y, depth.cm) %>%
  mutate(eventID_year = cumsum(precip_event)) # orders ppt events every year

# dataframe of single-day precipitation events
precip.events_single <- left_join(apex_sw_subset, precip.events) %>% 
  group_by(Y, depth.cm) %>%
  arrange(depth.cm, date) %>%
  filter(precip_event == 1 & 
           is.na(lag(eventID_year)) & # if the previous day isn't a ppt event
           is.na(lead(eventID_year))) %>% # if the next day isn't a ppt event
  mutate(single.event = "yes")

#### Calculating Volumetric Water Content Difference ###########################
## using apex_sw instead of apex_sw_subset b/c subset breaks for() loop on beginning and end dates of growing season

# new column that indicates single-day precipitation events
apex_sw <- left_join(apex_sw, precip.events_single) %>%
  select(-eventID_year) %>%
  mutate(single.event = ifelse(is.na(single.event), "no", "yes")) %>% # converting non-single ppt events to "no"
  arrange(depth.cm, date) # ordering and grouping data by soil depth then date

# initialize column for difference in VWC for single precip events
apex_sw$vwc.diff <- NA

for(i in 1:nrow(apex_sw)) {
  if(apex_sw$single.event[i] == "yes") {
    t1_date <- apex_sw$date[i]
    t0_date <- t1_date - days(1)
    t2_date <- t1_date + days(1)
    
    # Filter for vwc.daily values corresponding to t0 and t2 dates
    depth_value <- apex_sw$depth.cm[i]
    t0_vwc <- apex_sw$vwc.daily[apex_sw$date == t0_date & 
                               apex_sw$depth.cm == depth_value]
    t2_vwc <- apex_sw$vwc.daily[apex_sw$date == t2_date & 
                               apex_sw$depth.cm == depth_value]
    
    apex_sw$vwc.diff[i] <- t2_vwc - t0_vwc
  } else{
    apex_sw$vwc.diff[i] <- NA
  }
}

## adding new columns to subset
# filters out growing season and precip <= 2mm
apex_sw_subset.02 <- left_join(apex_sw_subset, apex_sw) %>%
  mutate(single.event = ifelse(is.na(single.event), "no", "yes")) # transforming column back to "yes/no" from "yes/NA"

#### Calculting Total Water Content ############################################
# filtering for days with single precipitation events
apex_sw_subset.03 <- apex_sw_subset.02 %>% filter(!is.na(vwc.diff))

# initialize wc.cm column with NA values (blank)
apex_sw_subset.03$wc.cm <- NA

# creates new column for conversion of volumetric water content (%) to water content (cm)
for(i in 2:nrow(apex_sw_subset.03)){
  if(apex_sw_subset.03$depth.cm[i] != 2.5) { # these soil depths are 10cm in depth
    apex_sw_subset.03$wc.cm[i] <- apex_sw_subset.03$vwc.diff[i] * 10
  } else if(apex_sw_subset.03$depth.cm[i] == 2.5) { # only depth that is 5cm
    apex_sw_subset.03$wc.cm[i] <- apex_sw_subset.03$vwc.diff[i] * 5
  }
}

# create new dataframe that shows the sum of the water content for each single precipitation event
apex_wc <- apex_sw_subset.03 %>%
  group_by(Pasture, date, PRCP) %>%
  summarize(daily.wc.cm = sum(wc.cm, na.rm = TRUE)) # summing the wc across all depths for each date*pasture
