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
  filter(precip_event == 1 & 
           is.na(lag(eventID_year)) & # if the previous day isn't a ppt event
           is.na(lead(eventID_year))) %>% # if the next day isn't a ppt event
  mutate(single.event = "yes")

# new column that indicates single-day precipitation events
apex_sw_subset.02 <- left_join(apex_sw_subset, precip.events_single) %>%
  select(-eventID_year) %>%
  mutate(single.event = ifelse(is.na(single.event), "no", "yes")) # converting non-single ppt events to "no"

# creates new column for the difference of T2 - T0 (water balance)
for (i in 2:nrow(apex_sw_subset.02)) {
  if (apex_sw_subset.02$single.event[i] == "yes") { # day of (T1)
    apex_sw_subset.02$vwc.diff[i] <- apex_sw_subset.02$vwc.daily[i+1] - # day after (T2)
      apex_sw_subset.02$vwc.daily[i-1] # day before (T0)
  } else {
    apex_sw_subset.02$vwc.diff[i] <- NA
  }
}

# subsetting for days with single precipitation events
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
  group_by(Pasture, date) %>%
  summarize(daily.wc.cm = sum(wc.cm, na.rm = TRUE)) # summing the wc across all depths for each date*pasture
