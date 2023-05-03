library(tidyverse)
library(lubridate)

setwd("D:/APEX data and scripts/Data/")

## Reference dataframe for pasture, pastureID, and ecological site
pastID_ecosite <- read.csv("PastureID_ecosite.csv")

## APEX output
apexsad <- read.delim("D:/APEX model/APEX1905_New - Copy/APEX1905_New/CONUNN_AGM.sad",
                     sep = "", dec = ".", skip = 8) %>%
  select(Y,M,D,ID,PRCP,X0.5CM:X75.85CM) %>% # import desired columns
  unique() # SW values repeat so only need unique values per date

# summarizing for soil water data (volumetric water content)
apex_sw <- apexsad %>%
  mutate(date = paste(Y, M, D, sep="-")) %>%
  mutate(date = ymd(date)) %>% # adding date column in date format
  pivot_longer(cols = X0.5CM:X75.85CM, names_to = "depth.range",
               values_to = "vwc.daily") %>% # making a longer dataframe that matches CPER data format
  mutate(depth.range = gsub("X","", depth.range), # removing extra characters
         vwc.daily = vwc.daily*100) %>% # converting to %
  filter(date >= "2018-01-01") # date start of CPER data

## Creating reference dataframe for soil depth
depth.range <- unique(apex_sw$depth.range)
depth.cm <- c(2.5, 10, 20, 30,40,50,60,70,80)

depth.ref <- data.frame(depth.range, depth.cm)

# adding in depths and pasture names
apex_sw <- merge(apex_sw, depth.ref, by = "depth.range") %>%
  merge(pastID_ecosite, by.x = "ID", by.y = "PastureID")

## CPER soil moisture field data (2018-2022)
cper_sw <- read.csv("19N-15E_soilmoisture_2018-2022.csv") %>%
  mutate(date = mdy(date),
         Y = year(date))

## CPER precipitation
cper_ppt <- read.csv("CPER PPT/CPER daily climate data.csv") %>%
  mutate(date = mdy(TimeStamp)) %>%
  filter(date >= "2018-01-01") %>%
  na.omit()

cper_ppt_month <- cper_ppt %>% group_by(MONTH, YEAR) %>% 
  summarize(mon.precip = sum(RainTotal_mm)) %>%
  mutate(date = ym(paste(YEAR,MONTH)))

## Plotting volumetric water content
plot_sw <- function(past_name, depth){
  
  # filtering for desired pasture and soil depth (cm)
  apex <- apex_sw %>% filter(Pasture == past_name & depth.cm == depth) 
  max.apex <- max(apex$vwc.daily)
  
  max.ppt <- max(cper_ppt$RainTotal_mm)
  
  coeff <- max.apex/max.ppt
  
  ppt <- cper_ppt %>% mutate(ppt_scale = RainTotal_mm*coeff)
  
  cper <- cper_sw %>% filter(site == past_name & depth.cm == depth) %>%
    na.omit(vwc.daily) %>% # removing missing data
    mutate(vwc.daily_scale = vwc.daily*(max.apex/max(vwc.daily)))
  
  plot <- ggplot() +
    geom_bar(data = ppt, aes(x = date, y = ppt_scale, 
                             fill = "Daily Precip"), alpha = 0.75, 
             stat = 'identity') +
    geom_line(data = apex, aes(x = date, y = vwc.daily, color = "APEX")) +
    geom_line(data = cper, aes(x = date, y = vwc.daily_scale, color = "CPER"))  +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    xlab("Time") +
    scale_y_continuous("Volumetric Water Content (%)", 
                       # limits = c(0,30),
                       sec.axis = sec_axis(~./coeff, name = "Daily Precipitation (mm)")) +
    scale_color_brewer(type = "div", palette = "Dark2") +
    scale_fill_manual(values = "#0072B2") + # added for fill color
    labs(color = "Data Source", fill = "Data Type") + # added for fill color
    ggtitle(paste("Pasture",past_name, ", soil depth =", depth, "cm")) +
    theme_bw()
  
  
  return(plot)
}

plot_sw(past_name = "19N", depth = 10)
