apex.dat2_past.eco <- merge(apexsad1, pastID_ecosite,
                   by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
  group_by(Date2,Pasture,Ecosite,CPNM) %>%
  summarize(STL = signif(mean(STL), 3), # STL = standing live
            STD = signif(mean(STD), 3), # STD = standing dead
            GZSL = mean(GZSLkg.ha), # grazing intake of standing live
            GZSD = mean(GZSDkg.ha)) %>% # grazing intake of standing dead
  mutate(bm = ((STL + STD)*1000) + GZSL + GZSD)

# apex.dat3_past.eco <- merge(apexsad3, pastID_ecosite,
#                    by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
#   group_by(Date2,Pasture,Ecosite,CPNM) %>%
#   summarize(STL = signif(mean(STL), 3), # STL = standing live
#             STD = signif(mean(STD), 3), # STD = standing dead
#             GZSL = mean(GZSLkg.ha), # grazing intake of standing live
#             GZSD = mean(GZSDkg.ha)) %>% # grazing intake of standing dead
#   mutate(bm = ((STL + STD)*1000) + GZSL + GZSD)

## 20 subareas
# apex.dat1_past.fg <- merge(apexsad1, pastID_ecosite,
#                            by.x = "ID", by.y = "PastureID",all.x = TRUE) %>%
#   mutate(bm = ((STL+STD)*1000) + GZSLkg.ha + GZSDkg.ha) %>%
#   group_by(Date2, Pasture, Ecosite, CPNM) %>%
#   summarize(fg.bm = sum(bm))
# 
# apex.dat1_past_wspg <- apex.dat1_past.fg %>% filter(Pasture == "7NW") %>%
#   filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
#   group_by(Pasture, Date2) %>%
#   summarize(past.bm = sum(fg.bm))

## field data
dat2_past.fg <- dat2 %>%
  group_by(Date2, Pasture, Ecosite, CPNM2) %>%
  summarize(fg.bm = sum(Biomass))

dat2_sub7_wspg <-  dat2_past.fg %>%
  filter(CPNM2 == "WSPG" | CPNM2 == "BOBU") %>%
  filter(Pasture == "7NW" & Ecosite == "Sandy") %>%
  group_by(Pasture, Date2) %>%
  summarize(w.bm = sum(fg.bm)/6)
dat2_sub27_wspg <-  dat2_past.fg %>%
  filter(CPNM2 == "WSPG" | CPNM2 == "BOBU") %>%
  filter(Pasture == "7NW" & Ecosite == "Loamy") %>%
  group_by(Pasture, Date2) %>%
  summarize(w.bm = sum(fg.bm)/2)  
dat2_sub28_wspg <-  dat2_past.fg %>%
  filter(CPNM2 == "WSPG" | CPNM2 == "BOBU") %>%
  filter(Pasture == "7NW" & Ecosite == "Salt Flats") %>%
  group_by(Pasture, Date2) %>%
  summarize(w.bm = sum(fg.bm)/3)  

dat2_past_wspg <- rbind(dat2_sub7_wspg, dat2_sub27_wspg, dat2_sub28_wspg) %>%
  group_by(Pasture, Date2) %>%
  summarize(past.bm = sum(w.bm), se = sd(w.bm)/sqrt(length(w.bm)))

## 34 subareas
apex.dat2_sub7_wspg <- apex.dat2_past.eco %>% 
  filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
  filter(Pasture == "7NW" & Ecosite == "Sandy") %>%
  group_by(Pasture, Date2) %>%
  summarize(w.bm = sum(bm)/6)
apex.dat2_sub27_wspg <- apex.dat2_past.eco %>% 
  filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
  filter(Pasture == "7NW" & Ecosite == "Loamy") %>%
  group_by(Pasture, Date2) %>%
  summarize(w.bm = sum(bm)/2)
apex.dat2_sub28_wspg <- apex.dat2_past.eco %>% 
  filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
  filter(Pasture == "7NW" & Ecosite == "Salt Flats") %>%
  group_by(Pasture, Date2) %>%
  summarize(w.bm = sum(bm)/3)


apex.dat2_past_wspg <- rbind(apex.dat2_sub7_wspg, apex.dat2_sub27_wspg, apex.dat2_sub28_wspg) %>%
  group_by(Pasture, Date2) %>%
  summarize(past.bm = sum(w.bm))


# apex.dat3_sub7_wspg <- apex.dat3_past.eco %>% 
#   filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
#   filter(Pasture == "7NW" & Ecosite == "Sandy") %>%
#   group_by(Pasture, Date2) %>%
#   summarize(w.bm = sum(bm)/6)
# apex.dat3_sub27_wspg <- apex.dat3_past.eco %>% 
#   filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
#   filter(Pasture == "7NW" & Ecosite == "Loamy") %>%
#   group_by(Pasture, Date2) %>%
#   summarize(w.bm = sum(bm)/2)
# apex.dat3_sub28_wspg <- apex.dat3_past.eco %>% 
#   filter(CPNM == "WSPG" | CPNM == "BOBU") %>%
#   filter(Pasture == "7NW" & Ecosite == "Salt Flats") %>%
#   group_by(Pasture, Date2) %>%
#   summarize(w.bm = sum(bm)/3)
# 
# 
# apex.dat3_past_wspg <- rbind(apex.dat3_sub7_wspg, apex.dat3_sub27_wspg, apex.dat3_sub28_wspg) %>%
#   group_by(Pasture, Date2) %>%
#   summarize(past.bm = sum(w.bm))

## plotting WSPG + BOBU
leg.title = "# of Subareas-Basal Area"
parm = "APEX vs. Cage Biomass:"
run1 = "20-constant"
run2 = "34-constant"
run3 = "34-peak biomass"


title <- paste("Biomass of Functional Group WSPG for Pasture 7NW (APEX vs field data)", sep = " ")


# napexcag1 <- apex.dat1_past_wspg %>% filter(Date2 > "2013-01-01")
napexcag2 <- apex.dat2_past_wspg %>% filter(Date2 > "2013-01-01")
# napexcag3 <- apex.dat3_past_wspg %>% filter(Date2 > "2013-01-01")

ggplot() +
  geom_line(data=napexcag1, aes(x = Date2, y = past.bm, color = run1),
            linetype="solid", size=1) +
  geom_line(data=napexcag2, aes(x = Date2, y = past.bm, color = run2),
            size=1) +
  geom_line(data=napexcag3, aes(x = Date2, y = past.bm, color = run3),
            size=1) +
  geom_point(data = dat2_past_wspg, aes(x=Date2, y=past.bm), color = 'black') +
  geom_errorbar(data = dat2_past_wspg,
                aes(x=Date2, ymin=pmax(0,(past.bm-se)), ymax=(past.bm+se)),
                width=30,colour="black", size = 1) +
  geom_vline(xintercept = as.numeric(apex.dat2_past$Date2[date_range]),
             linetype = "dashed", size = 1) +
  ylab("Warm-Season Biomass (kg per ha)") +
  xlab("Time") +
  # ylim(0,ymax) +
  theme_classic() +
  theme(text = element_text(family = "serif", size = 20)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_date(date_breaks = "1 year",date_labels = "%Y")
