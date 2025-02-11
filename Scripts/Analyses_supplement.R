  #load libraries
library(tidyverse)
library(vegan)
library(emmeans) #for posthoc comparisons of means
library(spaa) #for dist2list function 
library(cowplot)


#load community data
alldat_long = read.csv("Formatted_data/Seedbank_aboveground_merged_long.csv") %>%
  mutate(Site = as.factor(Site), Plot = as.factor(Plot), Treatment = as.factor(Treatment), Species = as.factor(Species),
         type = as.factor(type)) %>%
  mutate(site_byelev = factor(Site, levels = c("Blank", "Desert scrub", "Desert grassland", "Juniper savanna", "Ponderosa pine meadow", "Mixed conifer meadow"))) %>%
  mutate(trt_order = factor(Treatment, levels = c("Water Exclusion", "Control", "Water Addition", "Blank"))) %>%
  filter(Site != "Blank") %>% #filter blanks for now, but will need to figure that out
  #add elevations for sites 
  mutate(elevation = case_when(
    Site == "Desert scrub" ~ 1566,
    Site == "Desert grassland" ~ 1636,
    Site == "Juniper savanna" ~ 1930,
    Site == "Ponderosa pine meadow" ~ 2179,
    Site == "Mixed conifer meadow" ~ 2591
  )) %>%
  mutate(elevation = as.numeric(elevation)) %>%
  mutate(elevfact = as.factor(elevation)) %>%
  mutate(type = recode_factor(type, aboveground = "Aboveground", 
                              seedbank = "Seedbank"))


alldat_rel_wide = read.csv("Formatted_data/Relabun_sbandab_merged_wide.csv") %>%
  mutate(Site = as.factor(Site), Plot = as.factor(Plot), Treatment = as.factor(Treatment), type = as.factor(type)) %>%
  mutate(site_byelev = factor(Site, levels = c("Blank", "Desert scrub", "Desert grassland", "Juniper savanna", "Ponderosa pine meadow", "Mixed conifer meadow"))) %>%
  mutate(trt_order = factor(Treatment, levels = c("Water Exclusion", "Control", "Water Addition", "Blank")))%>%
  filter(Site != "Blank") %>%
  mutate(siteplot = paste0(Site, Plot)) %>%
  filter(siteplot != "Desert scrub1" ) %>% #no seedlings in Desert scrub 1 seedbank, so this will cause trouble with diversity and dissimilarity, remove here (but keep for abundance? which uses long format data)
  mutate(ID = paste(Site, Treatment, Plot, type, sep = "_")) %>% #create ID for use when making distance matrices and reshaping them
  select(-X, -siteplot) %>%
  #add elevations for sites 
  mutate(elevation = case_when(
    Site == "Desert scrub" ~ 1566,
    Site == "Desert grassland" ~ 1636,
    Site == "Juniper savanna" ~ 1930,
    Site == "Ponderosa pine meadow" ~ 2179,
    Site == "Mixed conifer meadow" ~ 2591
  )) %>%
  mutate(elevation = as.numeric(elevation))%>%
  mutate(elevfact = as.factor(elevation)) %>%
  mutate(type = recode_factor(type, aboveground = "Aboveground", 
                              seedbank = "Seedbank"))

alldat_raw_wide = read.csv("Formatted_data/totabun_sbandab_merged_wide.csv") %>%
  mutate(Site = as.factor(Site), Plot = as.factor(Plot), Treatment = as.factor(Treatment), type = as.factor(type))%>%
  mutate(site_byelev = factor(Site, levels = c("Blank", "Desert scrub", "Desert grassland", "Juniper savanna", "Ponderosa pine meadow", "Mixed conifer meadow"))) %>%
  mutate(trt_order = factor(Treatment, levels = c("Water Exclusion", "Control", "Water Addition", "Blank")))%>%
  filter(Site != "Blank") %>%
  mutate(siteplot = paste0(Site, Plot)) %>%
  filter(siteplot != "Desert scrub1" ) %>% #no seedlings in Desert scrub 1 seedbank, so this will cause trouble with diversity and dissimilarity, remove here (but keep for abundance? which uses long format data)
  mutate(ID = paste(Site, Treatment, Plot, type, sep = "_")) %>% #create ID for use when making distance matrices and reshaping them
  select(-X, -siteplot) %>%
  #add elevations for sites 
  mutate(elevation = case_when(
    Site == "Desert scrub" ~ 1566,
    Site == "Desert grassland" ~ 1636,
    Site == "Juniper savanna" ~ 1930,
    Site == "Ponderosa pine meadow" ~ 2179,
    Site == "Mixed conifer meadow" ~ 2591
  )) %>%
  mutate(elevation = as.numeric(elevation))%>%
  mutate(elevfact = as.factor(elevation)) %>%
  mutate(type = recode_factor(type, aboveground = "Aboveground", 
                              seedbank = "Seedbank"))




 
#load and format temperature data
#data from overwinter storage
winterfiles = list.files("Raw_data/raw ibutton data_overwinter", pattern= ".csv$")


winterdat = read_csv(paste0("Raw_data/raw ibutton data_overwinter/",winterfiles), #read all files in
                     id = "iid") %>%
  bind_rows() %>% #and bind by rows
  mutate(iid = c(substr(iid, nchar(iid)-7,nchar(iid)-6))) %>%#change file to unit number
  mutate(Date_Time = as.POSIXct(Date_Time, format ="%m/%d/%Y %H:%M")) %>%
  mutate(Date = as.Date(Date_Time)) %>%
  mutate(hour = round(hour(Date_Time),0)) %>%
  mutate(iid = as.factor(iid)) %>%
  select(-Unit)
write.csv(winterdat, file = "Formatted_data/winterandspringstoragetemps.csv")

#temperature data from growout in greenhouse
growoutfiles = list.files('Raw_data/raw hobo data_growout', pattern= ".csv$")
growoutfiles = growoutfiles[substr(growoutfiles, 1, 3) != "Air"]

growdat = read_csv(paste0("Raw_data/raw hobo data_growout/",growoutfiles), #read all files in
                   id = "iid", skip=2, col_names = c("num", "Date_Time", "Temp")) %>% #why isn't it skipping first row?
  bind_rows() %>% #and bind by rows
  mutate(iid = c(substr(iid, nchar(iid)-5,nchar(iid)-4))) %>%#change file to unit number
  mutate(Date_Time = as.POSIXct(Date_Time, format ="%m/%d/%Y %H:%M")) %>%
  mutate(Date = as.Date(Date_Time)) %>%
  mutate(hour = round(hour(Date_Time),0)) %>%
  mutate(iid = as.factor(iid)) 
 
write.csv(growdat, file = "Formatted_data/growouttemps.csv")

#format and plot overwinter temperatures

meanhourly = winterdat %>%
  group_by(Date, hour,iid) %>%
  summarize(meanhourlyT = mean(Value), sdhourlyT = sd(Value), n = length(Value))

meandaily_iid = meanhourly %>%
  group_by(Date, iid) %>%
  summarize(meandailyTid = mean(meanhourlyT ), nid = length(meanhourlyT ),
            maxdailyTid = max(meanhourlyT ),mindailyTid = min(meanhourlyT )) %>% #get mean, min and max for each button
  ungroup()
#summary(meandaily_iid)  

meandaily = meandaily_iid %>%
  group_by(Date) %>% #average across ibuttons
  summarize(meandailyT = mean(meandailyTid ), n = length(meandailyTid ),
            maxdailyT = max(maxdailyTid ),mindailyT = min(mindailyTid),
            sdmeandailyT = sd(meandailyTid, na.rm=T), sdmindailyT = sd(mindailyTid, na.rm=T),
            sdmaxdailyT = sd(maxdailyTid, na.rm=T)) %>%
  mutate(semean = meandailyT/sqrt(n), semin = mindailyT/sqrt(n),
         semax = maxdailyT/sqrt(n))

meanTplot = ggplot(meandaily, aes(x = Date, y= meandailyT)) + 
  geom_line(color = "black", linewidth = 0.5) + #geom_point() +
  #geom_errorbar(aes(ymin= meandailyT - semean, ymax= meandailyT + semean), width=0.1)+
  theme_classic(base_size = 20) +# theme(legend.text = element_text(size = 14))+
  labs(y="Temperature (°C)") #, title = "Overwinter storage temperatures")  


FigS1_overwintertemps = meanTplot +
  geom_line(data = meandaily, aes(x = Date, y= mindailyT), color = "blue")+
  geom_line(data = meandaily, aes(x = Date, y= maxdailyT), color = "red")
FigS1_overwintertemps #need to work on formatting
#ggsave("Plots/FigS1_overwintertemps.jpg", height = 8, width = 12)
 

#format and plot growout conditions

gmeanhourly = growdat %>%
  group_by(Date, hour,iid) %>%
  summarize(meanhourlyT = mean(Temp), sdhourlyT = sd(Temp), n = length(Temp))

gmeandaily_iid = gmeanhourly %>%
  group_by(Date, iid) %>%
  summarize(meandailyTid = mean(meanhourlyT ), nid = length(meanhourlyT ),
            maxdailyTid = max(meanhourlyT ),mindailyTid = min(meanhourlyT )) %>% #get mean, min and max for each button
  ungroup()

gmeandaily = gmeandaily_iid %>%
  group_by(Date) %>% #average across ibuttons
  summarize(meandailyT = mean(meandailyTid ), n = length(meandailyTid ),
            maxdailyT = max(maxdailyTid ),mindailyT = min(mindailyTid),
            sdmeandailyT = sd(meandailyTid, na.rm=T), sdmindailyT = sd(mindailyTid, na.rm=T),
            sdmaxdailyT = sd(maxdailyTid, na.rm=T)) %>%
  mutate(semean = meandailyT/sqrt(n), semin = mindailyT/sqrt(n),
         semax = maxdailyT/sqrt(n)) %>%
  drop_na()
summary(gmeandaily)
gmeanTplot = ggplot(gmeandaily, aes(x = Date, y= meandailyT)) + 
  geom_line(color = "black", size = 0.5) + #geom_point() +
  #geom_errorbar(aes(ymin= meandailyT - semean, ymax= meandailyT + semean), width=0.1)+
  theme_classic(base_size = 20) +# theme(legend.text = element_text(size = 14))+
  labs(y="Temperature (°C)") #, title = "Grow out temperatures")  


FigS2_growouttemps = gmeanTplot +
  geom_line(data = gmeandaily, aes(x = Date, y= mindailyT), color = "blue")+
  geom_line(data = gmeandaily, aes(x = Date, y= maxdailyT), color = "red")
FigS2_growouttemps #need to work on formatting
#ggsave("Plots/FigS2_growouttemps.jpg", height = 8, width = 12)

 

##Species richness and diversity calculation


#create dataframe that is just columns with relative abundances of each species
allsp_rel = alldat_rel_wide %>%
  #make rownames siteplot names
  column_to_rownames(var = "ID") %>%
  select(-Site, -Plot, -Treatment, -Year, -type, -site_byelev, -trt_order, -elevation, - elevfact) %>%
  #make NAs = 0 (NA means species wasn't in that plot/sample)
  mutate(across(everything(), ~replace_na(.x, 0))) 

#calculate richness and diversity metrics
alldat_rel_wide = alldat_rel_wide %>%
  filter(Site != "Blank") %>%
  mutate(rich = specnumber(allsp_rel)) %>%
  mutate(simp = diversity(allsp_rel, index = "simpson")) %>%
  mutate(shannon = diversity(allsp_rel, index = "shannon")) 
 
#create dataframes for aboveground vs seedbank communities, since area sampled affects richness and diversity (Vandvik et al. 2016)
alldat_rel_wide_ab = alldat_rel_wide %>%
  filter(type == "Aboveground") %>%
  droplevels()
summary(alldat_rel_wide_ab)

alldat_rel_wide_sb = alldat_rel_wide %>%
  filter(type == "Seedbank") %>%
  droplevels()
summary(alldat_rel_wide_sb)

#Simpson diversity
#aboveground
lm_simp_all_full_ab = lm(simp ~ Treatment*elevfact, data = alldat_rel_wide_ab)
summary(lm_simp_all_full_ab)                   
anova(lm_simp_all_full_ab)

#posthoc comparison
test(emmeans(lm_simp_all_full_ab, pairwise ~ Treatment, by=c("elevfact")) )
#aboveground:
#exclusion different from control at 1636 and 1930m, no differences with addition

#seedbank
lm_simp_all_full_sb = lm(simp ~ Treatment*elevfact, data = alldat_rel_wide_sb)
summary(lm_simp_all_full_sb)                   
anova(lm_simp_all_full_sb)

#posthoc comparison
test(emmeans(lm_simp_all_full_sb, pairwise ~ Treatment, by=c("elevfact")) )
#seedbank: no significant differences

 
#plot simpsons
simp_means = alldat_rel_wide %>%
  group_by(Site, site_byelev, Treatment, trt_order, type, elevation, elevfact) %>%
  summarize(meansimp = mean(simp), sdsimp = sd(simp), samplesize = n()) %>%
  mutate(sesimp = sdsimp/sqrt(samplesize))%>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    type == "Aboveground" & elevation == 1636 & Treatment == "Water Exclusion" ~ "*",
    type == "Aboveground" & elevation == 1930 & Treatment == "Water Exclusion" ~ "*"))%>%
  mutate(elevation_plotting = ifelse(type == "Aboveground", elevation - 30, elevation))

simp_plot_elevation = ggplot(simp_means, aes(x = elevation, y= meansimp, color = Treatment, 
                                             shape = type)) + 
  geom_point(position=position_dodge(10), size = 6) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meansimp - sesimp, ymax =  meansimp + sesimp), width = 0.2,
                position=position_dodge(10)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Simpson Diversity", color = "Treatment", linetype = "Community type", shape = "Community type") +
  geom_text(aes(label = siglabel), 
            position = position_dodge(), hjust = -1, size =12, show.legend = F) + 
  theme(#legend.justification = c(0.95, .05),legend.position = c(0.95,.05), 
    #     legend.key = element_rect(colour = NA, fill = NA),
    text = element_text(size = 20))

#### Fig. S3, simpsons diversity  ####
simp_plot_elevation
#ggsave("Plots/FigS3_Simpson diversity.jpg", height = 8, width = 12)



 

