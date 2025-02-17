#load libraries
library(tidyverse)
library(vegan)
library(emmeans) #for posthoc comparisons of means
library(spaa) #for dist2list function 
library(ggplot2) #for plotting
library(cowplot) #for plotting




#load data
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
summary(alldat_long)

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
summary(alldat_rel_wide)

alldat_wide = read.csv("Formatted_data/totabun_sbandab_merged_wide.csv") %>%
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
summary(alldat_wide)

#Total abundance
#calculate total abundance for each plot
abun_tots = alldat_long %>%
  group_by(Site, Treatment, type, Plot, site_byelev, trt_order, elevation, elevfact) %>%
  summarize(totabun = sum(totabun)) 

#linear models for above ground and seedbank separately
#aboveground
lm_abun_ab = lm(totabun ~ elevfact*Treatment, data = subset(abun_tots, type == "Aboveground"))
summary(lm_abun_ab)                   
anova(lm_abun_ab)

#Post hoc contrasts
emmeans(lm_abun_ab, pairwise ~ Treatment|elevfact) 
#at 1566, Water exclusion sig diff from control
#at 1636: no sig treatment diff
#at 1930: water exclusion diff from control 
#at 2179: water exclusion diff from control 
#at 2591: water exclusion diff from control     

lm_abun_sb = lm(totabun ~ elevfact * Treatment, data = subset(abun_tots, type == "Seedbank"))
summary(lm_abun_sb)                   
anova(lm_abun_sb)

#post hoc contrasts
emmeans(lm_abun_sb, pairwise ~ Treatment|elevfact) 
#only significant difference is at 2179, water exclusion different from control

#plot abundances, Figure 1
abun_means = abun_tots %>%
  group_by(Site, site_byelev, Treatment, trt_order, type, elevation, elevfact) %>%
  summarize(meanabun = mean(totabun), sdabun = sd(totabun), samplesize = n()) %>%
  mutate(seabun = sdabun/sqrt(samplesize))

abun_means_ab = abun_means %>%
  filter(type == "Aboveground") %>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    elevation  == 1566 & Treatment == "Water Exclusion" ~ "*",
    elevation  == 1930 & Treatment == "Water Exclusion" ~ "*",
    elevation  == 2179 & Treatment == "Water Exclusion" ~ "*",
    elevation  == 2591 & Treatment == "Water Exclusion" ~ "*"))

abun_means_sb = abun_means %>%
  filter(type == "Seedbank") %>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    elevation  == 2179 & Treatment == "Water Exclusion" ~ "*"))

abund_elevation_ab = ggplot(abun_means_ab, aes(x = elevation, y= meanabun, group = Treatment, color = Treatment)) + 
  geom_point(position=position_dodge(25), size = 6, shape = 1, stroke =1) + #geom_line(position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = meanabun - seabun, ymax = meanabun + seabun), width = 0.2,
                position=position_dodge(25)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4", "dodgerblue4", "firebrick4" )) +
  labs(x = "Elevation (m, asl)", y = "Mean total abundance\n (Mean cover)", color = "Treatment",
       title = "Above ground community") +
  geom_text(aes(label = siglabel), 
            position = position_dodge(25), hjust = -0.5, size =12, show.legend = F) + 
  theme(legend.justification = c(0.05, .95),legend.position = c(0.05,.95), 
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 20)) +
  xlim(c(1560,2610))

abund_elevation_sb = ggplot(abun_means_sb, aes(x = elevation, y= meanabun, group = Treatment, color = Treatment)) + 
  geom_point(position=position_dodge(25), size = 6, shape = 18) + #geom_line(position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = meanabun - seabun, ymax = meanabun + seabun), width = 0.2,
                position=position_dodge(25)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4", "dodgerblue4", "firebrick4" )) +
  labs(x = "Elevation (m, asl)", y = "Mean total abundance\n (Mean count)", color = "Treatment",
       title = "Seed bank community") +
  geom_text(aes(label = siglabel), 
            position = position_dodge(25), hjust = -0.5, size =12, show.legend = F) + 
  theme(legend.justification = c(0.05, .95),legend.position = c(0.05,.95), 
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 20)) +
  xlim(c(1560,2610))

#### Figure 1: abundances ####
plot_grid(abund_elevation_ab, abund_elevation_sb  , align = "hv", labels = c("A.", "B."), label_size=20)
#ggsave("./Plots/Fig1_Mean abundance.jpg", height = 8, width = 12)



#Relative abundance of plant functional types
#create PFT data frame.  Sum relative abundances by functional type
pftdat = alldat_long %>%
  filter(Site != "Blank") %>%
  group_by(Site, site_byelev,elevation ,elevfact, Plot, Treatment, trt_order, PFT, type) %>%
  summarize(pft_sum = sum(totabun, na.rm=T), pft_sum_rel = sum(relabun, na.rm=T)) %>%   
  ungroup() %>%
  mutate(trt_type = as.factor(paste(trt_order, type, sep="_")))

summary(pftdat)
str(pftdat)

#need to reshape relative abundances into wide form for MANOVA
pft_wide = pftdat %>%
  pivot_wider(
    id_cols = c(Site, site_byelev, elevfact, elevation, Plot, Treatment, trt_order, type),
    names_from = PFT, 
    values_from = pft_sum_rel) %>%
  mutate(A = replace_na(A, 0),
         PF = replace_na(PF),
         PG4 = replace_na(PG4, 0),
         S = replace_na(S, 0),
         PG3 = replace_na(PG3, 0))   

dep_vars = with(pft_wide, cbind(A, PG4, PF, PG3))

relabun_manova_mod = manova(dep_vars ~ Treatment*elevfact*type, data = pft_wide)
summary(relabun_manova_mod, test = c("Wilks")) 

#Since we have a significant manova, we can do oneway ANOVAs to test for pairwise differences
#Annuals 
annualdat = subset(pftdat, PFT == "A" ) %>%
  droplevels()

annual_mod = lm(pft_sum_rel ~ Treatment*elevfact*type, data = annualdat)
anova(annual_mod)
#plot(annual_mod)

test(emmeans(annual_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Seedbank")) ) #no significant differences
test(emmeans(annual_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Aboveground")) ) 
#Water exclusion different from control at 1566, no other differences

#perennial forbs
perennialdat = subset(pftdat, PFT == "PF" ) %>%
  droplevels()

perennial_mod = lm(pft_sum_rel ~ Treatment*elevfact*type, data = perennialdat)
anova(perennial_mod)
#plot(perennial_mod)

test(emmeans(perennial_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Seedbank")) ) 
#Water exclusion different from control at 1566 and 1930, no other differences
test(emmeans(perennial_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Aboveground")) ) 
#no  differences


#perennial C3 Grasses
c3dat = subset(pftdat, PFT == "PG3" ) %>%
  droplevels()

c3_mod = lm(pft_sum_rel ~ Treatment*elevfact*type, data = c3dat)
anova(c3_mod)
#plot(c3_mod)

test(emmeans(c3_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Seedbank")) ) 
#no  differences
test(emmeans(c3_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Aboveground")) ) 
#no  differences


#perennial C4 Grasses
c4dat = subset(pftdat, PFT == "PG4" ) %>%
  droplevels()

c4_mod = lm(pft_sum_rel ~ Treatment*elevfact*type, data = c4dat)
anova(c4_mod)
#plot(c4_mod)

test(emmeans(c4_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Seedbank")) ) 
#Exclusion different from control at 1930m
test(emmeans(c4_mod, pairwise ~ Treatment*type|elevfact, at = list(type = "Aboveground")) ) 
#Exclusion marginally different at 1636 and 1930


#plot relative abundance of functional types
meanpftrel = pftdat %>%
  group_by(Site, site_byelev, elevation, elevfact,Treatment, trt_order, PFT, 
           type, trt_type) %>%
  summarize(meanPFTrel = mean(pft_sum_rel, na.rm=T), sdPFTrel = sd(pft_sum_rel, na.rm=T), samplesize = n()) %>%
  mutate(sePFTrel = sdPFTrel/sqrt(samplesize))%>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant", siglabel_trans = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    PFT == "A" & type == "Aboveground" & elevation == 1566 & Treatment == "Water Exclusion" ~ "*",
    PFT == "PF" & type == "Seedbank" & elevation == 1566  & Treatment == "Water Exclusion" ~ "*",
    PFT == "PF" & type == "Seedbank" & elevation == 1930  & Treatment == "Water Exclusion" ~ "*",
    PFT == "PG4" &type == "Aboveground" & elevation == 1636  & Treatment == "Water Exclusion" ~ "+",
    PFT == "PG4" &type == "Aboveground" & elevation == 1930  & Treatment == "Water Exclusion" ~ "+",
    PFT == "PG4" &type == "Seedbank" & elevation == 1930 & Treatment == "Water Exclusion" ~ "*"))  %>%
  mutate(elevation_plotting = ifelse(type == "Aboveground", elevation - 20, elevation))

## panel by relative abundance ##
pft_elevation_Annuals = ggplot(subset(meanpftrel, PFT == "A"), aes(x = elevation_plotting, y= meanPFTrel, color = Treatment,
                                                                   shape = type)) + 
  geom_point(position=position_dodge(40), size = 4) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanPFTrel - sePFTrel, ymax =  meanPFTrel + sePFTrel), width = 0.2,
                position=position_dodge(40)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Relative abundance", color = "Treatment", 
       shape = "Community type") + ggtitle("Annual forbs")+
  geom_text(aes(label = siglabel), 
            position = position_dodge(), hjust = .8, size =8, show.legend = F) + 
  theme(text = element_text(size = 14))+
  ylim(0,1) + xlim(1500,2700)

pft_elevation_Perennials = ggplot(subset(meanpftrel, PFT == "PF"), aes(x = elevation_plotting, y= meanPFTrel, color = Treatment,
                                                                       shape = type)) + 
  geom_point(position=position_dodge(40), size = 4) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanPFTrel - sePFTrel, ymax =  meanPFTrel + sePFTrel), width = 0.2,
                position=position_dodge(40)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Relative abundance", color = "Treatment", 
       shape = "Community type") + ggtitle("Perennial forbs")+
  geom_text(aes(label = siglabel), 
            position = position_dodge(), hjust = .8, size =8, show.legend = F) + 
  theme(legend.justification = c(0.95, .99),legend.position = c(0.95,.99), 
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 14))   + 
  ylim(0,1)+ xlim(1500,2700)


pft_elevation_C3s = ggplot(subset(meanpftrel, PFT == "PG3"), aes(x = elevation_plotting, y= meanPFTrel, color = Treatment,
                                                                 shape = type)) + 
  geom_point(position=position_dodge(40), size = 4) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanPFTrel - sePFTrel, ymax =  meanPFTrel + sePFTrel), width = 0.2,
                position=position_dodge(40)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Relative abundance", color = "Treatment", 
       shape = "Community type") + ggtitle("Perennial C3 grasses")+
  geom_text(aes(label = siglabel), 
            position = position_dodge(), hjust = .8, size =8, show.legend = F) + 
  theme(legend.position="bottom", text = element_text(size = 14))+ guides(shape="none")+
  ylim(0,1)+ xlim(1500,2700)

pft_elevation_C4s = ggplot(subset(meanpftrel, PFT == "PG4"), aes(x = elevation_plotting, y= meanPFTrel, color = Treatment,
                                                                 shape = type)) + 
  geom_point(position=position_dodge(40), size = 4) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanPFTrel - sePFTrel, ymax =  meanPFTrel + sePFTrel), width = 0.2,
                position=position_dodge(40)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Relative abundance", color = "Treatment", 
       shape = "Community type") + ggtitle("Perennial C4 grasses")+
  geom_text(aes(label = siglabel), 
            position = position_dodge(), hjust = .8, size =8, show.legend = F) + 
  theme(legend.position="bottom", text = element_text(size = 14))+ guides(color = "none") + 
  ylim(0,1)+ xlim(1500,2700)


#### Figure 2: Relative abundance of plant functional types ####

plot_grid(pft_elevation_Annuals + theme(legend.position = "none"),
          pft_elevation_Perennials+ theme(legend.position = "none") ,
          pft_elevation_C3s, 
          pft_elevation_C4s, 
          nrow=2, labels = c("A.", "B.", "C.", "D."), label_size=20)
#ggsave("./Plots/Fig2_PFTrelabun_byPFT.jpg", height = 8, width = 12)


#Species richness and diversity calculation
#create dataframe that is just columns with relative abundances of each species
allsp_rel = alldat_rel_wide %>%
  #make rownames siteplot names
  column_to_rownames(var = "ID") %>%
  select(-Site, -Plot, -Treatment, -Year, -type, -site_byelev, -trt_order, -elevation, - elevfact) %>%
  #make NAs = 0 (NA means species wasn't in that plot/sample)
  mutate(across(everything(), ~replace_na(.x, 0))) 

#calculate richness and diversity
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

## Species richness  ##
#aboveground
lm_rich_all_ab = lm(rich ~ Treatment*elevfact, data = alldat_rel_wide_ab)
summary(lm_rich_all_ab)                   
anova(lm_rich_all_ab)

#posthoc comparison
test(emmeans(lm_rich_all_ab, pairwise ~ Treatment, by=c("elevfact")) )
#aboveground:
#exclusion sig different from control at 1930

#seedbank
lm_rich_all_sb = lm(rich ~ Treatment*elevfact, data = alldat_rel_wide_sb)
summary(lm_rich_all_sb)                   
anova(lm_rich_all_sb)

#posthoc comparison
test(emmeans(lm_rich_all_sb, pairwise ~ Treatment, by=c("elevfact")) )
#Seedbank:
#no significant differences


#Species diversity
#aboveground
lm_shannon_all_ab = lm(shannon ~ Treatment*elevfact, data = alldat_rel_wide_ab)
summary(lm_shannon_all_ab)                   
anova(lm_shannon_all_ab)

#posthoc comparison
test(emmeans(lm_shannon_all_ab, pairwise ~ Treatment, by=c("elevfact")) )
#aboveground:
#exclusion sig different from control at 1636m and 1930, no differences with addition

#seedbank
lm_shannon_all_sb = lm(shannon ~ Treatment*elevfact, data = alldat_rel_wide_sb)
summary(lm_shannon_all_sb)                   
anova(lm_shannon_all_sb)

#posthoc comparison
test(emmeans(lm_shannon_all_sb, pairwise ~ Treatment, by=c("elevfact")) )
#seedbank: marg sig difference between control and addition at 1636


#Plot richness and diversity
#richness
rich_means = alldat_rel_wide %>%
  group_by(Site, site_byelev, Treatment, trt_order, type, elevation, elevfact) %>%
  summarize(meanrich = mean(rich), sdrich = sd(rich), samplesize = n()) %>%
  mutate(serich = sdrich/sqrt(samplesize)) %>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    type == "Aboveground" & elevation  == 1930 & Treatment == "Water Exclusion" ~ "*"))%>%
  mutate(elevation_plotting = ifelse(type == "Aboveground", elevation - 30, elevation))

rich_plot_elevation = ggplot(rich_means, aes(x = elevation_plotting, y= meanrich, color = Treatment, 
                                             shape = type)) + 
  geom_point(position=position_dodge(10), size = 6) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanrich - serich, ymax =  meanrich + serich), width = 0.2,
                position=position_dodge(10)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Species richness", color = "Treatment", linetype = "Community type", shape = "Community type") +
  geom_text(aes(label = siglabel), hjust = 1, size =10, show.legend = F) + #,   position = position_dodge() 
  theme(legend.justification = c(0.05, .99),legend.position = c(0.05,.99), 
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 20))

#shannon diversity
shannon_means = alldat_rel_wide %>%
  group_by(Site, site_byelev, Treatment, trt_order, type, elevation, elevfact) %>%
  summarize(meanshannon = mean(shannon), sdshannon = sd(shannon), samplesize = n()) %>%
  mutate(seshannon = sdshannon/sqrt(samplesize))%>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    type == "Aboveground" &  elevation  == 1636  & Treatment == "Water Exclusion" ~ "*",
    type == "Seedbank" &  elevation  == 1636  & Treatment == "Water Addition" ~ "*",
    type == "Aboveground" & elevation  == 1930 & Treatment == "Water Exclusion" ~ "*"))  %>%
  mutate(elevation_plotting = ifelse(type == "Aboveground", elevation - 30, elevation))

shannon_plot_elevation = ggplot(shannon_means, aes(x = elevation_plotting, y= meanshannon, color = Treatment, 
                                                   shape = type)) + 
  geom_point(position=position_dodge(10), size = 6) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanshannon - seshannon, ymax =  meanshannon + seshannon), width = 0.2,
                position=position_dodge(10)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Shannon Diversity", color = "Treatment", linetype = "Community type", shape = "Community type") +
  geom_text(aes(label = siglabel), hjust = 1, size =12, show.legend = F) + #,   position = position_dodge()
  theme(legend.justification = c(0.05, .95),legend.position = c(0.05,.95), 
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 20))

#### Figure 3: Richness and diversity ####
plot_grid(rich_plot_elevation, shannon_plot_elevation+ theme(legend.position = "none")   , labels = c("A.", "B."), label_size=20)
#ggsave("./Plots/Fig3_richness_diversity.jpg", height = 8, width = 15)


#Species composition - NMDS
#create grouping data
all_groups = alldat_rel_wide %>%
  select(Site, Treatment, type, site_byelev, trt_order, elevation, elevfact)

set.seed(123)
all_NMS_dist = metaMDS(as.matrix(allsp_rel), distance = "bray", k=3, maxit= 999, trymax=500) #didn't converge well with k=2

goodness(all_NMS_dist)
stressplot(all_NMS_dist)

nmds_dist_scores = as.data.frame(scores(all_NMS_dist, "sites")) %>%
  mutate(site = all_groups$Site, treatment = all_groups$Treatment, type = all_groups$type,
         site_byelev = all_groups$site_byelev, trt_order = all_groups$trt_order,
         elevation = all_groups$elevation)

#plot it
nmds_dist_plot_trt = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y= NMDS2, shape = type)) +
  geom_point(aes(colour = trt_order), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "Treatment", shape= "Community", fill = "Elevation", linetype= "Elevation") 

#add convex hulls
site_ds <- nmds_dist_scores[nmds_dist_scores$site == "Desert scrub",][chull(nmds_dist_scores[nmds_dist_scores$site ==  "Desert scrub", c("NMDS1", "NMDS2")]), ]  
site_dg <- nmds_dist_scores[nmds_dist_scores$site == "Desert grassland", ][chull(nmds_dist_scores[nmds_dist_scores$site  == "Desert grassland", c("NMDS1", "NMDS2")]), ]  
site_js <- nmds_dist_scores[nmds_dist_scores$site == "Juniper savanna", ][chull(nmds_dist_scores[nmds_dist_scores$site == "Juniper savanna", c("NMDS1", "NMDS2")]), ]  
site_ppm <- nmds_dist_scores[nmds_dist_scores$site == "Ponderosa pine meadow", ][chull(nmds_dist_scores[nmds_dist_scores$site == "Ponderosa pine meadow", c("NMDS1", "NMDS2")]), ]  # hull values for site A
site_mcm <- nmds_dist_scores[nmds_dist_scores$site == "Mixed conifer meadow",][chull(nmds_dist_scores[nmds_dist_scores$site == "Mixed conifer meadow", c("NMDS1", "NMDS2")]), ]  

hull.data <- rbind(site_ds, site_dg, site_js, site_ppm, site_mcm)  %>%
  mutate(elevfact = as.factor(elevation))

elevs = sort(unique(all_groups$elevfact))

nmds_all = nmds_dist_plot_trt + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=elevfact, linetype = elevfact),color = "black", fill = NA) 

#### Figure 4: NMDS of species composition ####        
nmds_all +
  annotate("text", x=-1, y=-2, label=paste0("Desert scrubland, ",elevs[1],"m"),
           color="black", size = 5)+
  annotate("text", x=-0.5, y=-0.25, label=paste0("Desert grassland, ",elevs[2],"m"),
           color="black", size = 5)+
  annotate("text", x=1, y=1.5, label=paste0("Juniper savanna, ",elevs[3],"m"),
           color="black", size = 5)+
  annotate("text", x=1.75, y=0.8, label=paste0("Ponderosa pine meadow, ",elevs[4],"m"),
           color="black", size = 5)+
  annotate("text", x=1.9, y=-1.2, label=paste0("Mixed conifer meadow, ",elevs[5],"m"),
           color="black", size = 5)
#ggsave("./Plots/Fig4.NMDS.jpg", height = 6, width = 12)


#Species composition - permANOVA
all_bray_rel = vegdist(allsp_rel, method = 'bray')
permanova_bray = adonis2(all_bray_rel ~ Treatment*type*elevfact, perm= 999, data = all_groups)
summary(permanova_bray)
permanova_bray


#Species composition - Dissimilarity
all_bray_rel = vegdist(allsp_rel, method = 'bray')

all_bray_long = all_bray_rel%>%
  dist2list() %>% 
  separate(col, "_", into= c("site1", "treatment1", "plot1", "type1")) %>%
  separate(row, "_", into=(c("site2", "treatment2", "plot2", "type2"))) %>%
  mutate(sitetreatplot1 = paste(site1, treatment1, plot1, sep="_")) %>%
  mutate(sitetreatplot2 = paste(site2, treatment2, plot2, sep="_")) 

bray_long_type = all_bray_long %>%
  filter(sitetreatplot1 == sitetreatplot2) %>%  #this filters to dissimilarities within plots
  filter(type1 != type2) %>% #this filters to seedbank vs above ground
  filter(!duplicated(sitetreatplot1)) %>%
  rename(bray = value) %>%
  mutate(Site = as.factor(site1), Treatment = as.factor(treatment1), Plot = as.numeric(plot1), 
         type1 = as.factor(type1), type2= as.factor(type2)) %>%
  select(-sitetreatplot1, -sitetreatplot2, -site2, -treatment2, -plot2, -site1, -treatment1, -plot1) %>%
  select(Site, Treatment, Plot, type1, type2, bray) 

bray_long_type = bray_long_type %>% 
  mutate(site_byelev = factor(Site, levels = c("Desert scrub", "Desert grassland", "Juniper savanna", "Ponderosa pine meadow", "Mixed conifer meadow"))) %>%
  mutate(trt_order = factor(Treatment, levels = c("Water Exclusion", "Control", "Water Addition"))) %>%
  mutate(elevation = case_when(
    Site == "Desert scrub" ~ 1566,
    Site == "Desert grassland" ~ 1636,
    Site == "Juniper savanna" ~ 1930,
    Site == "Ponderosa pine meadow" ~ 2179,
    Site == "Mixed conifer meadow" ~ 2591
  )) %>%
  mutate(elevfact = as.factor(elevation))

dim(bray_long_type) #60 plots in study - 1 plot had no seeds growout, so 59 to compare 
summary(bray_long_type)
#write.csv(bray_long_type, file = "../Formatted_data/bray_relabun_abovevssb_withinplot.csv")   

lm_bray_full  = lm(bray ~ elevfact*Treatment, data = bray_long_type)
summary(lm_bray_full)                   
anova(lm_bray_full)

#posthoc comparison
test(emmeans(lm_bray_full, pairwise ~ Treatment|elevfact) )
#Control vs water exclusion:
# marg sig at 1636 (P=0.087), 1930m
#water addition:
#sig at 1566, 1930

#Plot dissimilarity
bray_means = bray_long_type %>%
  group_by(Site, site_byelev, Treatment, trt_order, elevation, elevfact) %>%
  summarize(meanbray = mean(bray), sdbray = sd(bray), samplesize = n()) %>%
  mutate(sebray = sdbray/sqrt(samplesize)) %>%
  #add labels from contrasts for graphing
  mutate(siglabel = "Nonsignificant") %>%
  mutate(siglabel = case_when(
    elevation == 1566 & Treatment == "Water Addition" ~ "*", 
    elevation == 1636  & Treatment == "Water Exclusion" ~ "+",
    elevation == 1930  & Treatment == "Water Exclusion" ~ "*",
    elevation == 1930  & Treatment == "Water Addition" ~ "*"))

bray_plot_elevation = ggplot(bray_means, aes(x = elevation, y= meanbray, color = Treatment)) + 
  geom_point(position=position_dodge(10), size = 6) + #geom_line(position=position_dodge(0.3), aes(linetype = type)) +
  geom_errorbar(aes(ymin = meanbray - sebray, ymax =  meanbray + sebray), width = 0.2,
                position=position_dodge(10)) + theme_bw() + 
  scale_color_manual(values = c("darkolivegreen4",  "dodgerblue4", "firebrick4"  )) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = "Elevation (m, asl)", y = "Bray Curtis Dissimilarity", color = "Treatment", linetype = "Community type", shape = "Community type") +
  geom_text(aes(label = siglabel), 
            position = position_dodge(), hjust = -1, size =12, show.legend = F) + 
  theme(legend.justification = c(0.95, .05),legend.position = c(0.95,.05), 
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 20))

#### Fig. 5: Bray-Curtis distance #####
bray_plot_elevation
#ggsave("../Plots/fig5_braydissim.jpg", height = 8, width = 12)
