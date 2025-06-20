#### load libraries ####
library(dplyr)
library(tidyverse)
library(ggplot2)


#### load data ####
seeddat = read.csv("Raw_data/Seedbank_traytally.csv") #this file has corrected species names
summary(seeddat)
str(seeddat)
head(seeddat)


#NPP data from Seth, updated through 2021 
plantdat = read.csv("Raw_data/Abovegrounddata_2021.csv") #this file has corrected species names and PFTs, and only 2021 data (to match seedbank data)
summary(plantdat)
str(plantdat)
head(plantdat)

#NPP data has plant functional type (PFT) data, so create a species and PFT dataframe to populate for seedbank data

plantspPFTlist =  plantdat %>%
  select(Species, PFT) %>%
  distinct()

#metadata for plots and treatments
plotstreats = read.csv("Raw_data/Treatments_GRAMPS.csv") #treatment table is in repo, this is a relational pathway
head(plotstreats)

plotstreats = plotstreats %>%
  select(-quad) %>%
  distinct() #this removes extra rows at BP where 2 quadrats were censused for demography, but not for this study
dim(plotstreats)

#### format seedbank data ####
#join to get treatment and make treatment a factor
dim(seeddat)
seeddat = left_join(seeddat, plotstreats, by = c("Site", "Plot"))%>%
          mutate(Treatment = ifelse(Site == "Blank", "Blank", Treatment)) %>% #add "blank" to treatments
          mutate(Treatment = as.factor(Treatment), Site = as.factor(Site))

dim(seeddat)   
head(seeddat)
summary(seeddat)


#join with species and PFT list and make PFT a factor
seeddat = left_join(seeddat, plantspPFTlist) %>%
          mutate(PFT = as.factor(PFT))
          
dim(seeddat) 
head(seeddat)
summary(seeddat)  

unique(as.factor(seeddat$PFT)) #have some NA's for PFTs, which may mean they are not in NPP dataset, at least not the 2021 data

filter(seeddat, is.na(PFT) == T) 

#Adding PFTs 
seeddat = seeddat %>%
  mutate(PFT = factor(if_else(Species == "AMDI", "PF", as.character(PFT)))) %>% 
  mutate(PFT = factor(if_else(Species == "BASC", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "CHAL", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "CHSE", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "MEMO", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "BLTR", "PG4", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "OXLA", "PF", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "ERSP", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "MEMO", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "MEAL", "A", as.character(PFT)))) %>%
  mutate(PFT = as.factor(PFT))

#correct numbers per to emergence in blank, control trays deployed during the growout
#get sum of seedlings in control plots, these likely represent contamination in the soil (common weeds near Flagstaff)
seeddat[seeddat$Site == "Blank",]

blanktots = seeddat %>%
  filter(Site == "Blank") %>%
  group_by(Species, PFT) %>%
  summarize(blanktot = mean(Tally))


#sum by tray
seedtots_long = seeddat %>% 
  group_by(Site, Plot, Treatment, Species, PFT) %>%
  summarise(tottally = sum(Tally, na.rm=T)) %>% #for seedbank data, abundances = counts (will be cover for emergent community)
  ungroup() 

#join with blanktotals 
dim(seedtots_long)
seedtots_long = left_join(seedtots_long, blanktots)
dim(seedtots_long)
head(seedtots_long)
summary(seedtots_long)

#correct for contamination
seedtots_long = seedtots_long %>%
  mutate(blanktot = if_else(is.na(blanktot) == T, 0, blanktot)) %>%
  mutate(totabun = tottally - blanktot) %>%  #raw abundance is actually corrected total abundance
  mutate(totabun = if_else(totabun <0, 0, totabun)) #if that subtraction made it less than zero, correct to zero

summary(seedtots_long)

#Calculate relative abundances
seedtots_long = seedtots_long %>%
  mutate(Site = as.factor(Site), Plot = as.factor(Plot), Treatment = as.factor(Treatment), Species = as.factor(Species)) %>%
  #calculate relative abundance at the plot level
  group_by(Site, Plot, Treatment) %>%
  mutate(relabun = ifelse(totabun >0, totabun/sum(totabun),0)) %>%
  ungroup() %>%
  mutate(Year = 2021, type = as.factor("seedbank")) #add year and community type

summary(seedtots_long)
head(seedtots_long)
dim(seedtots_long)

write.csv(seedtots_long, file = "Formatted_data/seedbank_rawandrelabun_long.csv")

#create Wide dataframes
sb_sprelabun_wide = seedtots_long %>%
  select(-totabun, -PFT) %>% 
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = relabun)
summary(sb_sprelabun_wide)
str(sb_sprelabun_wide)

write.csv(sb_sprelabun_wide, file = "Formatted_data/seedbank_relabun_wide.csv")

sb_sptotabun_wide = seedtots_long %>%
  select(-relabun, -PFT) %>% 
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = totabun)
summary(sb_sptotabun_wide)
str(sb_sptotabun_wide)
write.csv(sb_sptotabun_wide, file = "Formatted_data/seedbank_totabun_wide.csv")

#### format NPP/above ground data ####

#remember here there are 2 quadrats/plot.   
plant_quad = plantdat %>%
  rename(treat = Treatment) %>%
  mutate(Treatment = recode(treat, "C" = "Control", "W" = "Water Addition", "D" = "Water Exclusion")) %>%
  mutate(type = "Aboveground") %>%
  #rename(raw_canopy_cov = Canopy_cover, raw_NPP = NPP_g_m2) %>%
  mutate(Site = as.factor(Site), Plot = as.factor(Plot), Treatment = as.factor(Treatment), Species = as.factor(Species)) %>%
  group_by(Year, Site, Plot, Treatment, Species, PFT, Quadrat) %>% #calculate cover at the quadrat level (per email from Seth 10/19/22)
  summarize(quadcover = sum(Canopy_cover*Count)) %>%
  ungroup()

plantdat_plot = plant_quad %>%
  #calculate raw canopy cover and then relative abundance based on cover at the species level by plot 
  group_by(Year, Site, Plot, Treatment, Species) %>%
  mutate(canopycov_plot = sum(quadcover, na.rm=T)/2) %>% #take mean of canopy cover across 2 NPP quadrats 
  ungroup() %>%
  select(-Quadrat, -quadcover) %>%
  distinct() #removes duplicate values because of quadrat


plantdat_long = plantdat_plot %>%
  group_by(Year, Site, Plot, Treatment) %>%
  mutate(relabun = canopycov_plot/sum(canopycov_plot)) %>%
  mutate(type = as.factor("aboveground")) %>%
  ungroup()

summary(plantdat_long)            

write.csv(plantdat_long, file = "Formatted_data/aboveground_NPP_rawandrelabun_long.csv")

#create Wide dataframes
ab_sprelabun_wide = plantdat_long %>%
  select(-canopycov_plot, -PFT) %>%
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = relabun)
summary(ab_sprelabun_wide)
str(ab_sprelabun_wide)

write.csv(ab_sprelabun_wide, file = "Formatted_data/aboveground_NPP_relabun_wide.csv")

#### merge aboveground and seedbank data ####

#use just collection year for emergent
plantdat_long = plantdat_long %>%
    rename(totabun = canopycov_plot) #rename raw abundance to match seedbank data

summary(plantdat_long)

length(unique(seedtots_long$Species)) #36 unique species in seedbank data
length(unique(plantdat_long$Species)) #72 unique species in above ground data

length(intersect(unique(seedtots_long$Species), unique(plantdat_long$Species))) #27 species in common in both emergent and seedbank data

# Merge dataframes: long format #
names(seedtots_long)
names(plantdat_long)


dim(seedtots_long)
dim(plantdat_long)

alldat_long = full_join(seedtots_long, plantdat_long) 

dim(alldat_long)
dim(seedtots_long) + dim(plantdat_long) #dimensions look good

summary(alldat_long)

write.csv(alldat_long, file = "Formatted_data/Seedbank_aboveground_merged_long.csv")

## Wide format dataframes ##
all_sprelabun_wide = alldat_long %>%
  select(-totabun, -PFT) %>%
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = relabun)
summary(all_sprelabun_wide)
str(all_sprelabun_wide)

write.csv(all_sprelabun_wide, file = "Formatted_data/Relabun_sbandab_merged_wide.csv")

all_sptotabun_wide = alldat_long %>%
  select(-relabun, -PFT) %>%
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = totabun)
summary(all_sptotabun_wide)
str(all_sptotabun_wide)

write.csv(all_sptotabun_wide, file = "Formatted_data/Totabun_sbandab_merged_wide.csv")

#### Adding native vs non-native (introduced) ####
summary(alldat_long)
nativeint = read.csv("Raw_data/species list_allNPPandSB_NI.csv") %>%
  select(Species, Native.Introduced) 

summary(nativeint)

dim(alldat_long)

#join dataframes to get native vs introduced  
alldat_long =  left_join(alldat_long, nativeint) %>%
    mutate(Species = as.factor(Species), Native.Introduced = as.factor(Native.Introduced)) 
dim(alldat_long)
summary(alldat_long)

alldat_long[is.na(alldat_long$Native.Introduced)==T,] #only species that were not identified in above ground data 

write.csv(alldat_long, file = "Formatted_data/Seedbank_aboveground_merged_long.csv")


#### Create and save species list #### 
splist = alldat_long %>%
  select(Species, PFT, Native.Introduced) %>%
  distinct()
dim(splist)
head(splist)

write.csv(splist, file = "Formatted_data/species list_allNPPandSB.csv")
