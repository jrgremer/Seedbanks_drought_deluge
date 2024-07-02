#### load libraries ####
library(dplyr)
library(tidyverse)
library(ggplot2)


#### load data ####
seeddat = read.csv("Raw_data/Seedbank_traytally.csv") #this file has corrected species names
summary(seeddat)
str(seeddat)
head(seeddat)


#NPP data from Seth, updated through 2021 (no 2022 yet)
plantdat = read.csv("Raw_data/Abovegrounddata_2022.csv") #this file has corrected species names and PFTs
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
          mutate(Treatment = as.factor(Treatment))

dim(seeddat)   
head(seeddat)
summary(seeddat)

filter(seeddat, is.na(Treatment) == T) #only NAs are the blank trays

#join with species and PFT list and make PFT a factor
seeddat = left_join(seeddat, plantspPFTlist) %>%
          mutate(PFT = as.factor(PFT))
          
dim(seeddat) 
head(seeddat)
summary(seeddat)  

unique(as.factor(seeddat$PFT)) #have some NA's for PFTs, which may mean they are not in NPP dataset

filter(seeddat, is.na(PFT) == T) 

#Adding PFTs for those I can
seeddat = seeddat %>%
  mutate(PFT = factor(if_else(Species == "AMDI", "PF", as.character(PFT)))) %>% 
  mutate(PFT = factor(if_else(Species == "BASC", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "CHAL", "A", as.character(PFT)))) %>%
  mutate(PFT = factor(if_else(Species == "CHSE", "A", as.character(PFT)))) %>%
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
  summarise(tottally = sum(Tally, na.rm=T)) %>% #for seedbank data, raw abundances = counts (will be cover for emergent community)
  ungroup() 

#join with blanktotals and correct for contamination
dim(seedtots_long)
seedtots_long = left_join(seedtots_long, blanktots)
dim(seedtots_long)
head(seedtots_long)

#correct for contamination
seedtots_long = seedtots_long %>%
  mutate(blanktot = if_else(is.na(blanktot) == T, 0, blanktot)) %>%
  mutate(totabun = tottally - blanktot) %>%  #raw abundance is actually corrected total abundance
  mutate(totabun = if_else(totabun <0, 0, totabun)) #if that subtraction made it less than zero, correct to zero

summary(seedtots_long)

#Calculate relative abundances
seedtots_long = seedtots_long %>%
  mutate(Treatment = ifelse(Site == "Blank", "Blank", Treatment)) %>% #add "blank" to treatments
  mutate(Site = as.factor(Site), Plot = as.factor(Plot), Treatment = as.factor(Treatment), Species = as.factor(Species)) %>%
  #calculate relative abundance at the plot level
  group_by(Site, Plot, Treatment) %>%
  mutate(relabun = totabun/sum(totabun)) %>%
  ungroup() %>%
  mutate(Year = 2021, type = as.factor("seedbank"))

summary(seedtots_long)
head(seedtots_long)
dim(seedtots_long)

#write.csv(seedtots_long, file = "Formatted_data/seedbank_rawandrelabun_long.csv")

#create Wide dataframes
sb_sprelabun_wide = seedtots_long %>%
  select(-totabun, -PFT) %>% 
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = relabun)
summary(sb_sprelabun_wide)
str(sb_sprelabun_wide)

#write.csv(sb_sprelabun_wide, file = "Formatted_data/seedbank_relabun_wide.csv")

sb_sptotabun_wide = seedtots_long %>%
  select(-relabun, -PFT) %>% 
  pivot_wider(
    id_cols = c(Site, Plot, Treatment, Year, type), 
    names_from = Species, 
    values_from = totabun)
summary(sb_sptotabun_wide)
str(sb_sptotabun_wide)
#write.csv(sb_sptotabun_wide, file = "Formatted_data/seedbank_totabun_wide.csv")



#Remember to create a species list