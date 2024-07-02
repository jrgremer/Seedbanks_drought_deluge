#### load libraries ####
library(dplyr)
library(tidyverse)
library(ggplot2)


#### load data ####
seeddat = read.csv("Raw_data/Seedbankgrowout_rawdata.csv")
summary(seeddat)
str(seeddat)
head(seeddat)


#NPP data from Seth, updated through 2021 (no 2022 yet)
plantdat = read.csv("Raw_data/Abovegrounddata_2022.csv")
summary(plantdat)
str(plantdat)
head(plantdat)

#metadata for plots and treatments
plotstreats = read.csv("Raw_data/Treatments_GRAMPS.csv") #treatment table is in repo, this is a relational pathway
head(plotstreats)

plotstreats = plotstreats %>%
  select(-quad) %>%
  distinct() #this removes extra rows at BP where 2 quadrats were censused for demography, but not for this study
dim(plotstreats)

#### format seedbank data:  cleaning, adding plots and treatments, wide and long formats ####

#make sure Tray = plot number sampled
seeddat = seeddat %>%
  rename(Plot = Tray)



#Remember to create a species list