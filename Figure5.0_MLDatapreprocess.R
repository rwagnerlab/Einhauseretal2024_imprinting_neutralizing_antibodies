##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#preprocesses and cleans data for the ml models
library(tidyverse)

rm(list = ls())

load("./adjustment/workspace_adjusted_wide_simple.Rda")
data <- exportwide_simple
rm(exportwide_simple)

visits <- c("v1","v2","v4","v5")
sera_group <- c("FA","FD","FO","UA","UD","UO")
variants <- c("D614G","Alpha","Delta","BA1","BA2","BA5")

#clean up data
data %>%
  as.data.frame()%>%
  filter(data$group %in% sera_group)-> data_clean
data_clean %>%
  filter(data_clean$visit %in% visits)%>%
  select(-BQ,-JN,-XBB)  -> data_clean
  
#rename cols
colnames(data_clean)<- c("ID_study","Visit","Variant","Vaccination","Group",variants)

#change the vaccinated to v and u for vaccinated and unvaccinated
data_clean$Vaccination <- as.character(data_clean$Vaccination)
data_clean[data_clean$Vaccination == 1,"Vaccination"] <-"V"
data_clean[data_clean$Vaccination == 0,"Vaccination"] <-"U"
data_clean$Vaccination <- as.factor(data_clean$Vaccination)

#clear all NAs  
data_clean<-na.omit(data_clean)
#rearrange order of columns
new_order_of_columns<-c("ID_study","Visit","Group","Variant",variants,"Vaccination")
data_clean <- data_clean %>% select(all_of(new_order_of_columns))

save(data_clean,file = "./adjustment/MachineLearn_data_clean.Rda")
