##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses
#generates antigenic maps
library(Racmacs)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(patchwork)
library(png)
rm(list = ls())
options(RacOptimizer.num_cores = parallel::detectCores())

#load in the data
load("./adjustment/workspace_adjusted_wide_simple.Rda")
exportwide<-exportwide_simple
rm(exportwide_simple)
#adjust colors  and shape to the like############
# generate sera colors and shapes based on infection variant and vaccine status
exportwide <- exportwide %>%
  mutate(color = case_when(
    group == "FO" ~ '#F56E07',
    group == "FD" ~ '#2980B9',
    group == "FA" ~ '#b000b5',
    group == "UO" ~ '#faa911',
    group == "UD" ~ '#4792ba',
    group == "UA" ~ '#9B59B6',
    TRUE ~ NA_character_
  ))
exportwide <- exportwide %>%
  mutate(shape = case_when(
    vaccine == 1 ~ "22",
    vaccine == 0 ~ "21",
    TRUE ~ NA_character_
  ))

#set the virus colors
virus.col <- c('#C0392B', '#9B59B6', '#2980B9',
               '#E67E22','#F1C40F', '#1ABC9C', '#34495e','#5f495e','#7f495e')
# replace NAs with stars *
exportwide <- exportwide %>%
  mutate_at(vars(D614G:JN), ~replace(., is.na(.), "*"))


#set opacity for vaccinated and unvaccinated submaps
vaccinated <- unique(exportwide$vaccine)
exportwide <- exportwide %>%
  mutate(Vaxopacity = case_when(
    vaccine == 1 ~ "0.6",
    vaccine == 0 ~ "0",
    TRUE ~ NA_character_
  ))%>%
  mutate(Unvaxopacity = case_when(
    vaccine == 1 ~ "0",
    vaccine == 0 ~ "0.6",
    TRUE ~ NA_character_
  ))

#split data to visits ##########
visits <- as.character(unique(exportwide$visit))
visit_data<- as.list(NA)
for( i in 1:length(visits)){
  visit_data[[i]] <- subset(exportwide, visit == visits[i])
}



#now make the maps

maps <- as.list(NA)

for(i in 1:length(visits)){
  #create the titer table
  
  visit_dat <- visit_data[[i]]
  titertab <- as.data.frame(visit_dat[,6:14])
  
  row.names(titertab) <- visit_dat$ID_study
  sera.col <- visit_dat$color
  sera.shape <- visit_dat$shape
  titertab_transposed <- as.data.frame(t(titertab))
  #create the map 
  map <- acmap( titer_table = titertab_transposed)
  map <- optimizeMap( map = map, number_of_dimensions = 2, number_of_optimizations = 1000, minimum_column_basis = "none", options = list(ignore_disconnected = TRUE))
  
  #improve optics
  srGroups(map)<-as.character(visit_dat$group)
  agSize(map) <- 40
  srSize(map) <- 8
  names(virus.col) <- agNames(map)
  names(sera.col) <- srNames(map)
  names(sera.shape) <- srNames(map)
  #antigen colors
  agFill(map) <-virus.col
  agOutline(map)<-virus.col
  agOpacity(map) <- 0.6
  
  #sera colors and shape
  srOutline(map)<-"black"
  srFill(map)<-sera.col
  srShape(map)<-sera.shape
  # hide either vaccinated or unvaccinated from the map
  for( j in 1:(length(vaccinated)+1)){
    if(j == 1){
      sera.opac = 0.6
      k = i
    }
    else if(j == 2){
      sera.opac <- visit_dat$Vaxopacity
      names(sera.opac)<- srNames(map)
      k = i + length(visits)
    }else{
      sera.opac <- visit_dat$Unvaxopacity
      names(sera.opac)<- srNames(map)
      k = i + 2*length(visits)
    }
    
    srOpacity(map) <- as.numeric(sera.opac)
   
    #save the generated map to a list
    maps[[k]] <- map
  }
 
  #clean up workspace
  rm(map)
  rm(titertab)
  rm(titertab_transposed)
  rm(visit_dat)
}

#realign all maps to v2 for congruent optics
for(i in 1:12){
  if( i != 2){
    maps[[i]] <- realignMap(maps[[i]], maps[[2]], translation = TRUE, scaling = FALSE)
  }
}




#save maps
if(!dir.exists("./fig3_mapping/")){dir.create("./fig3_mapping/")}
for(i in 1:12){
save.acmap(maps[[i]], filename = paste0("./fig3_mapping/map",i,".ace"))
}

  


