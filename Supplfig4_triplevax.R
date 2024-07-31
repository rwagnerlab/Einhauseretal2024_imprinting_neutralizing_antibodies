##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
### script to generate the arranged data for the triplevax plots (suppl Fig 3)

library(writexl)
library(tidyverse)
library(readxl)
rm(list = ls())

load("./adjustment/workspace_adjusted_long_simple.Rda")
#read in the 3x vaccinated IDs
vaxIDs<- read_xlsx("./3xvax.xlsx")
data <- subset(export_simple, export_simple$group == "FO")

# Add column to check if ID is in vaxIDs aka triple vaccinated
data <- data %>%
  mutate(InVaxIDs = if_else(ID_study %in% vaxIDs$ID, TRUE, FALSE))

#subset for only the variants of interest
VOI <- unique(data$virus)

virus_dataframes <- list()

# Iterate over each virus type
for (virus in VOI) {
  # Filter the dataframe for the current virus
  virus_df <- data %>% filter(virus == !!virus)
  
  # Create a new dataframe with columns for each visit combined with 3xvax
  reshaped_df <- virus_df %>%
    unite("visit_vax", InVaxIDs, visit, sep = "_") %>%
    spread(key = visit_vax, value = adjusted_neutralization)
  
  # Add the reshaped dataframe to the list
  virus_dataframes[[virus]] <- reshaped_df
}

# create output directory
if(!dir.exists("./Supplfig3_triplevax")){dir.create("./Supplfig3_triplevax")}
for(i in VOI){
  write_xlsx(virus_dataframes[[i]], path = paste0("./Supplfig3_triplevax/", i, ".xlsx"))
}
