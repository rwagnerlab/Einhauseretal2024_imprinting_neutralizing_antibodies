##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
rm(list=ls())
library(dplyr)
library(writexl)
library(tidyverse)
#load data and define necessary variables
load("./adjustment/workspace_adjusted_long_simple.Rda")
workspace_long<- export_simple
rm(export_simple)
visits <- unique(as.character(workspace_long$visit))
steps <- c(0:2561)
`%notin%` <- Negate(`%in%`)

#load function
source("./functions/mag-breadth-function.R")
#calculate the magnitude-breadth curve

# split the data into the visits
visit_data <- vector("list", length(visits))
for( i in 1:length(visits)){
  visit_data[[i]]<-subset(workspace_long,workspace_long$visit == visits[i])
}


# Apply the function to each visit
visit_mb <- vector("list", length(visits))

for( i in 1:length(visits)){
  visit_mb[[i]] <- magbreadth(df = visit_data[[i]])
}

#now the same for mb without bq

#subset the workspace to not contain any bq values
workspace_nobq <- subset(workspace_long, workspace_long$virus %notin% c("BQ","JN","XBB"))

visit_data_nobq <- vector("list", length(visits))
for( i in 1:length(visits)){
  visit_data_nobq[[i]]<-subset(workspace_nobq,workspace_nobq$visit == visits[i])
}

# Pre-allocate list to store results for mb without bq
# Apply the function to each visit
visit_mb_nobq <- vector("list", length(visits))

for( i in 1:length(visits)){
  visit_mb_nobq[[i]] <- magbreadth(df = visit_data_nobq[[i]])
}

#export the lists with the mb data
path = "./fig2_mb/"
dir.create(path)

for(i in 1:length(visits)){
  mb_export <- cbind(steps,visit_mb[[i]])
  mb_export_nobq <- cbind(steps,visit_mb_nobq[[i]])
  write_xlsx(mb_export, path = paste(path,"/mb_",visits[i],".xlsx", sep = ""))
  write_xlsx(mb_export_nobq, path = paste(path,"/mb_",visits[i],"_nobq.xlsx", sep = ""))
}

######aucs##########
#for statistical comparison of mb aucs, auc of curve = mean auc of single curves , auc of single curve is equal to mean as shown in huang et al
workspace_wide <- workspace_long %>%
  pivot_wider( names_from = virus, values_from = adjusted_neutralization)

workspace_wide$auc_nobq <- rowMeans(workspace_wide[, c("D614G", "Alpha", "Delta", "BA1", "BA2", "BA5")], na.rm = TRUE)
workspace_wide$auc_bq <- rowMeans(workspace_wide[, c("D614G", "Alpha", "Delta", "BA1", "BA2", "BA5", "BQ", "JN", "XBB")], na.rm = TRUE)

#subset the visits
visit_aucs <- vector("list", length(visits))
visit_aucs_nobq <- vector("list", length(visits))
for (i in 1: length(visits)){
  visit_data_auc<-subset(workspace_wide, workspace_wide$visit == visits[i])
  visit_data_auc<- subset(visit_data_auc, select = c(ID_study ,group, auc_bq))
  visit_data_auc <- visit_data_auc %>%
    pivot_wider(names_from = group, values_from = auc_bq)
  visit_aucs[[i]]<-visit_data_auc
  
  visit_data_auc_nobq<-subset(workspace_wide, workspace_wide$visit == visits[i])
  visit_data_auc_nobq<- subset(visit_data_auc_nobq, select = c(ID_study ,group, auc_nobq))
  visit_data_auc_nobq <- visit_data_auc_nobq %>%
    pivot_wider(names_from = group, values_from = auc_nobq)
  visit_aucs_nobq[[i]]<-visit_data_auc_nobq
}

#export the aucs
for(i in 1:length(visits)){
  auc_export<-visit_aucs_nobq[[i]]
  # Specify the desired order of columns
  desired_order <- c("ID_study", "FA", "UA", "FD", "UD", "FO", "UO")
  # Rearrange the columns
  auc_export_reordered <- auc_export[, desired_order]
  #write excel
  write_xlsx(auc_export_reordered, path = paste(path,"/auc_mb_",visits[i],"_nobq.xlsx", sep = ""))
}

