##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
rm(list=ls())

library(dplyr)
library(writexl)
library(readxl)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
#load function
source("./functions/mag-breadth-function.R")


#load data and define necessary variables
load("./adjustment/workspace_adjusted_long_simple.Rda")

#load the 3x vaccinated information
vaxIDs<- read_xlsx("./3xvax.xlsx")

export_simple <- export_simple%>%
  mutate(InVaxIDs = if_else(ID_study %in% vaxIDs$ID, TRUE, FALSE))

export_simple <- export_simple%>%
  mutate(grouptriple = if_else(InVaxIDs == TRUE, paste0(group,"3x"), paste0(group,"2x")))

#overwrite group with the triplevaccinated information

export_simple$group <- export_simple$grouptriple
export_simple<- export_simple%>%
  select(-InVaxIDs, -grouptriple)
save(export_simple, file = "./adjustment/triplevax.Rda")

workspace_long<- export_simple
rm(export_simple)
rm(vaxIDs)

visits <- unique(as.character(workspace_long$visit))
steps <- c(0:2561)


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

#export the lists with the mb data and fuse the original FO values to it
path <- "./supplfig5_triplevaxmb"
if(!dir.exists(path)){dir.create(path)}
path_mbmainfig <-"./fig2_mb/"
desired_ordering <- c("steps","FA","FD","FO","UA","UD","UO","O2xV","O3xV")
for(i in 1:length(visits)){
  mb_export <- cbind(steps,visit_mb[[i]])
  mb_export_nobq <- cbind(steps,visit_mb_nobq[[i]])
  # read the main figure mb curves and paste the FO curve here 
  mainfig <-read_xlsx(paste0(path_mbmainfig,"mb_",visits[i],".xlsx"))
  mainfig_nobq <-read_xlsx(paste0(path_mbmainfig,"mb_",visits[i],"_nobq.xlsx"))
  mb_export <- cbind(mb_export,mainfig$FO)
  mb_export_nobq <- cbind(mb_export_nobq, mainfig_nobq$FO)
  colnames(mb_export) <- c("steps","FA","FD","O2xV","O3xV","UA","UD","UO","FO")
  mb_export <- mb_export[,desired_ordering]
  colnames(mb_export_nobq) <- c("steps","FA","FD","O2xV","O3xV","UA","UD","UO","FO")
  mb_export_nobq <- mb_export_nobq[,desired_ordering]
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
aucs_fused<-as.data.frame(matrix(NA,nrow = length(visit_aucs_nobq[[1]]$FO3x)))
for(i in 1:length(visits)){
  auc_export<-visit_aucs_nobq[[i]]
  # Specify the desired order of columns
  desired_order <- c( "FO2x","FO3x")
  # Rearrange the columns
  auc_export_reordered <- auc_export[, desired_order]
  colnames(auc_export_reordered) <- c(paste0("O2xV-",visits[i]),paste0("O3xV-",visits[i]))
  aucs_fused <- cbind(aucs_fused,auc_export_reordered)
}
aucs_fused<- aucs_fused%>%
  select(-V1)
#write excel
write_xlsx(aucs_fused, path = paste(path,"/auc_mb_all_nobq.xlsx", sep = ""))
