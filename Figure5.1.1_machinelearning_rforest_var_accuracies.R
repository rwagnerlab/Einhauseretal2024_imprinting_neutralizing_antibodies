##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#calculates random forests
library(party)
library(tidyverse)
library(writexl)
library(readxl)
library(ggraph)
library(neuralnet)
library(cito)
library(randomForest)

#clean up, set up environment--------------
rm(list = ls())
mainDir<- getwd()
subdir <- paste0(mainDir,"/Fig5_machinelearning/forest/")

load(paste0(subdir,"vaccinated/vaccinated_allforests.Rda"))
load(paste0(subdir,"vaccinated/vaccinated_alltest.Rda"))

visits <-c("v1","v2","v4","v5")
variants <- c("A","D","O")
variant_accuracies <- vector("list",4)

for(i in 1:length(visits)){
  testsets <- vaccinated_alltest[[i]]
  models <- vaccinated_allforests[[i]]
  varaccuracies <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 3))
  colnames(varaccuracies) <- c("A", "D", "O")
  for(j in 1:100){
    test <- testsets[[j]]
    model<- models[[j]]
    #prediction with the model
    predictions <-predict(model,newdata = test)
    #prediction vs true value
    confusion_matrix <- table(predictions,test$Variant)
    rownames(confusion_matrix) <- c("A", "D", "O")
    colnames(confusion_matrix) <- c("A", "D", "O")
    
    # Extracting the diagonal elements (True Positives)
    TP <- diag(confusion_matrix)
    # Calculating the sum of rows (True Positives + False Negatives)
    row_sums <- rowSums(confusion_matrix)
    
    # Calculating subaccuracies
    subaccuracies <- TP / row_sums
    
    # Displaying the subaccuracies
    names(subaccuracies) <- rownames(confusion_matrix)
    varaccuracies[j,]<-subaccuracies
  }
  variant_accuracies[[i]]<-varaccuracies
}
save(variant_accuracies,file = paste0(subdir,"vaccinated/variant_accuracies.Rda"))

all_vaccinated_variant_accuracies<-variant_accuracies[[1]]
for (i in 2:length(visits)) { all_vaccinated_variant_accuracies<- cbind(all_vaccinated_variant_accuracies,variant_accuracies[[i]])}
write_xlsx(all_vaccinated_variant_accuracies,path = paste0(subdir,"vaccinated/variant_accuracies.xlsx"))

#same for unvaccinated models
load(paste0(subdir,"unvaccinated/unvaccinated_allforests.Rda"))
load(paste0(subdir,"unvaccinated/unvaccinated_alltest.Rda"))

unvax_variant_accuracies <- vector("list",4)

for(i in 1:length(visits)){
  uvx_testsets <- unvaccinated_alltest[[i]]
  uvx_models <- unvaccinated_allforests[[i]]
  uvx_varaccuracies <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 3))
  colnames(uvx_varaccuracies) <- c("A", "D", "O")
  for(j in 1:100){
    test <- uvx_testsets[[j]]
    model<- uvx_models[[j]]
    #prediction with the model
    predictions <-predict(model,newdata = test)
    #prediction vs true value
    confusion_matrix <- table(predictions,test$Variant)
    rownames(confusion_matrix) <- c("A", "D", "O")
    colnames(confusion_matrix) <- c("A", "D", "O")
    
    # Extracting the diagonal elements (True Positives)
    TP <- diag(confusion_matrix)
    # Calculating the sum of rows (True Positives + False Negatives)
    row_sums <- rowSums(confusion_matrix)
    
    # Calculating subaccuracies
    subaccuracies <- TP / row_sums
    
    # Displaying the subaccuracies
    names(subaccuracies) <- rownames(confusion_matrix)
    uvx_varaccuracies[j,]<-subaccuracies
  }
  unvax_variant_accuracies[[i]]<-varaccuracies
}
save(unvax_variant_accuracies,file = paste0(subdir,"unvaccinated/unvx_variant_accuracies.Rda"))
#export excel file
all_unvaccinated_variant_accuracies<-unvax_variant_accuracies[[1]]
for (i in 2:length(visits)) { all_unvaccinated_variant_accuracies<- cbind(all_unvaccinated_variant_accuracies,unvax_variant_accuracies[[i]])}
write_xlsx(all_unvaccinated_variant_accuracies,path = paste0(subdir,"unvaccinated/unvax_variant_accuracies.xlsx"))
