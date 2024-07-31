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

source(paste0(mainDir,"/functions/test_train_dataset_generator.R"))
source(paste0(mainDir,"/functions/bruteforest.R"))

load(file = paste0(mainDir,"/adjustment/MachineLearn_data_clean.Rda"))


subdir <- paste0(mainDir,"/Fig5_machinelearning/forest/")

if(!dir.exists(subdir)){dir.create(subdir)}
setwd(subdir)
dir.create("./global",showWarnings = FALSE)
dir.create("./vaccinated",showWarnings = FALSE)
dir.create("./unvaccinated",showWarnings = FALSE)
#for debugging 
#setwd(mainDir)

#set the number of iterations
iterations = 99

#split data to vaccinated and unvaccinated-------
data_clean %>%
  filter(data_clean$Vaccination == "V" )%>%
  as.data.frame()->Vaccinated

data_clean %>%
  filter(data_clean$Vaccination == "U" )%>%
  as.data.frame()->Unvaccinated



#Calculation for Vaccinated--------------

#split data by visits

vax_data_by_visit<-as.list(1)
visits<-unique(data_clean$Visit)
for(i in 1:length(visits)){
  Vaccinated %>%
    filter(Vaccinated$Visit == visits[i])%>%
    as.data.frame() -> vax_data_by_visit[[i]]
}

#vaccinated models---------
best_test_visits<-as.list(1)
best_train_visits<-as.list(1)
best_accuracy_visits<-as.list(1)
best_prediction_visits<-as.list(1)
best_table_visits<-as.list(1)
best_forest_visits<-as.list(1)
vaccinated_allaccuracies <- as.data.frame(matrix(NA,nrow = iterations+1,ncol = length(visits)))
vaccinated_allforests<-as.list(1)
vaccinated_visitforests <- as.list(1)
vaccinated_alltest <- as.list(1)
vaccinated_alltrain <- as.list(1)
vaccinated_visittest <- as.list(1)
vaccinated_visittrain <- as.list(1)
vaccinated_alloob <- as.data.frame(matrix(NA,nrow = iterations+1,ncol = length(visits)))
for(v in 1:4){
  dat=vax_data_by_visit[[v]]
  seed_1 = 123
  print(paste("========================Visit", v))
  #generate test and train datasets, calculate initial forest
  train_and_test_dataset(dat = dat, seed_var = seed_1)
  train<-train_dataset[,4:10]
  test<-test_dataset[,4:10]
  bruteforest(train,test)
  best_test<-test_dataset
  best_train<-train_dataset
  best_accuracy<-accuracy
  best_prediction<-predictio
  best_table<-tabl
  best_forest<-forest
  vaccinated_allaccuracies[1,v]<-accuracy
  vaccinated_alloob[1,v]<-forest$err.rate[500,1]
  vaccinated_visitforests[[1]]<-forest
  vaccinated_visittest[[1]] <- test
  vaccinated_visittrain[[1]] <- train
  for(k in 1:iterations){
    seed_2 = seed_1 + k
    print(paste("iteration:", k))
    train_and_test_dataset(dat = dat, seed_var = seed_2)
    train2<-train_dataset[,4:10]
    test2<-test_dataset[,4:10]
    bruteforest(train2,test2)
    if(best_accuracy<accuracy){
      best_test<-test_dataset
      best_train<-train_dataset
      best_accuracy<-accuracy
      best_prediction<-predictio
      best_table<-tabl
      best_forest<-forest
    }
    vaccinated_allaccuracies[k+1,v]<-accuracy
    vaccinated_alloob[k+1,v]<-forest$err.rate[500,1]
    vaccinated_visitforests[[k+1]]<-forest
    vaccinated_visittest[[k+1]] <- test2
    vaccinated_visittrain[[k+1]] <- train2
  }
  best_test_visits[[v]]<-best_test
  best_train_visits[[v]]<-best_train
  best_accuracy_visits[[v]]<-best_accuracy
  best_prediction_visits[[v]]<-best_prediction
  best_table_visits[[v]]<-best_table
  best_forest_visits[[v]]<-best_forest
  vaccinated_allforests[[v]]<-vaccinated_visitforests
  vaccinated_alltest[[v]] <- vaccinated_visittest
  vaccinated_alltrain[[v]] <- vaccinated_visittrain
}

save(best_test_visits,file = "./vaccinated/best_test_visits.Rda")
save(best_train_visits,file = "./vaccinated/best_train_visits.Rda")
save(best_accuracy_visits,file = "./vaccinated/best_accuracy_visits.Rda")
save(best_prediction_visits,file = "./vaccinated/best_prediction_visits.Rda")
save(best_table_visits,file = "./vaccinated/best_table_visits.Rda")
save(best_forest_visits,file = "./vaccinated/best_forest_visits.Rda")
save(vaccinated_allaccuracies, file = "./vaccinated/vaccinated_allaccuracies.Rda")
save(vaccinated_allforests, file = "./vaccinated/vaccinated_allforests.Rda")
save(vaccinated_alltest, file = "./vaccinated/vaccinated_alltest.Rda")
save(vaccinated_alltrain, file = "./vaccinated/vaccinated_alltrain.Rda")
save(vaccinated_alloob, file = "./vaccinated/vaccinated_alloob.Rda")

write_xlsx(vaccinated_alloob, path = "./vaccinated/alloob.xlsx")
write_xlsx(vaccinated_allaccuracies, path = "./vaccinated/allaccuracies.xlsx")
#unvaccinated-----------------
unvax_data_by_visit<-as.list(1)
visits<-unique(data_clean$Visit)
for(i in 1:length(visits)){
  Unvaccinated %>%
    filter(Unvaccinated$Visit == visits[i])%>%
    as.data.frame() -> unvax_data_by_visit[[i]]
}

best_test_visits_unvax<-as.list(1)
best_train_visits_unvax<-as.list(1)
best_accuracy_visits_unvax<-as.list(1)
best_prediction_visits_unvax<-as.list(1)
best_table_visits_unvax<-as.list(1)
best_forest_visits_unvax<-as.list(1)
unvaccinated_allaccuracies <- as.data.frame(matrix(NA,nrow = iterations+1,ncol = length(visits)))
unvaccinated_allforests<-as.list(1)
unvaccinated_visitforests <- as.list(1)
unvaccinated_alltest <- as.list(1)
unvaccinated_alltrain <- as.list(1)
unvaccinated_visittest <- as.list(1)
unvaccinated_visittrain <- as.list(1)
unvaccinated_alloob <- as.data.frame(matrix(NA,nrow = iterations+1,ncol = length(visits)))

for(v in 1:4){
  seed_1 = 123
  dat=unvax_data_by_visit[[v]]    
  print("========================Visit")
  #generate test and train datasets, calculate initial forest
  train_and_test_dataset(dat = dat, seed_var = seed_1)
  train<-train_dataset[,4:10]
  test<-test_dataset[,4:10]
  bruteforest(train,test)
  best_test<-test_dataset
  best_train<-train_dataset
  best_accuracy<-accuracy
  best_prediction<-predictio
  best_table<-tabl
  best_forest<-forest
  unvaccinated_allaccuracies[1,v]<-accuracy
  unvaccinated_visittrain[[1]]<-train_dataset
  unvaccinated_visittest[[1]]<-test_dataset
  unvaccinated_visitforests[[1]]<-forest
  unvaccinated_alloob[1,v]<-forest$err.rate[500,1]
  
  for(k in 1:iterations){
    seed_2 = seed_1 + k
    print("iteration:")
    train_and_test_dataset(dat = dat, seed_var = seed_2)
    train2<-train_dataset[,4:10]
    test2<-test_dataset[,4:10]
    bruteforest(train2,test2)
    if(best_accuracy<accuracy){
      best_test<-test_dataset
      best_train<-train_dataset
      best_accuracy<-accuracy
      best_prediction<-predictio
      best_table<-tabl
      best_forest<-forest
    }
    unvaccinated_allaccuracies[k+1,v]<-accuracy
    unvaccinated_visittrain[[1+k]]<-train_dataset
    unvaccinated_visittest[[1+k]]<-test_dataset
    unvaccinated_visitforests[[1+k]]<-forest
    unvaccinated_alloob[k+1,v]<-forest$err.rate[500,1]
  }
  best_test_visits_unvax[[v]]<-best_test
  best_train_visits_unvax[[v]]<-best_train
  best_accuracy_visits_unvax[[v]]<-best_accuracy
  best_prediction_visits_unvax[[v]]<-best_prediction
  best_table_visits_unvax[[v]]<-best_table
  best_forest_visits_unvax[[v]]<-best_forest
  unvaccinated_allforests[[v]]<-vaccinated_visitforests
  unvaccinated_alltest[[v]] <- vaccinated_visittest
  unvaccinated_alltrain[[v]] <- vaccinated_visittrain
}
save(best_test_visits_unvax,file = "./unvaccinated/best_test_visits_unvax.Rda")
save(best_train_visits_unvax,file = "./unvaccinated/best_train_visits_unvax.Rda")
save(best_accuracy_visits_unvax,file = "./unvaccinated/best_accuracy_visits_unvax.Rda")
save(best_prediction_visits_unvax,file = "./unvaccinated/best_prediction_visits_unvax.Rda")
save(best_table_visits_unvax,file = "./unvaccinated/best_table_visits_unvax.Rda")
save(best_forest_visits_unvax,file = "./unvaccinated/best_forest_visits_unvax.Rda")
save(unvaccinated_allforests, file = "./unvaccinated/unvaccinated_allforests.Rda")
save(unvaccinated_alltrain, file = "./unvaccinated/unvaccinated_alltrain.Rda")
save(unvaccinated_alltest, file = "./unvaccinated/unvaccinated_alltest.Rda")
save(unvaccinated_allaccuracies, file = "./unvaccinated/unvaccinated_allaccuracies.Rda")
save(unvaccinated_alloob, file = "./unvaccinated/unvaccinated_alloob.Rda")

write_xlsx(unvaccinated_alloob, path = "./unvaccinated/alloob.xlsx")
write_xlsx(unvaccinated_allaccuracies, path = "./unvaccinated/allaccuracies.xlsx")
#load("./unvaccinated/best_test_visits_unvax.Rda")
#load("./unvaccinated/best_train_visits_unvax.Rda")
#load("./unvaccinated/best_accuracy_visits_unvax.Rda")
#load("./unvaccinated/best_prediction_visits_unvax.Rda")
#load("./unvaccinated/best_table_visits_unvax.Rda")
#load("./unvaccinated/best_forest_visits_unvax.Rda")

#final analysis-------------

variantaccuracyvax<-as.data.frame(matrix(NA,nrow=4,ncol = 3))
colnames(variantaccuracyvax)<-(c("A","D","O"))
for(i in 1:4){
  dat<-as.data.frame(best_table_visits[[i]])
  variantaccuracyvax[i,1]<-(dat[1,3]/sum(dat[1:3,3]))
  variantaccuracyvax[i,2]<-(dat[5,3]/sum(dat[4:6,3]))
  variantaccuracyvax[i,3]<-(dat[9,3]/sum(dat[7:9,3]))
}

variantaccuracyunvax<-as.data.frame(matrix(NA,nrow=4,ncol = 3))
colnames(variantaccuracyunvax)<-(c("A","D","O"))
for(i in 1:4){
  dat<-as.data.frame(best_table_visits_unvax[[i]])
  variantaccuracyunvax[i,1]<-(dat[1,3]/sum(dat[1:3,3]))
  variantaccuracyunvax[i,2]<-(dat[5,3]/sum(dat[4:6,3]))
  variantaccuracyunvax[i,3]<-(dat[9,3]/sum(dat[7:9,3]))
}
save(variantaccuracyvax,file = "./variantaccuracyvax.Rda")
save(variantaccuracyunvax,file ="./variantaccuracyunvax.Rda")
best_forest_visits_unvax[[4]]

save(vaccinated_allaccuracies, file = "./vaccinated_allaccuracies.Rda")
save(unvaccinated_allaccuracies, file = "./unvaccinated_allaccuracies.Rda")

write_xlsx(vaccinated_allaccuracies, path = "./vaccinated_allaccuracies.xlsx")
write_xlsx(unvaccinated_allaccuracies, path = "./unvaccinated_allaccuracies.xlsx")
