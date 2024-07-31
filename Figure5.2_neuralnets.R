##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
rm(list = ls())
library(tidyverse)
library(parallel)
library(foreach)
library(neuralnet)

mainDir<-getwd()
subdir<-paste0(mainDir,"/Fig5_machinelearning/neural")
if(!dir.exists(subdir)){dir.create(subdir)}


setwd(file.path(mainDir))
source("./functions/train_test_dataset_list.R")
load(file = "./adjustment/MachineLearn_data_clean.Rda")
setwd(file.path(subdir))
nmodeler<-function(x){
  x<-as.data.frame(x)
  x<-x[,4:10]
  model<-neuralnet(Variant~., data = x, hidden = c(5), rep = 1, linear.output = FALSE, stepmax = 1e9, threshold = 0.1)
  return(model)
}

#split dataset to vaccinated and unvaccinated---------

data_clean %>%
  filter(data_clean$Vaccination == "V" )%>%
  as.data.frame()->Vaccinated

data_clean %>%
  filter(data_clean$Vaccination == "U" )%>%
  as.data.frame()->Unvaccinated

vax_data_by_visit<-as.list(1)
visits<-unique(data_clean$Visit)
for(i in 1:length(visits)){
  Vaccinated %>%
    filter(Vaccinated$Visit == visits[i])%>%
    as.data.frame() -> vax_data_by_visit[[i]]
}
unvax_data_by_visit<-as.list(1)
visits<-unique(data_clean$Visit)
for(i in 1:length(visits)){
  Unvaccinated %>%
    filter(Unvaccinated$Visit == visits[i])%>%
    as.data.frame() -> unvax_data_by_visit[[i]]
}
#multicore------------
clust<-makeCluster(4)
clusterCall(clust,function() library(dplyr))
clusterExport(clust,c('vax_data_by_visit','train_and_test_dataset_list'))
vax<-parLapply(clust,vax_data_by_visit,train_and_test_dataset_list)
stopCluster(clust)

#vaccinated visits------------
accuracies_vax<-as.list(1)
allmodels_vax<-as.list(1)
for(j in 1:4){
  #j=2
  
  visit<-vax[[j]]
  train<-visit[[1]]
  test<-visit[[2]]
  clust<-makeCluster(4)
  clusterCall(clust,function() library(neuralnet))
  clusterExport(clust,c('test','train','nmodeler'))
  models<-parLapply(clust,train,nmodeler)
  stopCluster(clust)
  
  clust<-makeCluster(4)
  clusterCall(clust,function() library(dplyr))
  clusterCall(clust,function() library(neuralnet))
  clusterCall(clust,function() library(foreach))
  clusterExport(clust,c('test','models'))
  accuracy_l<-foreach(i=1:100) %dopar% {
    dat<-test[[i]]
    dat<-as.data.frame(dat)
    dat<-dat[,4:10]
    model<-models[[i]]
    names<-colnames(dat)
    position<-which(names == "Variant")
    pred <- predict(model, newdata = dat[,-position])
    labels <- c("A", "D", "O")
    prediction_label <- data.frame(max.col(pred)) %>%     
      mutate(pred=labels[max.col.pred.]) %>%
      select(2) %>%
      unlist()
    return(mean(prediction_label == dat$Variant))}
  stopCluster(clust)
  accuracies_vax[[j]]<-accuracy_l
  allmodels_vax[[j]]<-models
}
  
save(accuracies_vax, file = "./accuracies_vax.rds")
save(vax, file= "./datasets_vax.rds")
save(allmodels_vax, file ="./allmodels_vax.rds")



#multicore unvax
clust<-makeCluster(4)
clusterCall(clust,function() library(dplyr))
clusterExport(clust,c('unvax_data_by_visit','train_and_test_dataset_list'))
unvax<-parLapply(clust,unvax_data_by_visit,train_and_test_dataset_list)
stopCluster(clust)

#unvaccinated visits
accuracies_unvax<-as.list(1)
allmodels_unvax<-as.list(1)
for(j in 1:4){
  visit<-unvax[[j]]
  train<-visit[[1]]
  test<-visit[[2]]
  clust<-makeCluster(10)
  clusterCall(clust,function() library(neuralnet))
  clusterExport(clust,c('test','train','nmodeler'))
  models<-parLapply(clust,train,nmodeler)
  stopCluster(clust)
  
  clust<-makeCluster(5)
  clusterCall(clust,function() library(dplyr))
  clusterCall(clust,function() library(neuralnet))
  clusterCall(clust,function() library(foreach))
  clusterExport(clust,c('test','models'))
  accuracy_l<-foreach(i=1:10) %dopar% {
    dat<-test[[i]]
    dat<-as.data.frame(dat)
    dat<-dat[,4:10]
    model<-models[[i]]
    names<-colnames(dat)
    position<-which(names == "Variant")
    pred <- predict(model, newdata = dat[,-position])
    labels <- c("A", "D", "O")
    prediction_label <- data.frame(max.col(pred)) %>%     
      mutate(pred=labels[max.col.pred.]) %>%
      select(2) %>%
      unlist()
    return(mean(prediction_label == dat$Variant))}
  stopCluster(clust)
  accuracies_unvax[[j]]<-accuracy_l
  allmodels_unvax[[j]]<-models
}
save(accuracies_unvax, file = "./accuracies_unvax.rds")
save(unvax, file= "./datasets_unvax.rds")
save(allmodels_unvax, file ="./allmodels_unvax.rds")




#extract all the model accuracies
allaccuracies_vaccinated <- as.data.frame(matrix(NA,nrow = 100, ncol = 4))
for(k in 1:4){
  current_models <- allmodels_vax[[k]]
  visitdata <- vax[[k]]
  testdata <- visitdata[[2]]
  clust<-makeCluster(5)
  clusterCall(clust,function() library(dplyr))
  clusterCall(clust,function() library(neuralnet))
  clusterCall(clust,function() library(foreach))
  clusterExport(clust,c('test','models'))
  accuracy_list<-foreach(i=1:100) %dopar% {
    dat<-testdata[[i]]
    dat<-as.data.frame(dat)
    dat<-dat[,4:10]
    model<-current_models[[i]]
    names<-colnames(dat)
    position<-which(names == "Variant")
    pred <- predict(model, newdata = dat[,-position])
    labels <- c("A", "D", "O")
    prediction_label <- data.frame(max.col(pred)) %>%     
      mutate(pred=labels[max.col.pred.]) %>%
      select(2) %>%
      unlist()
    return(mean(prediction_label == dat$Variant))}
  stopCluster(clust)
  for(i in 1:100){
    allaccuracies_vaccinated[i,k] <- accuracy_list[[i]]
  }
}

#now for the unvaccinated
allaccuracies_unvaccinated <- as.data.frame(matrix(NA,nrow = 100, ncol = 4))
for(k in 1:4){
  current_models <- allmodels_unvax[[k]]
  visitdata <- unvax[[k]]
  testdata <- visitdata[[2]]
  clust<-makeCluster(5)
  clusterCall(clust,function() library(dplyr))
  clusterCall(clust,function() library(neuralnet))
  clusterCall(clust,function() library(foreach))
  clusterExport(clust,c('test','models'))
  accuracy_list<-foreach(i=1:100) %dopar% {
    dat<-testdata[[i]]
    dat<-as.data.frame(dat)
    dat<-dat[,4:10]
    model<-current_models[[i]]
    names<-colnames(dat)
    position<-which(names == "Variant")
    pred <- predict(model, newdata = dat[,-position])
    labels <- c("A", "D", "O")
    prediction_label <- data.frame(max.col(pred)) %>%     
      mutate(pred=labels[max.col.pred.]) %>%
      select(2) %>%
      unlist()
    return(mean(prediction_label == dat$Variant))}
  stopCluster(clust)
  for(i in 1:100){
    allaccuracies_unvaccinated[i,k] <- accuracy_list[[i]]
  }
}

save(allaccuracies_unvaccinated,file = "./unvaccinated_allaccuracies.Rda")
save(allaccuracies_vaccinated,file = "./vaccinated_allaccuracies.Rda")

writexl::write_xlsx(allaccuracies_unvaccinated,path = "./allaccuracies_unvaccinated.xlsx")
writexl::write_xlsx(allaccuracies_vaccinated,path = "./allaccuracies_vaccinated.xlsx")
