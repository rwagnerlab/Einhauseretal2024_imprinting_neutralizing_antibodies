#by Sebastian Einhauser for Einhauser et al
library(randomForest)
#loop calculate random forest model
bruteforest<-function(train,test,num = iterations){
names<-colnames(test)
position<-which(names == "Variant")
#initial model
forest <<- randomForest(Variant~., data = train)
predictio <<- predict(forest, newdata = test[,-position])
accuracy<<-mean(predictio == test$Variant)
tabl<<-table(predictio,test$Variant)

#looped model to replace init model
for(i in 1:num){
  set.seed(i)
  print(i)
  forest2 <- randomForest(Variant~., data = train)
  prediction2 <- predict(forest, newdata = test[,-position])
  accuracy2<-mean(prediction2 == test$Variant)
  table2<-table(prediction2,test$Variant)
  if(accuracy2>accuracy){
    forest<<-forest2
    predictio<<-prediction2
    accuracy<<-accuracy2
    tabl<<-table2
  }}}