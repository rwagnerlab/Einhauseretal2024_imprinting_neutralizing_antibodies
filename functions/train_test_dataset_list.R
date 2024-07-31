##by Sebastian Einhauser for Einhauser et al. 
train_and_test_dataset_list<-function(dat){
  seed_var = 123
  visit_sets<-as.list(1)
  #split into learning and testing dataset 10x
  train_list<-as.list(1)
  test_list<-as.list(1)
  for(x in 1:100){
    #create samples for variant
    seed_var2 = seed_var + x 
    set.seed(seed_var2)
    
    variants<-sort(as.character(unique(dat$Variant)))
    variant_freq<-as.data.frame(table(dat$Variant))
    n_samples_variant<-round((variant_freq$Freq)/2)
    drawn_variant<-as.list(1)
    dat<-as.data.frame(dat)
    #sample the test and training dataset by randomly splitting in half
    for(i in 1:length(variants)){
      dat%>%
        filter(Variant == variants[i])-> interim
      indices<-sample(nrow(interim), size = n_samples_variant[i])
      drawn_variant[[i]]<-interim[indices,1]
    }
    IDs<-c(drawn_variant[[1]],drawn_variant[[2]],drawn_variant[[3]])
    dat%>%
      filter(dat$ID_study %in% IDs)-> test_data
    dat%>%
      filter(!(dat$ID_study %in% IDs))-> train_data
    
    
    #expand the the train_data set to 1000 samples with equal fractions
    train_expansion<-as.list(1)
    for(i in 1:length(variants)){
      variant_number<-333-sum(train_data$Variant==variants[i])
      train_data%>%
        filter(Variant == variants[i])-> interim
      indices<-sample(nrow(interim), size = variant_number, replace = TRUE)
      train_expansion[[i]]<-interim[indices,]
    }
    
    
    expansion<-rbind(train_expansion[[1]],train_expansion[[2]],train_expansion[[3]])
    
    for(i in 1:nrow(expansion)){
      for(j in 5:10){
        distribution<-rnorm(50,mean = expansion[i,j],sd=(expansion[i,j]/10))
        expansion[i,j]<-distribution[runif(n=1,min=1,max=50)]
      }
    }
    train_dataset<-rbind(train_data,expansion)
    test_dataset<-test_data
    train_list[[x]]<-train_dataset
    test_list[[x]]<-test_dataset
    }
    visit_sets[[1]]<-train_list
    visit_sets[[2]]<-test_list
    return(visit_sets)
  }
