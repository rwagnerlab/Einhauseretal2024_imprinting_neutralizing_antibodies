##by Sebastian Einhauser for Einhauser et al. 
magbreadth <- function(df)
{
  df <- subset(df, select = c(group, adjusted_neutralization))
  df <- na.omit(df)
  # 1. Count the number of neutralization values per group
  group_counts <- table(df$group)
  
  # List of numbers from 1 to 2561
  numbers <- 0:2561
  
  # 2. Loop over the list of numbers and count neutralization values >= current number for each group
  results <- list()
  for (num in numbers) 
  {
    counts <- sapply(split(df$adjusted_neutralization, df$group), function(neutralization_values) {
      sum(neutralization_values >= num)
    })
    results[[as.character(num)]] <- counts
  }
  
  # 3.fuse the generated lists together to a data frame
  mb_counts <-as.data.frame(results[[1]])
  for( i in 2: length(numbers))
  {
   mb_counts <- cbind(mb_counts, as.data.frame(results[[i]]))    
  }
  colnames(mb_counts) <- as.character(numbers)
  
  mb_counts <- as.data.frame(t(mb_counts))
  
  #transform the counts into percentages
  groups<-unique(df$group)
  for ( i in 1: length(groups))
  {
    mb_counts[,i]<-100*(mb_counts[,i]/group_counts[i])
  }
return(mb_counts)
}


