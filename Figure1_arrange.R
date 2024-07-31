##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#script to transform data layout for figure 1
library(writexl)
library(tidyverse)

rm(list = ls())

load("./adjustment/workspace_adjusted_long_simple.Rda")


mapping <- c("FA" = 1, "UA" = 2, "FD" = 3, "UD" = 4, "FO" = 5, "UO" = 6)

export_simple$groupnum <- mapping[as.character(export_simple$group)]
export_simple$visgroup <- paste(export_simple$groupnum,export_simple$visit, sep = "")

df_sorted <- export_simple %>%
  arrange(ID_study) %>%
  arrange(visgroup)


viruses<- unique(df_sorted$virus)

dataprocessed<-as.list(NA)
for( i in 1:length(viruses)){
  interim <- subset(df_sorted, virus == viruses[i])
  interim <- interim %>%
      pivot_wider(names_from = visgroup, values_from = adjusted_neutralization)
  dataprocessed[[i]] <- interim
}

dir.create("./fig1")

for( i in 1:length(viruses)){
  
    interim <- dataprocessed[[i]]
    write_xlsx(interim, paste0("./fig1/",viruses[i],".xlsx"))
}

