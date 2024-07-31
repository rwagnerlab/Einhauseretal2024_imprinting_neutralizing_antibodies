#script to transform data layout for suppl figure 1
##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
library(writexl)
library(tidyverse)

rm(list = ls())

load("./adjustment/workspace_adjusted_long_simple.Rda")


mapping <- c("FA" = 1, "UA" = 2, "FD" = 3, "UD" = 4, "FO" = 5, "UO" = 6)
VOI <- c("BQ","JN","XBB")#

export_simple$groupnum <- mapping[as.character(export_simple$group)]
export_simple <- subset(export_simple, export_simple$visit == "v5")

export_simple_VOI <- subset(export_simple, export_simple$virus %in% VOI)
export_simple_VOI <- export_simple_VOI %>%
  pivot_wider(names_from = virus, values_from = adjusted_neutralization)

sorted_VOI <- export_simple_VOI[order(export_simple_VOI$groupnum),]
sorted_VOI <- na.omit(sorted_VOI)
sorted_VOI <- sorted_VOI %>%
  pivot_wider(names_from = group, values_from = VOI)
if(!dir.exists("./supplfig1_BQJNXBB")){dir.create("./supplfig1_BQJNXBB")}
write_xlsx(sorted_VOI,"./supplfig1_BQJNXBB/values.xlsx")

