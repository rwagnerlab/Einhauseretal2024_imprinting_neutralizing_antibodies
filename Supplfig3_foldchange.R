##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
### script to generate the fold change values

library(writexl)
library(tidyverse)
library(readxl)
rm(list = ls())

load("./adjustment/workspace_adjusted_long_simple.Rda")

# Calculate fold-change values with Alpha neut as reference for FA and UA and so on
data_foldchange <- export_simple %>%
  group_by(group) %>%
  mutate(
    Reference = case_when(
      group %in% c("FA", "UA") & virus == "Alpha" ~ adjusted_neutralization,
      group %in% c("FD", "UD") & virus == "Delta" ~ adjusted_neutralization,
      group %in% c("FO", "UO") & virus == "BA1" ~ adjusted_neutralization,
      TRUE ~ NA_real_
    )
  ) %>%
  fill(Reference, .direction = "downup") %>%
  mutate(FoldChange = adjusted_neutralization / Reference) %>%
  ungroup() %>%
  select(-Reference)



# Calculate median, 25th percentile, and 75th percentile of fold-change for each group, visit, and virus
summary_stats <- data_foldchange %>%
  group_by(group, visit, virus) %>%
  summarise(
    MedianFoldChange = median(FoldChange, na.rm = TRUE),
    Q1FoldChange = quantile(FoldChange, 0.25, na.rm = TRUE),
    Q3FoldChange = quantile(FoldChange, 0.75, na.rm = TRUE)
  ) %>%
  ungroup()
summary_stats <- na.omit(summary_stats)

#split summary stats into the single groups
groups <- unique(summary_stats$group)
stats <- vector('list', length = length(groups))
for (i in 1:length(groups)) {
  stat <- subset(summary_stats, summary_stats$group == groups[i])
  stat <- stat %>%
    select(-group)
  
  stat_wide <- stat %>%
    pivot_wider(names_from = visit, values_from = c(MedianFoldChange, Q1FoldChange, Q3FoldChange))
  
  # Get the current column names
  cols <- colnames(stat_wide)
  
  # Create a new order of columns
  visits <- unique(summary_stats$visit)
  new_order <- c("virus")
  for (visit in visits) {
    new_order <- c(new_order, paste0("MedianFoldChange_", visit), paste0("Q1FoldChange_", visit), paste0("Q3FoldChange_", visit))
  }
  
  # Reorder columns
  stat_wide <- stat_wide %>%
    select(all_of(new_order))
  
  # Sort rows by Virus
  stat_wide <- stat_wide %>%
    arrange(factor(virus, levels = c("D614G", "Alpha", "Delta", "BA1", "BA2", "BA5", "BQ")))
  
  stats[[i]] <- stat_wide
}



#export
if(!dir.exists("./supplfig2_foldchange")){
  dir.create("./supplfig2_foldchange/")
}
for(i in 1:length(groups)){
  write_xlsx(stats[[i]],paste0("./supplfig2_foldchange/", groups[i],".xlsx"))
}
