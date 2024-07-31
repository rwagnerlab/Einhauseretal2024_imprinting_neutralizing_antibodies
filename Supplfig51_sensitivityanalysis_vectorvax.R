##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#sensitivity analysis to determine the impact of vector vaccines
library(readxl)
library(dplyr)
library(writexl)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggsignif)
rm(list=ls())
load("./adjustment/workspace_adjusted_wide_simple.Rda")

vaccine_type<-read_excel("./vaccinations.xlsx")

#adjust the format to workable data
vaccine_type <- vaccine_type %>%
  mutate(across(5:21, ~ case_when(
    is.na(.) ~ "0",
    . == "x" ~ "1",
    TRUE ~ .
  )))

# Convert the specified columns to numeric
vaccine_type <-  vaccine_type %>%
  mutate(across(5:21, as.numeric))

#get all the vector vaccines
# Select columns whose names contain "a", "j", or both
selected_cols <- grep("[AJ]", names(vaccine_type), value = TRUE)

# Calculate the sum across rows for selected columns
vaccine_type$vector<- rowSums(vaccine_type[selected_cols])

#all vaccines
vaccine_type<-vaccine_type %>%
  mutate(sumall = rowSums(across(5:21)))
#substract the vector vaccines
vaccine_type$no_vector <- vaccine_type$sumall - vaccine_type$vector

#take care of weird column name
names<-colnames(vaccine_type)
names[1]<-"ID_Proband"
colnames(vaccine_type)<- names

#now select the dataof interest
vaccine_type <- vaccine_type %>%
  select(all_of(c("ID_Proband", "vector", "no_vector")))


merged_data <- exportwide_simple %>%
  left_join(vaccine_type, by = c("ID_study" = "ID_Proband"))

merged_data <- merged_data %>%
  pivot_longer(cols = c("D614G","Alpha","Delta","BA1","BA2","BA5","BQ","XBB","JN"), names_to = "virus", values_to = "neutralization")

merged_data <- na.omit(merged_data)
#count the vector group
unique(merged_data$group)
countingset<-subset(merged_data, virus == "D614G" & visit == "v1" & vaccine == 1)
sum(countingset$vector)

merged_data<- subset(merged_data, vaccine == 1)
#start plotting and statistics

visit_labels <- c(v1 = "Visit 1", v2 = "Visit 2", v4 = "Visit 4", v5 = "Visit 5")

figureDir<-paste0(getwd(),"/sensitivityanalysis_vector")
if(!dir.exists(figureDir)){dir.create(figureDir)}
png(paste0(figureDir,"/sensitivity_analysis_vector.png"), res = 300, width = 3000, height = 1600)
ggplot(merged_data, aes(x = virus, y = neutralization, fill = factor(vector))) +
  geom_boxplot() +
  facet_wrap(~ visit, scales = "free", labeller = labeller(visit = visit_labels)) +
  labs(
    x = "Neutralized Virus",
    y = "Neutralization IC 50",
    fill = "Vector-Vaccine",
    title = "Impact of Vector-Vaccines on Variant Specific Neutralization"
  ) +
  scale_y_log10() +
  theme_minimal()
dev.off()


# Define function to perform Wilcoxon test
wilcoxon_test <- function(df) {
  test_result <- tryCatch({
    test <- wilcox.test(neutralization ~ vector, data = df, exact = FALSE)
    p_value <- test$p.value
    list(p_value = p_value)
  }, error = function(e) {
    list(p_value = NA)
  })
  
  return(test_result)
}

# Group by visit and virus, then perform unpaired Wilcoxon test for each group
test_results <- merged_data %>%
  group_by(visit, virus) %>%
  nest() %>%
  mutate(test_result = map(data, wilcoxon_test)) %>%
  unnest_wider(test_result)

# Print test results
print(test_results)

write_xlsx(test_results, path = paste0(figureDir,"/multitests.xlsx"))

#second option perform tests on pooled neutralization for all viruses
test_results2 <- merged_data %>%
  group_by(visit) %>%
  nest() %>%
  mutate(test_result = map(data, wilcoxon_test)) %>%
  unnest_wider(test_result)
write_xlsx(test_results2, path = paste0(figureDir,"/pooledtests.xlsx"))


#plot for the pooled analysis
png(paste0(figureDir,"/sensitivity_analysis_vector_pooled.png"), res = 300, width = 3000, height = 1600)
ggplot(merged_data, aes(x = visit, y = neutralization, fill = factor(vector))) +
  geom_boxplot() +
  labs(
    x = "Visit",
    y = "Neutralization IC50",
    fill = "Vector-Vaccine",
    title = "Impact of Vector-Vaccines on Overall Neutralization"
  ) +
  scale_y_log10() +
  theme_minimal()

dev.off()

# Add y_position for the significance annotations
max_values <- merged_data %>%
  group_by(visit, virus) %>%
  summarize(max_value = max(neutralization)) %>%
  ungroup()

test_results <- test_results %>%
  left_join(max_values, by = c("visit", "virus")) %>%
  mutate(y_position = max_value * 1.1)  # Adjust y_position as needed

max_values2 <- merged_data %>%
  group_by(visit) %>%
  summarize(max_value = max(neutralization)) %>%
  ungroup()
test_results2 <- test_results2 %>%
  left_join(max_values2, by = c("visit")) %>%
  mutate(y_position = max_value *1.1) #adjust if needed



#generate one large plot with both subplots
merged_data <- merged_data %>%
  mutate(
    vector = factor(vector, levels = c(0, 1)),
    visit = factor(visit),
    virus = factor(virus)
  ) %>%
  drop_na(vector, neutralization, visit, virus)
# Define custom legend labels
vector_labels <- c("0" = "No Vector", "1" = "Vector")

# Generate the first plot
plot_a <- ggplot(merged_data, aes(x = factor(vector), y = neutralization, fill = factor(vector))) +
  geom_boxplot() +
  facet_wrap(~ visit, scales = "free", labeller = labeller(visit = visit_labels)) +
  labs(
    x = "Vector-Vaccine",
    y = "Neutralization IC50",
    fill = "Vector-Vaccine",
    title = "Impact of Vector-Vaccines on Overall Neutralization"
  ) +
  scale_fill_manual(values = c("0" = "#2fbfc4", "1" = "#e96c6c"), labels = vector_labels) +
  scale_y_log10() +
  theme_minimal() +
  geom_signif(
    comparisons = list(c("0", "1")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    textsize = 3,
    vjust = 0
  )

# Generate the second plot
plot_b <- ggplot(merged_data, aes(x = virus, y = neutralization, fill = factor(vector))) +
  geom_boxplot() +
  facet_wrap(~ visit, scales = "free", labeller = labeller(visit = visit_labels)) +
  labs(
    x = "Neutralized Virus",
    y = "Neutralization IC50",
    fill = "Vector-Vaccine",
    title = "Impact of Vector-Vaccines on Variant Specific Neutralization"
  ) +
  scale_fill_manual(values = c("0" = "#2fbfc4", "1" = "#e96c6c"), labels = vector_labels) +
  scale_y_log10() +
  theme_minimal()  
  

# Combine the plots
combined_plot <- plot_grid(plot_a, plot_b, ncol = 1, align = 'v', rel_heights = c(1.5, 1))

# Add labels "A" and "B"
labeled_plot <- ggdraw() +
  draw_plot(combined_plot, 0, 0, 1, 1) +
  draw_plot_label(label = c("A", "B"), x = c(0, 0), y = c(1, 0.4), hjust = -0.5, vjust = 1.5, size = 14)

# Save the combined and labeled plot to a file
png(paste0(figureDir,"/combined_sensitivity_analysis_vector.png"), res = 300, width = 4500, height = 4000)
print(labeled_plot)
dev.off()
