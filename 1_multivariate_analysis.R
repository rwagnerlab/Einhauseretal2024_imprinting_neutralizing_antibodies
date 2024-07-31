##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#multivariate analysis including smoking status, vaccine type, age, sex, time between infection and number of vaccinations to determine impact of infection variant on neutralization
library(readxl)
library(dplyr)
library(writexl)
library(tidyverse)
library(car)
library(vegan)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(plotly)
library(reshape2)
library(AER)
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

merged_data$vector[merged_data$vaccine == 0] <- NA

#now add the smoking, age , gender , number of  vaccination and time between vaccination information
df<-read_xlsx("./Age-sex-smoking.xlsx")
df$Gender<-as.factor(df$Gender)
df$Age<-as.numeric(df$Age)
df$Number_of_vaccinations<-as.numeric(df$Number_of_vaccinations)
df$Nicotine_consumption<-as.factor(df$Nicotine_consumption)
df$Timebetweenvaccination_infection<-as.numeric(df$Timebetweenvaccination_infection)
df$Variant_of_infection<-NULL
df$Group<-NULL
df$Year_of_birth<-NULL

#merge again
merged_data <- merged_data %>%
  left_join(df, by = c("ID_study" = "ID_Proband"))


#pivot long
merged_data <- merged_data %>%
  pivot_longer(cols = c("D614G","Alpha","Delta","BA1","BA2","BA5","BQ","XBB","JN"), names_to = "virus", values_to = "neutralization")

# Omit rows with NA values in the 'neutralization' column
merged_data_clean <- merged_data[complete.cases(merged_data$neutralization), ]

#tried log transform for normality but using a nonpar test instead
#merged_data_clean$neutralization<-log(merged_data_clean$neutralization)

# Reshape to wide format
merged_data_wide <- merged_data_clean %>%
  pivot_wider(names_from = visit, values_from = neutralization, names_prefix = "neutralization_visit")

shapiro.test(merged_data_wide$neutralization_visitv1)



# Prepare data for PERMANOVA it doesn't handle NAs so two seperate models are needed
merged_data_vax <- subset(merged_data_wide, merged_data_wide$vaccine == 1)
merged_data_unvax <- subset(merged_data_wide, merged_data_wide$vaccine == 0)

merged_data_vax <- na.omit(merged_data_vax)
#get rid of all nas for unvaccinated individuals by nonapplicable variables
merged_data_unvax$vector<-NULL
merged_data_unvax$no_vector<-NULL
merged_data_unvax$Timebetweenvaccination_infection<-NULL
merged_data_unvax$Number_of_vaccinations<-NULL
merged_data_unvax <- na.omit(merged_data_unvax)
# PERMANOVA requires a distance matrix.
distance_matrix_vax <- dist(merged_data_vax[, c("neutralization_visitv1", "neutralization_visitv2", "neutralization_visitv4", "neutralization_visitv5")], method = 'euclidean')

# Fit the PERMANOVA model
permanova_model_vax <- adonis2(distance_matrix_vax ~ Age + Gender + Nicotine_consumption + variant + Timebetweenvaccination_infection + Number_of_vaccinations + vector +virus, 
                          data = merged_data_vax)




# PERMANOVA requires a distance matrix.
distance_matrix_unvax <- dist(merged_data_unvax[, c("neutralization_visitv1", "neutralization_visitv2", "neutralization_visitv4", "neutralization_visitv5")], method = 'euclidean')

# Fit the PERMANOVA model
permanova_model_unvax <- adonis2(distance_matrix_unvax ~ Age + Gender + Nicotine_consumption + variant + virus, 
                               data = merged_data_unvax)

if(!dir.exists("./multivariate")){dir.create("./multivariate")}
# Print the summary of PERMANOVA results


#post hoc testing

variant <- merged_data_vax$variant
# Perform pairwise comparisons
pairwise_results_vax <- pairwise.adonis(distance_matrix_vax, factors = variant)


variant_unvax <- merged_data_unvax$variant
pairwise_results_unvax <- pairwise.adonis(distance_matrix_unvax, factors = variant_unvax)


###adjusted pairwise analysis 
##for vaccinated
vaccinated<-merged_data_clean[merged_data_clean$vaccine == 1,]
vaccinated$logneutralization <- log2(vaccinated$neutralization)
tobit_model <- tobit(logneutralization ~ visit + vector + Gender + Age + Timebetweenvaccination_infection + Number_of_vaccinations + Nicotine_consumption + virus, 
                       data = vaccinated, left = log2(1), right = log2(2561))

# Predict values using the Tobit model
fitted_values <- predict(tobit_model, newdata = vaccinated)
antilog_fitted <- 2^fitted_values
antilog_fitted <- pmin(pmax(antilog_fitted, 1), 2561)
actual_values <- vaccinated$neutralization
vaccinated$residuals <- actual_values - antilog_fitted

vaccinated<- na.omit(vaccinated)
vaccinated$neutralization<-NULL
vaccinated$logneutralization<-NULL

vaccinated_wide <- vaccinated %>%
  pivot_wider(names_from = visit, values_from = residuals, names_prefix = "residuals_")
# distance matrix adjusted
# Remove rows with any NA values in residuals
vaccinated_wide <- vaccinated_wide[complete.cases(vaccinated_wide[, c("residuals_v1", "residuals_v2", "residuals_v4", "residuals_v5")]), ]
# Make sure the variant column is properly referenced
variant <- vaccinated_wide$variant
#calc distances and do pairwise adonis
distance_matrix_vax_adj <- dist(vaccinated_wide[, c("residuals_v1", "residuals_v2", "residuals_v4", "residuals_v5")], method = 'euclidean')
pairwise_results_vax_adj <- pairwise.adonis(distance_matrix_vax_adj, factors = variant)


##for unvaccinated
unvaccinated<-merged_data_clean[merged_data_clean$vaccine == 0,]
unvaccinated$logneutralization <- log2(unvaccinated$neutralization)
tobit_model_unvax <- tobit(logneutralization ~ visit + Gender + Age + Nicotine_consumption + virus, 
                     data = unvaccinated, left = log2(1), right = log2(2561))

# Predict values using the Tobit model
fitted_values_uvx <- predict(tobit_model_unvax, newdata = unvaccinated)
antilog_fitted_uvx <- 2^fitted_values_uvx
antilog_fitted_uvx <- pmin(pmax(antilog_fitted_uvx, 1), 2561)
actual_values_uvx <- unvaccinated$neutralization
unvaccinated$residuals <- actual_values_uvx - antilog_fitted_uvx


unvaccinated$neutralization<-NULL
unvaccinated$logneutralization<-NULL

unvaccinated_wide <- unvaccinated %>%
  pivot_wider(names_from = visit, values_from = residuals, names_prefix = "residuals_")
# distance matrix adjusted
# Remove rows with any NA values in residuals
unvaccinated_wide <- unvaccinated_wide[complete.cases(unvaccinated_wide[, c("residuals_v1", "residuals_v2", "residuals_v4", "residuals_v5")]), ]
# Make sure the variant column is properly referenced
variant_uvx <- unvaccinated_wide$variant
distance_matrix_unvax_adj <- dist(unvaccinated_wide[, c("residuals_v1", "residuals_v2", "residuals_v4", "residuals_v5")], method = 'euclidean')
pairwise_results_unvax_adj <- pairwise.adonis(distance_matrix_unvax_adj, factors = variant_uvx)
###########export results
sink("./multivariate/modelssummary.txt")
print("Vaccinated Multivariate-Model:")
print(permanova_model_vax)
print("-------------------------------------------------------")
print("Unvaccinated Multivariate-Model:")
print(permanova_model_unvax)
print("========================================================")
print("Post Hoc testing")
# Print pairwise comparison results
print("Vaccinated Model unajdusted:")
print(pairwise_results_vax)
print("Unvaccinated Model unadjusted:")
print(pairwise_results_unvax)
sink()

sink("./multivariate/modelssummary_adj.txt")
print("vaccinated Model adjusted:")
print(pairwise_results_vax_adj)
print("Unvaccinated Model adjusted:")
print(pairwise_results_unvax_adj)
sink()
#ordination plots
# Compute PCoA
pcoa_result <- cmdscale(distance_matrix_vax, k = 3, eig = TRUE)

# Create the data frame for plotting
pcoa_data <- data.frame(
  SampleID = merged_data_vax$ID_study,
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  PC3 = pcoa_result$points[, 3],
  Variant = merged_data_vax$variant,
  Age = merged_data_vax$Age,
  Gender = merged_data_vax$Gender,
  Nicotine_consumption = merged_data_vax$Nicotine_consumption
)

# Plot PCoA
ggplot(pcoa_data, aes(x = PC1, y = PC2, color = Variant)) +
  geom_point(size = 3, alpha = 0.3) +
  labs(title = "PCoA of Neutralization Data",
       x = paste("PC1 (", round(pcoa_result$eig[1]/sum(pcoa_result$eig) * 100, 1), "%)", sep = ""),
       y = paste("PC2 (", round(pcoa_result$eig[2]/sum(pcoa_result$eig) * 100, 1), "%)", sep = "")) +
  theme_minimal()

# Create the PCoA plot with additional covariates
ggplot(pcoa_data, aes(x = PC1, y = PC2, color = Variant, size = Age, shape = Gender, fill = Nicotine_consumption)) +
  geom_point(alpha = 0.6) +  # Adjust alpha for transparency
  scale_shape_manual(values = c(16, 17)) +  # Custom shapes for gender (e.g., 16 for female, 17 for male)
  scale_fill_manual(values = c("never" = "blue","formerly" = "grey", "yes" = "red")) +  # Custom fill colors for nicotine consumption
  labs(title = "PCoA of Neutralization Data",
       x = paste("PC1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)", sep = ""),
       y = paste("PC2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)", sep = "")) +
  theme_minimal() +
  theme(legend.position = "right")  # Optionally adjust legend position


# Create a 3D scatter plot with plotly
# Create a 3D scatter plot with plotly
plot_ly(data = pcoa_data, 
        x = ~PC1, 
        y = ~PC2, 
        z = ~PC3, 
        color = ~Variant,  # Use color to differentiate between variants
        size = ~Age,  # Optionally adjust size based on age
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(size = 5, opacity = 0.6)) %>%  # Simplified marker settings
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')),
         title = '3D PCoA Plot of Neutralization Data')

# Compute NMDS
nmds_result <- metaMDS(distance_matrix_vax, k = 2, trymax = 100)

# Create a data frame for plotting
nmds_data <- data.frame(SampleID = rownames(nmds_result$points),
                        NMDS1 = nmds_result$points[, 1],
                        NMDS2 = nmds_result$points[, 2],
                        Variant = merged_data_vax$variant,
                        Age = merged_data_vax$Age,
                        Gender = merged_data_vax$Gender,
                        Nicotine_consumption = merged_data_vax$Nicotine_consumption)

# Plot NMDS
ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = Variant)) +
  geom_point(size = 3) +
  labs(title = "NMDS of Neutralization Data",
       x = "NMDS1",
       y = "NMDS2") +
  theme_minimal()




