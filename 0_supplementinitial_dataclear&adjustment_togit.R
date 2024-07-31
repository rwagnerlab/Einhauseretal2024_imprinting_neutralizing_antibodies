##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#clear workspace and load packages
rm(list=ls())

library(readxl)
library(dplyr)
library(tidyverse)
library(censReg)
library(geepack)
library(AER)
library(ggplot2)
library(patchwork)
library(tidyr)


`%notin%` <- Negate(`%in%`)
#make directories
if(!dir.exists("./adjustment")){dir.create("./adjustment")}
###read data#####
#read in data and assign variable scales
data <- as.data.frame(read_excel("./NT_Covako_complete_long.xlsx"))
colnames(data) <- c('ID_num','ID_study','visit','D614G','Alpha','Delta','BA1','BA2','BA5','BQ','XBB','JN','variant','vaccine','group','study_center')
data$visit <- as.factor(data$visit)

#make neutralization data numeric
for(i in 4:12){
  data[,i]<-as.numeric(data[,i])
}
data$variant<-as.factor(data$variant)
data$vaccine<-as.factor(data$vaccine)
data$group<-as.factor(data$group)

data$ID_num<-as.factor(data$ID_num)
#######prepare the data########
#subset the groups of interest and make long for the adjustment model
#subset to remove unused groups
workspace <- subset(data,!(group %in% c("PO", "PA", "PD"))) 
#remove unused factor levels
workspace$group <- as.character(workspace$group)
workspace$group <- as.factor(workspace$group)
#make long
workspace_long <- workspace %>% pivot_longer(
  cols = c('D614G','Alpha','Delta','BA1','BA2','BA5','BQ','XBB','JN'),
  names_to = "virus",
  values_to = "neutralization"
)
##############
#workspace_long$virus <- as.factor(workspace_long$virus)

#read in data to retrieve the days_visit##########
days<-as.data.frame(read_excel("./Covako_Visits.xlsx"))
days<-days[,c(3,16:19)]

colnames(days)<-c('ID_study','v1','v2','v4','v5')
days$v1<-as.numeric(days$v1)
days$v2<-as.numeric(days$v2)
days$v4<-as.numeric(days$v4)
days$v5<-as.numeric(days$v5)

days<- days %>%
  pivot_longer(
    cols = c('v1','v2','v4','v5'),
    names_to = 'visit',
    values_to = 'days_visit'
  ) 

#merge
merged<-merge(workspace_long, days , by = c("ID_study","visit"), all = FALSE)

#print(is.na(merged$days_visit))
sum(is.na(merged$days_visit))

####take care of missing days_visit by imputing with the mean days for that visit.####
mean_days_per_visit <- merged %>%
  group_by(visit) %>%
  summarise(mean_days = mean(days_visit, na.rm = TRUE))

# Join the mean_days back to the original dataframe
data_with_means <- left_join(merged, mean_days_per_visit, by = "visit")
data_with_means <- as.data.frame(data_with_means)
# Replace NA values in the days column with the mean value for the corresponding visit
data_with_means <- data_with_means %>%
  mutate(days_visit = ifelse(is.na(days_visit), mean_days, days_visit)) %>%
  select(-mean_days) # Remove the mean_days column as no longer needed

#clean up
workspace_long<-data_with_means
rm(data_with_means)
rm(merged)
rm(days)
rm(data)
rm(mean_days_per_visit)
##prepare data for fit########

##########log2 transform neutralization data
workspace_long$log_neutralization <- log2(workspace_long$neutralization)

#create helper variables for proper study center adjustment
workspace_long <- workspace_long%>%
  pivot_wider(names_from = study_center,values_from = study_center)%>%
  rename_at(vars(`1`:`6`), ~ paste0("studycenter_", .))

workspace_long <- workspace_long %>%
  mutate_at(vars(studycenter_1:studycenter_6), ~ replace_na(., 0))%>%
  mutate_at(vars(studycenter_1:studycenter_6), ~ ifelse(. > 0, 1, .))


workspace_long$days_visit_unscaled <- workspace_long$days_visit

# Group by visit and scale the days_visit_unscaled variable
workspace_long <- workspace_long %>%
  group_by(visit) %>%
  mutate(days_visit_scaled = as.numeric(scale(days_visit_unscaled)))

#clean up duplicates
workspace_long$days_visit<-workspace_long$days_visit_scaled
workspace_long<-workspace_long%>%
  select(-days_visit_scaled)

#bring the vaccination data in the right format
workspace_long$vaccine <- as.character(workspace_long$vaccine)
workspace_long$vaccine <- as.factor(workspace_long$vaccine)
workspace_long$vaccine <- as.numeric(workspace_long$vaccine)
workspace_long$vaccine <- (workspace_long$vaccine -2)*-1
#clear NAs in neutralization
model_data<- na.omit(workspace_long)


########multiple tobit models########
#since neutralization is left and right censored a tobit model is used for adjustments
#since effects might not be the same for each visit, either one complex model is needed
#or: 4 separate models one for each visit.
#split the data to the 4 visits
visits <- unique(model_data$visit)
model_data_visits <- as.list(NA)
model_coefs<- as.list(NA)
tobit_models <- as.list(NA)
#loop over the visits and create the model
for(i in 1:length(visits)){
  #subset the data
  visit_data <- as.data.frame(subset(model_data, visit == visits[i]))
  
  #fit a model
  ##exclude studycenter 5 ( most participants) to remove exclusive redundancy and correct other studycenter effects dependent on center 5
  tobit_model <- tobit(log_neutralization ~ group*days_visit + virus + studycenter_1 + studycenter_2 + studycenter_3 + studycenter_4 + studycenter_6,
                         left = log2(1),
                         right = log2(2561),
                         data = visit_data)
  #extract the coefficients
  coefs<-coef(tobit_model)
  #predict neutralization values
  visit_data$predicted <- predict(tobit_model, visit_data)
  #save to lists
  model_coefs[[i]]<-coefs
  tobit_models[[i]]<-tobit_model
  model_data_visits[[i]]<-visit_data
  
  #clean up
  rm(visit_data)
  rm(tobit_model)
}
#summarize the models
summary(tobit_models[[1]])
summary(tobit_models[[2]])
summary(tobit_models[[3]])
summary(tobit_models[[4]])


########multiple tobit plots and checks#############

# check out visit models now

model_data_visit_combined<- rbind(
  model_data_visits[[1]],
  model_data_visits[[2]],
  model_data_visits[[3]],
  model_data_visits[[4]]
)

#plot effect of days visit on neutralization
ggplot(model_data_visit_combined, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = visit), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = visit), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs. Days Visit by Group") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "Visit") +
  scale_x_continuous(limits = c(min(model_data_visit_combined$days_visit), max(model_data_visit_combined$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed

#plot the model vs the measured data
ggplot(model_data_visit_combined,aes(x = log_neutralization, y = predicted)) +
  geom_point(aes(color = visit), alpha = 0.1) +
  geom_smooth( aes(y = predicted, color = visit), method = "glm") +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs Model Prediction") +
  theme_minimal()+
  labs(x = "Log Neutralization", y = "Prediction", color = "Visit") +
  scale_x_continuous(limits = c(0, 13)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed

#predicted vs measured all
cor(model_data_visit_combined$log_neutralization,model_data_visit_combined$predicted)
plot(model_data_visit_combined$log_neutralization,model_data_visit_combined$predicted)
abline(0,1, col = "red")



########single tobit#########
#fit a single model over all visits as a competitor

tobit_model_large <- tobit(log_neutralization ~ group*days_visit + virus + studycenter_1*visit+studycenter_2*visit+studycenter_3*visit+studycenter_4*visit+studycenter_6*visit,
                           left = log2(1),
                           right = log2(2561),
                           data = model_data)

summary(tobit_model_large)
model_data$predicted <- predict(tobit_model_large, newdata = model_data)

#compare single visit model with the large model
sorted_model_data_visit_combined<- model_data_visit_combined[order(model_data_visit_combined$ID_num,model_data_visit_combined$virus), ]
sorted_model_data<-model_data[order(model_data$ID_num,model_data$virus), ]
plot(sorted_model_data$predicted,sorted_model_data_visit_combined$predicted)
abline(0,1, col = "red")
#plot effect of days visit on neutralization
ggplot(model_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = visit), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = visit), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs. Days Visit by Group") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "Visit") +
  scale_x_continuous(limits = c(0, 300)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed
#plot the model vs the measured data
ggplot(model_data,aes(x = log_neutralization, y = predicted)) +
  geom_point(aes(color = visit), alpha = 0.1) +
  geom_smooth( aes(y = predicted, color = visit), method = "glm") +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs Model Prediction")+
  labs(x = "Log Neutralization", y = "Prediction", color = "Visit") +
  scale_x_continuous(limits = c(0, 13)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed
#Log-likelihood: -9660 worse compared to appx. -2500 for single models
#plot overall model vs measured
ggplot(model_data,aes(x = log_neutralization, y = predicted))+  
  geom_point(aes(color = visit), alpha = 0.1) +
  geom_smooth( aes(y = predicted, color = visit), method = "glm") +
  ggtitle("Log Neutralization vs Model Prediction")+
  labs(x = "Log Neutralization", y = "Prediction", color = "Visit") +
  scale_x_continuous(limits = c(0, 13)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))
cor(model_data$log_neutralization,model_data$predicted)

#single visit model effects figure ##########
plot_data <- model_data_visits[[1]]
pv1<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs. Days Visit by Group") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "Visit") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed
plot_data <- model_data_visits[[2]]
pv2<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs. Days Visit by Group") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "Visit") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed
plot_data <- model_data_visits[[3]]
pv3<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs. Days Visit by Group") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "Visit") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed
plot_data <- model_data_visits[[4]]
pv4<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Log Neutralization vs. Days Visit by Group") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "Visit") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))     # Adjust y-axis limits as needed
large_plot <- pv1/pv2/pv3/pv4+ plot_layout(ncol = 2)

png(filename = "./adjustment/adjustment_effects.png", res = 300, width = 3200, height = 1800)
large_plot
dev.off()



####adjust the data for studycenter effects and the visit time differences########
###use the single visit tobit models for adjustment
#adjust workspace by visit and append again
adjustedworkspace<-as.list(NA)

for(i in 1:length(visits)){
  visit_data <- subset(workspace_long, visit == visits[i])
  coef_estimates <- model_coefs[[i]]
  
  # Adjust the response variable for study center and days_visit effects
  
  visit_data <- visit_data %>%
    mutate(
      adjusted_neutralization = log_neutralization -
        coef_estimates["days_visit"] * days_visit -
        coef_estimates["studycenter_1"] * as.numeric(studycenter_1)-
        coef_estimates["studycenter_2"] * as.numeric(studycenter_2)-
        coef_estimates["studycenter_3"] * as.numeric(studycenter_3)-
        coef_estimates["studycenter_4"] * as.numeric(studycenter_4)-
        coef_estimates["studycenter_6"] * as.numeric(studycenter_6)
    )
  
  # Adjust for interaction terms
  for (grp in unique(visit_data$group)) {
    interaction_term <- paste0("group", grp, ":days_visit")
    if (interaction_term %in% names(coef_estimates)) {
      visit_data <- visit_data %>%
        mutate(
          adjusted_neutralization = adjusted_neutralization -
            coef_estimates[interaction_term] * (days_visit * (group == grp))
        )
    }
  }
  adjustedworkspace[[i]]<-visit_data
}  


#bind to a single data frame
adjustedworkspaceall<-rbind(adjustedworkspace[[1]],adjustedworkspace[[2]],adjustedworkspace[[3]],adjustedworkspace[[4]])

##############antilog the adjusted neutralization############
adjustedworkspaceall$adjusted_neutralization_log <- adjustedworkspaceall$adjusted_neutralization
adjustedworkspaceall <- adjustedworkspaceall %>%
    mutate(
      adjusted_neutralization = 2^adjusted_neutralization
    )


#set assay constraints of 1 and 2561 to the adjusted values
min(adjustedworkspaceall$adjusted_neutralization, na.rm = TRUE)
max(adjustedworkspaceall$adjusted_neutralization, na.rm = TRUE)
sum(is.na(adjustedworkspaceall$adjusted_neutralization))

adjustedworkspaceall$adjusted_neutralization  <- pmax(adjustedworkspaceall$adjusted_neutralization, 1)
adjustedworkspaceall$adjusted_neutralization  <- pmin(adjustedworkspaceall$adjusted_neutralization, 2561)

min(adjustedworkspaceall$adjusted_neutralization, na.rm = TRUE)
max(adjustedworkspaceall$adjusted_neutralization, na.rm = TRUE)
sum(is.na(adjustedworkspaceall$adjusted_neutralization))

plot(adjustedworkspaceall$neutralization,adjustedworkspaceall$adjusted_neutralization)                  
cor(model_data_visit_combined$log_neutralization,model_data_visit_combined$predicted, method = "spearman")

#plots to compare adjusted with measured by visit and group##########

plot_data <- subset(adjustedworkspaceall, adjustedworkspaceall$visit=="v1")
p1<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 1: Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

plot_data <- subset(adjustedworkspaceall, adjustedworkspaceall$visit=="v2")
p2<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 2:Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

plot_data <- subset(adjustedworkspaceall, adjustedworkspaceall$visit=="v4")
p4<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 4:Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

plot_data <- subset(adjustedworkspaceall, adjustedworkspaceall$visit=="v5")
p5<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 5:Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

large_plot_measuredvsadjust <- p1/p2/p4/p5+ plot_layout(ncol = 2)
png("./adjustment/adjustedvsmeasured.png", res = 300, width = 3200, height = 1800)
  large_plot_measuredvsadjust
dev.off()
adjustedworkspaceall %>%
  filter(group == "UA", visit == "v5")%>%
  summarize(max_adjusted_neutralization = max(adjusted_neutralization, na.rm = TRUE))

#########save adjusted data########
#save the created and cleaned data
export <- adjustedworkspaceall
export$log_neutralization <- NULL
export$adjusted_neutralization_log <- NULL
export$days_visit <- NULL
export$neutralization <- NULL
export$ID_num <- NULL
export$days_visit_unscaled <- NULL
export$studycenter_1 <- NULL
export$studycenter_2 <- NULL
export$studycenter_3 <- NULL
export$studycenter_4 <- NULL
export$studycenter_5 <- NULL
export$studycenter_6 <- NULL
#make wide again if needed later
exportwide <- export %>%
  pivot_wider(names_from = virus, values_from = adjusted_neutralization)

#save the dataframes
save(exportwide, file = "./adjustment/workspace_adjusted_wide.Rda")
save(export, file = "./adjustment/workspace_adjusted_long.Rda")



#applying a more simplified model less dependent on super small groups just using the vaccination status instead of the full group#############
model_coefs_simple<- as.list(NA)
tobit_models_simple <- as.list(NA)
model_data_visits_simple <- as.list(NA)
visits <- unique(model_data$visit)


#loop over the visits and create the model
for(i in 1:length(visits)){
  #subset the data
  visit_data <- as.data.frame(subset(model_data, visit == visits[i]))
  visit_data$predicted <- NULL
  #fit a model
  ##exclude studycenter 5 ( most participants) to remove exclusive redundancy and correct other studycenter effects dependent on center 5
  tobit_model <- tobit(log_neutralization ~ vaccine *days_visit + virus + studycenter_1 + studycenter_2 + studycenter_3 + studycenter_4 + studycenter_6,
                       left = log2(1),
                       right = log2(2561),
                       data = visit_data)
 
  #extract the coefficients
  coefs<-coef(tobit_model)
  #predict neutralization values
  visit_data$predicted <- predict(tobit_model, visit_data)
  #save to lists
  model_coefs_simple[[i]]<-coefs
  tobit_models_simple[[i]]<-tobit_model
  model_data_visits_simple[[i]]<-visit_data
  
  #clean up
  rm(visit_data)
  rm(tobit_model)
}
#summarize the models
summary(tobit_models_simple[[1]])
summary(tobit_models_simple[[2]])
summary(tobit_models_simple[[3]])
summary(tobit_models_simple[[4]])

########multiple simpler tobit plots and checks#############

# check out visit models now

model_data_visits_simple_combined<- rbind(
  model_data_visits_simple[[1]],
  model_data_visits_simple[[2]],
  model_data_visits_simple[[3]],
  model_data_visits_simple[[4]]
)
#output model summaries to a text file
sink("./adjustment/models_summary.txt")
print("overall R2 pearson:")
cor(model_data_visits_simple_combined$log_neutralization, model_data_visits_simple_combined$predicted)
print("single model summaries")
summary(tobit_models_simple[[1]])
summary(tobit_models_simple[[2]])
summary(tobit_models_simple[[3]])
summary(tobit_models_simple[[4]])
sink()

###use the simpler single visit tobit models for adjustment
#adjust workspace by visit and append again
adjustedworkspace_simple<-as.list(NA)


for(i in 1:length(visits)){
  visit_data <- subset(workspace_long, visit == visits[i])
  coef_estimates <- model_coefs_simple[[i]]
  
  # Adjust the response variable for study center and days_visit effects
  
  visit_data <- visit_data %>%
    mutate(
      adjusted_neutralization = log_neutralization -
        coef_estimates["days_visit"] * days_visit -
        coef_estimates["studycenter_1"] * as.numeric(studycenter_1)-
        coef_estimates["studycenter_2"] * as.numeric(studycenter_2)-
        coef_estimates["studycenter_3"] * as.numeric(studycenter_3)-
        coef_estimates["studycenter_4"] * as.numeric(studycenter_4)-
        coef_estimates["studycenter_6"] * as.numeric(studycenter_6)-
        coef_estimates["vaccine:days_visit"] * days_visit * as.numeric(vaccine)
    )
  
  #Adjust for interaction terms
#  for (grp in unique(visit_data$group)) {
#    interaction_term <- paste0("group", grp, ":days_visit")
#    if (interaction_term %in% names(coef_estimates)) {
#      visit_data <- visit_data %>%
#        mutate(
#          adjusted_neutralization = adjusted_neutralization -
#            coef_estimates[interaction_term] * (days_visit * (group == grp))
#        )
#    }
#  }
  adjustedworkspace_simple[[i]]<-visit_data
}  


#bind to a single data frame
adjustedworkspaceall_simple<-rbind(adjustedworkspace_simple[[1]],adjustedworkspace_simple[[2]],adjustedworkspace_simple[[3]],adjustedworkspace_simple[[4]])

############## finally antilog the adjusted neutralization############
adjustedworkspaceall_simple$adjusted_neutralization_log <- adjustedworkspaceall_simple$adjusted_neutralization
adjustedworkspaceall_simple <- adjustedworkspaceall_simple %>%
  mutate(
    adjusted_neutralization = 2^adjusted_neutralization_log
  )


#set assay constraints of 1 and 2561 to the adjusted values
min(adjustedworkspaceall_simple$adjusted_neutralization, na.rm = TRUE)
max(adjustedworkspaceall_simple$adjusted_neutralization, na.rm = TRUE)
sum(is.na(adjustedworkspaceall_simple$adjusted_neutralization))

adjustedworkspaceall_simple$adjusted_neutralization  <- pmax(adjustedworkspaceall_simple$adjusted_neutralization, 1)
adjustedworkspaceall_simple$adjusted_neutralization  <- pmin(adjustedworkspaceall_simple$adjusted_neutralization, 2561)


#adjustment plot for the simple model###########

plot_data <- subset(adjustedworkspaceall_simple, adjustedworkspaceall_simple$visit=="v1")
ps1<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 1: Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

plot_data <- subset(adjustedworkspaceall_simple, adjustedworkspaceall_simple$visit=="v2")
ps2<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 2:Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

plot_data <- subset(adjustedworkspaceall_simple, adjustedworkspaceall_simple$visit=="v4")
ps4<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 4:Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

plot_data <- subset(adjustedworkspaceall_simple, adjustedworkspaceall_simple$visit=="v5")
ps5<-ggplot(plot_data, aes(x = neutralization, y = adjusted_neutralization)) +
  geom_point(aes(color = group), alpha = 0.1) +
  geom_smooth(aes(y = adjusted_neutralization, color = group), method = "loess", se = TRUE) +
  facet_wrap(~group) +
  ggtitle("Visit 5:Adjusted Neutralization vs. Measured Neutralization") +
  theme_minimal() +
  labs(x = "Measured Neutralization", y = "Adjusted Neutralization", color = "Group") +
  scale_x_continuous(limits = c(0, 2562)) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 2562))     # Adjust y-axis limits as needed

large_plot_measuredvsadjust_simple <- ps1/ps2/ps4/ps5+ plot_layout(ncol = 2)
png("./adjustment/adjustedvsmeasured_simplemodel.png", res = 300, width = 3200, height = 1800)
large_plot_measuredvsadjust_simple
dev.off()
#effects plot for the simplified model###########
plot_data <- subset(model_data_visits_simple_combined, visit == "v1")
p1<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = as.factor(vaccine)), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = as.factor(vaccine)), method = "lm", se = TRUE) +
  ggtitle("Visit 1: Log Neutralization vs. Days Visit by vaccine") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "vaccine") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))
plot_data <- subset(model_data_visits_simple_combined, visit == "v2")
p2<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = as.factor(vaccine)), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = as.factor(vaccine)), method = "lm", se = TRUE) +
  ggtitle("Visit 2: Log Neutralization vs. Days Visit by vaccine") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "vaccine") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))
plot_data <- subset(model_data_visits_simple_combined, visit == "v4")
p4<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = as.factor(vaccine)), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = as.factor(vaccine)), method = "lm", se = TRUE) +
  ggtitle("Visit 4: Log Neutralization vs. Days Visit by vaccine") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "vaccine") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))
plot_data <- subset(model_data_visits_simple_combined, visit == "v5")
p5<-ggplot(plot_data, aes(x = days_visit, y = log_neutralization)) +
  geom_point(aes(color = as.factor(vaccine)), alpha = 0.1) +
  geom_smooth(aes(y = predicted, color = as.factor(vaccine)), method = "lm", se = TRUE) +
  ggtitle("Visit 5: Log Neutralization vs. Days Visit by vaccine") +
  theme_minimal() +
  labs(x = "Days Visit", y = "Log Neutralization", color = "vaccine") +
  scale_x_continuous(limits = c(0.9*min(plot_data$days_visit), 1.1*max(plot_data$days_visit))) +  # Adjust x-axis limits as needed
  scale_y_continuous(limits = c(0, 13))

large_plot_vax <- p1/p2/p4/p5+ plot_layout(ncol = 2)

png("./adjustment/adjustedeffects_simplemodel.png", res = 300, width = 3200, height = 1800)

large_plot_vax

dev.off()

png("./adjustment/MeasuredvsModel.png", res = 300, width = 3200, height = 1800)
ggplot(model_data_visits_simple_combined,aes(x = log_neutralization, y = predicted))+
  geom_point(aes(color = visit), alpha = 0.1)+
  geom_smooth(aes(y = predicted, color = visit), method = "lm", se = TRUE)+
  ggtitle("Log Neutralization vs. Model Prediction") +
  theme_minimal() +
  labs(x = "Log Neutralization", y = "Model predicted log neutralization", color = "Visit") +
  scale_x_continuous(limits = c(0, 13))+
  scale_y_continuous(limits = c(0, 13))
dev.off()
# finally export the simple model data###########

export_simple <- adjustedworkspaceall_simple
export_simple$log_neutralization <- NULL
export_simple$adjusted_neutralization_log <- NULL
export_simple$days_visit <- NULL
export_simple$neutralization <- NULL
export_simple$ID_num <- NULL

export_simple$studycenter_1 <- NULL
export_simple$studycenter_2 <- NULL
export_simple$studycenter_3 <- NULL
export_simple$studycenter_4 <- NULL
export_simple$studycenter_5 <- NULL
export_simple$studycenter_6 <- NULL

export_simple$days_visit_unscaled <- NULL
#make wide again if needed later
exportwide_simple <- export_simple %>%
  pivot_wider(names_from = virus, values_from = adjusted_neutralization)

#save the dataframes
save(exportwide_simple, file = "./adjustment/workspace_adjusted_wide_simple.Rda")
save(export_simple, file = "./adjustment/workspace_adjusted_long_simple.Rda")


#visualize the studycenter effect
adjustedworkspaceall_simple <- adjustedworkspaceall_simple %>%
  mutate(studycenter = case_when(
    studycenter_1 == 1 ~ "studycenter_1",
    studycenter_2 == 1 ~ "studycenter_2",
    studycenter_3 == 1 ~ "studycenter_3",
    studycenter_4 == 1 ~ "studycenter_4",
    studycenter_5 == 1 ~ "studycenter_5",
    studycenter_6 == 1 ~ "studycenter_6",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(studycenter))

#make  plots
png("./adjustment/studycenter_adjusted.png",res = 300, width = 1600, height = 1400)
ggplot(adjustedworkspaceall_simple, aes(x = studycenter, y = adjusted_neutralization, fill = visit)) +
  geom_boxplot() +
  facet_wrap(~visit) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Adjusted Neutralization by Study Center and Visit",
       x = "Study Center",
       y = "Adjusted Neutralization")
dev.off()
png("./adjustment/studycenter_unadjusted.png",res = 300, width = 1600, height = 1400)
ggplot(adjustedworkspaceall_simple, aes(x = studycenter, y = neutralization, fill = visit)) +
  geom_boxplot() +
  facet_wrap(~visit) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Unadjusted Neutralization by Study Center and Visit",
       x = "Study Center",
       y = "Unadjusted Neutralization")
dev.off()


###studycentereffects on the group model
#visualize the studycenter effect
adjustedworkspaceall <- adjustedworkspaceall %>%
  mutate(studycenter = case_when(
    studycenter_1 == 1 ~ "studycenter_1",
    studycenter_2 == 1 ~ "studycenter_2",
    studycenter_3 == 1 ~ "studycenter_3",
    studycenter_4 == 1 ~ "studycenter_4",
    studycenter_5 == 1 ~ "studycenter_5",
    studycenter_6 == 1 ~ "studycenter_6",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(studycenter))

#make  plots
png("./adjustment/studycenter_adjusted_complex.png",res = 300, width = 1600, height = 1400)
ggplot(adjustedworkspaceall, aes(x = studycenter, y = adjusted_neutralization, fill = visit)) +
  geom_boxplot() +
  facet_wrap(~visit) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Adjusted Neutralization by Study Center and Visit",
       x = "Study Center",
       y = "Adjusted Neutralization")
dev.off()
png("./adjustment/studycenter_unadjusted_complex.png",res = 300, width = 1600, height = 1400)
ggplot(adjustedworkspaceall, aes(x = studycenter, y = neutralization, fill = visit)) +
  geom_boxplot() +
  facet_wrap(~visit) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Unadjusted Neutralization by Study Center and Visit",
       x = "Study Center",
       y = "Unadjusted Neutralization")
dev.off()
