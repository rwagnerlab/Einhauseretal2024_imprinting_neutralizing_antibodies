##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
library(tidyverse)
library(ggplot2)
library(mgcv)

library(tidymv)
library(gridExtra)
library(gratia)

library(readxl)
library(writexl)
library(coin)

rm(list = ls())

#load data
load("./adjustment/workspace_adjusted_long_simple.Rda")
vaxtime <- read_xlsx("./vaxtime.xlsx")

# Perform the left join to add the Timevax column to export_simple
merged_df <- export_simple %>%
  left_join(vaxtime %>% select(ID, Visit, Timevax), by = c("ID_study" = "ID", "visit" = "Visit"))
merged_df$Timevax <- as.numeric(merged_df$Timevax)
btis<- subset(merged_df, merged_df$vaccine == 1)
btis$group <- droplevels(btis$group)
btis$virus <- as.factor(btis$virus)
btis<-btis[complete.cases(btis),]

#keep only the homologous neutralization
btis <- btis %>%
  filter((group == "FO" & virus == "BA1") |
           (group == "FD" & virus == "Delta") |
           (group == "FA" & virus == "Alpha"))

btis$group<- as.character(btis$group)
btis$group[btis$group == "FO"] <- "Vaccinated + Omicron BTI"
btis$group[btis$group == "FA"] <- "Vaccinated + Alpha BTI"
btis$group[btis$group == "FD"] <- "Vaccinated + Delta BTI"
btis$group <- as.factor(btis$group)
#make models

btis_v1 <- btis[btis$visit == "v1",]
btis_v1$virus <- droplevels(btis_v1$virus)
modelhomv1<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax) + s(Timevax, by = virus), data = btis_v1)
btis_v2 <- btis[btis$visit == "v2",]
btis_v2$virus <- droplevels(btis_v2$virus)
modelhomv2<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax)+ s(Timevax, by = virus), data = btis_v2)
btis_v4 <- btis[btis$visit == "v4",]
btis_v4$virus <- droplevels(btis_v4$virus)
modelhomv4<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax)+ s(Timevax, by = virus), data = btis_v4)
btis_v5 <- btis[btis$visit == "v5",]
btis_v5$virus <- droplevels(btis_v5$virus)
modelhomv5<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax)+ s(Timevax, by = virus), data = btis_v5)

#check summaries
summary(modelhomv1)
summary(modelhomv2)
summary(modelhomv4)
summary(modelhomv5)

#create plots
phom1 <- plot_smooths(model = modelhomv1, series = Timevax, comparison = group, virus)+theme(legend.position = "bottom")+ggtitle("Visit 1")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))
phom2 <- plot_smooths(model = modelhomv2, series = Timevax, comparison = group, virus)+theme(legend.position = "bottom")+ggtitle("Visit 2")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))
phom4 <- plot_smooths(model = modelhomv4, series = Timevax, comparison = group, virus)+theme(legend.position = "bottom")+ggtitle("Visit 4")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))
phom5 <- plot_smooths(model = modelhomv5, series = Timevax, comparison = group, virus)+theme(legend.position = "bottom")+ggtitle("Visit 5")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))
if(!dir.exists("./supplfig4_gams")){dir.create("./supplfig4_gams")}
png("./supplfig4_gams/gamhom.png", res = 300, height = 2600, width = 4000)
grid.arrange(phom1,phom2,phom4,phom5)
dev.off()


modelsimplev1<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax), data = btis_v1)
ps1<- plot_smooths(model = modelsimplev1, series = Timevax, comparison = group)+theme(legend.position = "bottom")+ggtitle("Visit 1")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))

modelsimplev2<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax), data = btis_v2)
ps2<- plot_smooths(model = modelsimplev2, series = Timevax, comparison = group)+theme(legend.position = "bottom")+ggtitle("Visit 2")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))

modelsimplev4<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax), data = btis_v4)
ps4<- plot_smooths(model = modelsimplev4, series = Timevax, comparison = group)+theme(legend.position = "bottom")+ggtitle("Visit 4")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))

modelsimplev5<-gam(adjusted_neutralization ~ s(Timevax, by = group)+ s(Timevax), data = btis_v5)
ps5<- plot_smooths(model = modelsimplev5, series = Timevax, comparison = group)+theme(legend.position = "bottom")+ggtitle("Visit 5")+ylab("IC50 vs. breakthrough variant")+xlab("Days from last vaccine till BTI")#+ylim(c(0,3500))

result_1 <- cor.test(btis_v1$adjusted_neutralization, btis_v1$Timevax, method = "spearman")
result_2 <- cor.test(btis_v2$adjusted_neutralization, btis_v2$Timevax, method = "spearman")
result_4 <- cor.test(btis_v4$adjusted_neutralization, btis_v4$Timevax, method = "spearman")
result_5 <- cor.test(btis_v5$adjusted_neutralization, btis_v5$Timevax, method = "spearman")

sink("./supplfig4_gams/gamssummary.txt")
print("Visit 1-----------------------------------------------------------")
print(result_1)
summary(modelsimplev1)
print("Visit 2-----------------------------------------------------------")
print(result_2)
summary(modelsimplev2)
print("Visit 4-----------------------------------------------------------")
print(result_4)
summary(modelsimplev4)
print("Visit 5-----------------------------------------------------------")
print(result_5)
summary(modelsimplev5)
sink()
png("./supplfig4_gams/gamhom_simple.png", res = 300, height = 2600, width = 4000)
grid.arrange(ps1,ps2,ps4,ps5)
dev.off()
