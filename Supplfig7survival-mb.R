##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
rm(list=ls())
library(dplyr)
library(writexl)
library(tidyverse)
library(gridExtra)
library("survival")
library("survminer")
library("Rcpp")
library(reshape2)
`%notin%` <- Negate(`%in%`)
#load data and define necessary variables
load("./adjustment/workspace_adjusted_long_simple.Rda")
workspace_long<- export_simple
rm(export_simple)

workspace_long$events <- ifelse(workspace_long$adjusted_neutralization <= 2560, 1, 0)
workspace_long<-na.omit(workspace_long)
workspace_nobq<- subset(workspace_long, virus %notin% c("BQ","XBB","JN"))

visits<- unique(workspace_nobq$visit)
# split the data into the visits
visit_data <- vector("list", length(visits))
for( i in 1:length(visits)){
  visit_data[[i]]<-subset(workspace_nobq,workspace_nobq$visit == visits[i])
}
plots<- vector("list",length(visits))
diff<-vector("list",length(visits))
pairwisediff<-vector("list",length(visits))
for( i in 1:length(visits)){
  single_visit <- visit_data[[i]]
  fit <- survfit(Surv(adjusted_neutralization, events) ~ group, data = single_visit)
  print(fit)               
  
  if ( i == 1){
  ggsurv<-ggsurvplot(fit,
             pval = FALSE, conf.int = TRUE,
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             linetype = 1, # Change line type by groups
             surv.median.line = "hv", # Specify median survival
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette = c("#9127e3", "#2E61FF","#bd6f04","#f127e3", "#2Eb1FF","#f0a748"))
  }else{
    ggsurv<-ggsurvplot(fit,
                       pval = FALSE, conf.int = TRUE,
                       legend = "none",
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       linetype = 1, # Change line type by groups
                       surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_bw(), # Change ggplot2 theme
                       palette = c("#9127e3", "#2E61FF","#bd6f04","#f127e3", "#2Eb1FF","#f0a748"))
  }
  
  # Customize axis labels
  ggsurv$plot <- ggsurv$plot + labs(title = paste0("Visit ", visits[i] ) ,x = "Neutralization Magnitude", y = "Neutralization Breadth", color = "Group", fill = "Group")+
    scale_color_manual(values = c("#9127e3", "#2E61FF", "#bd6f04", "#f127e3", "#2Eb1FF", "#f0a748"),
                       labels = c("Vaccinated Alpha", "Vaccinated Delta", "Vaccinated Omicron", "Unvaccinated Alpha", "Unvaccinated Delta", "Unvaccinated Omicron")) +
    scale_fill_manual(values = c("#9127e3", "#2E61FF", "#bd6f04", "#f127e3", "#2Eb1FF", "#f0a748"), 
                      labels = c("Vaccinated Alpha", "Vaccinated Delta", "Vaccinated Omicron", "Unvaccinated Alpha", "Unvaccinated Delta", "Unvaccinated Omicron"))
  
  plots[[i]] <- ggsurv$plot
  
  surv_diff <- survdiff(Surv(adjusted_neutralization, events) ~ group, data = single_visit)
  diff[[i]]<-surv_diff
  pairwise_diff<-pairwise_survdiff(Surv(adjusted_neutralization, events) ~ group, data = single_visit)
  pvals<-pairwise_diff$p.value
  write_xlsx(pvals,path = paste0("./fig2_mb/pairwise_pval",visits[i],".xlsx"))
  pairwisediff[[i]]<-pairwise_diff
}



for ( i in 1:length(visits)){
  #i=1
pvalue_table <- pairwisediff[[i]]
pvalue_table<- as.data.frame(pvalue_table$p.value)
pvalue_table <-pvalue_table %>%
  rownames_to_column() %>%
  gather(colname, pvalue, -rowname)
# Plot heatmap of p-values
pvalue_heatmap <- ggplot(pvalue_table, aes(x = rowname, y = colname, fill = pvalue)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red",      # Color for low values
    mid = "white",    # Color for the midpoint
    high = "white",   # Color for high values
    midpoint = 0.05, # Midpoint of the gradient
    limits = c(0, 0.05), # limits for consitent scaling
    na.value = "grey70"  # Color for NA values
  )+
  geom_text(aes(label = round(pvalue, 4)), size = 4) +
  theme_minimal() +
  labs(title = paste("Log-rank p-value Heatmap for Visit", visits[i]), x = "", y = "")

plots[[i + length(visits)]] <- pvalue_heatmap
}
# Arrange the plots in a 2x2 grid
png(filename = "./fig2_mb/logrank_pval.png",res = 600, width = 16000, height = 6000)
grid.arrange(grobs = plots, ncol = 4, nrow = 2)
dev.off()
png(filename = "./fig2_mb/logrank_pval_alternative.png",res = 600, width = 8000, height = 10000)
grid.arrange(grobs = plots, ncol = 2, nrow = 4)
dev.off()
