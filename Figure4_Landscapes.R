#this script uses functions and code from "BA.2 and BA.5 omicron differ immunologically from both BA.1 omicron and pre-omicron variants" by Roessler, A., Netzl, A., et al. 2022
#available on repository: https://github.com/acorg/roessler_netzl_et_al2022/
#DOI: 10.5281/zenodo.7341691..
##adapted by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(ablandscapes)
library(titertools)

library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)
library(gdata)
getwd()
set.seed(100)
#help(package="ablandscapes")
#help(package="meantiter")
#help(package="Racmacs")
source("./functions/remove_reactivity_bias.R")
source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")
'%notin%' <- Negate('%in%')
figure_dir <- "./landscapes/figures"
if(!dir.exists("./landscapes/singles")){dir.create("./landscapes/singles")}
if(!dir.exists(figure_dir)){dir.create(figure_dir)}

sr_colors <- read.csv("./sr_group_colors.csv", sep = ";", row.names = "SerumGroup")


# read in data 
load("./adjustment/workspace_adjusted_long_simple.Rda")

visits <- unique(export_simple$visit)

#loop over the visits
for(i in 1: length(visits)){
  #for debugging: i = 1
  map_long <- subset(export_simple, export_simple$visit == visits[i])  
  # Read the base map
  map <- read.acmap(paste0("./fig3_mapping/map",i,".ace"))
  map <- removeAntigens(map, c("XBB","JN","BQ"))
  
  lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
  #ags_to_fit_lndscp <- agNames(map)
  ags_to_fit_lndscp <- c("D614G", "Alpha", "Delta","BA1", "BA2", "BA5")
  ag.colors <- c('#C0392B', '#9B59B6', '#2980B9',
                 '#E67E22','#F1C40F', '#1ABC9C', '#34495e','#5f495e','#7f495e')
  agFill(map) <- ag.colors[1:length(ags_to_fit_lndscp)] 
  
  map_long %>%
    ungroup() %>%
    select(adjusted_neutralization, virus, ID_study, group) -> titerdata
  
  colnames(titerdata) <-c("titer", "ag_name", "sr_name", "sr_group")
  
  # remove non titrated sera and variants not measured for all visists
  titerdata <- titerdata %>%
    filter(titer != "*")%>%
    filter(ag_name %in% ags_to_fit_lndscp)
  
  
  titerdata %>%
    group_by(
      sr_group
    ) -> titerdata
  
  titerdata %>%
    group_map(
      get_titertable
    ) -> titertables
  
  lndscp_fits <- lapply(
    titertables,
    function(titertable) {
      
      ablandscape.fit(
        titers = titertable[,ags_to_fit_lndscp],
        bandwidth = 1,
        degree = 1,
        method = "cone",
        error.sd = 1,
        acmap = map,
        control = list(
          optimise.cone.slope = TRUE
        )
      )
      
    }
  )
  
  titertables_groups <- group_data(titerdata)
  
  # Add impulses
  titerdata %>%
    group_by(
      sr_group,
      ag_name
    ) %>%
    summarize(gmt = titertools::gmt(titer, dilution_stepsize = 2)["mean", "estimate"]) %>%
    # manually set GMT's that are lower than that to LOD2
    mutate(gmt = ifelse(gmt < log2(0.8), log2(0.8), gmt))-> gmt_data
  
  # angle for html page
  if(i == 1){
    angle <- list(
      rotation = c(-1.3937, -0.0056, -0.3038), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
      translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
      zoom = 2.5    # zoom = 1.1646 # higher is more zoomed out
   )
  }else{
    angle <- list(
      rotation = c(-1.3937, -0.0056, -0.3038), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
      translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
      zoom = 2.0    # zoom = 1.1646 # higher is more zoomed out
    )
  }

  #define
  titertables_groups$sr_group
  va_indices <- sapply(c("FA"),function(x){
    match(x, titertables_groups$sr_group)
  })
  vd_indices <- sapply(c("FD"),function(x){
    match(x, titertables_groups$sr_group)
  })
  vo_indices <- sapply(c("FO"),function(x){
    match(x, titertables_groups$sr_group)
  })
  ua_indices <- sapply(c("UA"),function(x){
    match(x, titertables_groups$sr_group)
  })
  ud_indices <- sapply(c("UD"),function(x){
    match(x, titertables_groups$sr_group)
  })
  uo_indices <- sapply(c("UO"),function(x){
    match(x, titertables_groups$sr_group)
  })
  indices <- list("FA" = va_indices, "FD" = vd_indices,
                  "FO" = vo_indices, "UA" = ua_indices, "UD" = ud_indices,"UO" = uo_indices) 
  
  lndscp_list <- list()
  data3js <- base_plot_data3js(map, lndscp_fits, agNames(map), lims, agNames(map))
  
  #optional plot single landscapes
  
  #for(ind in 1:length(indices)){
  #  lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups[indices[[ind]],], lndscp_fits[indices[[ind]]], map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors,
  #                                          hide_buttons = FALSE)
  #  
  #  lndscp <-r3js(
  #    lndscp_3js,
  #    rotation = angle$rotation,
  #    zoom = angle$zoom
  #  )
  #  
  #  lndscp<-light3js(
  #    lndscp,
  #    position = NULL,
  #    intensity = 0.5,
  #    type = "ambient",
  #    col = "white"
  #  )
  #  lndscp_list[[names(indices)[ind]]] <- lndscp
  #  saveWidget(lndscp, file = print(paste0("./landscapes/singles/", names(indices)[ind],"visit",visits[i],".html")), title = print(paste0(names(indices)[ind],"-Visit",visits[i])), selfcontained = TRUE )
  #  webshot(url = print(paste0("./landscapes/singles/", names(indices)[ind],"visit",visits[i],".html")),file = print(paste0("./landscapes/singles/",names(indices)[ind],"-Visit",visits[i],".png")),vwidth = 4000, vheight = 2000, zoom = 2)
  #}
  
  #plot landscape with everything
  lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups, lndscp_fits, map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors,
                                          hide_buttons = FALSE) 
  #format legend
  #leg <- as.factor(titertables_groups$sr_group)
  #lndscp_3js <- legend3js(
  # data3js = lndscp_3js,
  # legend = levels(leg),
  # fill = sr_colors$Color)
  
  #set lighting
  lndscp_3js<-light3js(
    lndscp_3js,
    position = c(5,5,5),
    intensity = 0.35,
    type = "ambient",
    col = "white"
  )
  lndscp_3js<-light3js(
    lndscp_3js,
    position = c(5,1000,5),
    intensity = 0.45,
    type = "directional",
    col = "white"
  )
  lndscp <-r3js(
    lndscp_3js,
    rotation = angle$rotation,
    zoom = angle$zoom
  )
  
  
  lndscp_list[["all"]] <- lndscp
  saveWidget(lndscp,paste0("./landscapes/figures/visit",visits[i],".html"), title = paste0("Visit",visits[i]), selfcontained = TRUE )
  #save3js(data3js = lndscp_3js,paste0("./landscapes/figures/visit",visits[i],".html"), title = paste0("Visit",visits[i]))
  saveRDS(lndscp_list, paste0("landscapes/all",visits[i],"_landscapes.rds"))
  webshot(url = paste0("./landscapes/figures/visit",visits[i],".html"),file = paste0("./landscapes/figures/visit",visits[i],".png"),vwidth = 4000, vheight = 2000, zoom = 2)
}

