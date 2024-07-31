##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
library(Racmacs)
library(png)

rm(list = ls())
options(RacOptimizer.num_cores = parallel::detectCores())


# List all files in the directory
path = c("./fig3_mapping/")
files <- list.files(path)

#first four maps are the full maps
maps<-vector("list", length = 4)
for(i in 1:4){
  maps[[i]] <- read.acmap(paste0(path,"map",i,".ace"))
}

bs_maps<-vector("list", length = 4)
for(i in 1:4){
  map <- maps[[i]]
  map<-keepBestOptimization(map)
  if(i<4){
    map<-removeAntigens(map,c("XBB","JN"))
  }
  bmap <-bootstrapMap(
    map = map,
    method = "resample",
    bootstrap_repeats = 1000,
    bootstrap_ags = TRUE,
    bootstrap_sr = FALSE,
    reoptimize = TRUE,
    optimizations_per_repeat = 100
  )
  bs_maps[[i]]<-bmap
  rm(map)
  rm(bmap)
}

blobmaps <- vector("list", length = 4)
for(i in 1:4){
  map<-bs_maps[[i]]
  map<-bootstrapBlobs(
    map,
    conf.level = 0.68,
    smoothing = 6,
    gridspacing = 0.25,
    antigens = TRUE,
    sera = FALSE,
    method = "ks"
  )
  blobmaps[[i]]<-map
  rm(map)
}

map1 <- blobmaps[[1]]
map2 <- blobmaps[[2]]
map3 <- blobmaps[[3]]
map4 <- blobmaps[[4]]
visits<- c("1","2","4","5")
for(i in 1:4){
  png(file =paste0("./fig3_mapping/antigenicmap_bootstrap",visits[i],".png"), res = 400, width=3000, height=3000)
    map <- blobmaps[[i]]
    srSize(map) <- 5
    
    srOpacity(map) <- 0.2
    plot(map, plot_ags = TRUE, plot_labels = "antigens", plot_blobs = TRUE, plot_sr = FALSE)
    title(paste0("Visit ",visits[i]," Bootstrap"))
  dev.off()
}

plot(map2, plot_ags = TRUE, plot_labels = "antigens", plot_blobs = TRUE, plot_sr = FALSE)
title("Visit 2 Bootstrap")
plot(map3, plot_ags = TRUE, plot_labels = "antigens", plot_blobs = TRUE, plot_sr = FALSE)
title("Visit 3 Bootstrap")
plot(map4, plot_ags = TRUE, plot_labels = "antigens", plot_blobs = TRUE, plot_sr = FALSE)
title("Visit 4 Bootstrap")

for(i in 1:4){
  save.acmap(blobmaps[[i]],filename = paste0(path,"blobmap",i,".ace"))
}

