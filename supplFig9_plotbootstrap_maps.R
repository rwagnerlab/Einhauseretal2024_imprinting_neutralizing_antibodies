##by Sebastian Einhauser for Einhauser et al. 
#Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
library(Racmacs)
library(png)

rm(list = ls())


m1<-read.acmap("./fig3_mapping/blobmap1.ace")
m2<-read.acmap("./fig3_mapping/blobmap2.ace")
m3<-read.acmap("./fig3_mapping/blobmap3.ace")
m4<-read.acmap("./fig3_mapping/blobmap4.ace")




png(file="./fig3_mapping/antigenicmapfinal.png", res = 400, width=16000, height=9000)
par(mfrow = c(3,4))
plot(m1, plot_ags = TRUE,plot_blobs = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 1 All", cex.main = 3)
plot(m2, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 2 All", cex.main = 3)
plot(m3, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
title( main ="Visit 4 All", cex.main = 3)
plot(m4, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
title( main ="Visit 5 All", cex.main = 3)