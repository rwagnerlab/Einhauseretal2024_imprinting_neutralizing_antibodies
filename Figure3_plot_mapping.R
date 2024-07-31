##by Sebastian Einhauser for Einhauser et al. 
##Longitudinal effects of SARS-CoV-2 breakthrough infection on imprinting of neutralizing antibody responses 
#plots antigenic maps calculateed previously
library(Racmacs)
library(png)

rm(list = ls())


m1<-read.acmap("./fig3_mapping/map1.ace")
m2<-read.acmap("./fig3_mapping/map2.ace")
m3<-read.acmap("./fig3_mapping/map3.ace")
m4<-read.acmap("./fig3_mapping/map4.ace")

m5<-read.acmap("./fig3_mapping/map5.ace")
m6<-read.acmap("./fig3_mapping/map6.ace")
m7<-read.acmap("./fig3_mapping/map7.ace")
m8<-read.acmap("./fig3_mapping/map8.ace")

m9<-read.acmap("./fig3_mapping/map9.ace")
m10<-read.acmap("./fig3_mapping/map10.ace")
m11<-read.acmap("./fig3_mapping/map11.ace")
m12<-read.acmap("./fig3_mapping/map12.ace")


png(file="./fig3_mapping/antigenicmapfinal.png", res = 400, width=16000, height=9000)
par(mfrow = c(3,4))
plot(m1, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 1 All", cex.main = 3)
plot(m2, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 2 All", cex.main = 3)
plot(m3, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 4 All", cex.main = 3)
plot(m4, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 5 All", cex.main = 3)

plot(m5, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 1 Vaccinated", cex.main = 3)
plot(m6, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 2 Vaccinated", cex.main = 3)
plot(m7, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 4 Vaccinated", cex.main = 3)
plot(m8, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 5 Vaccinated", cex.main = 3)

plot(m9, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 1 Unvaccinated", cex.main = 3)
plot(m10, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 2 Unvaccinated", cex.main = 3)
plot(m11, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 4 Unvaccinated", cex.main = 3)
plot(m12, plot_ags = TRUE, plot_sr = TRUE, plot_labels = 'antigens', indicate_outliers = "arrowheads", xlim = c(-3,3), ylim = c(-3,3))
title( main ="Visit 5 Unvaccinated", cex.main = 3)

dev.off()






