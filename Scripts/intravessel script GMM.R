## **GMM Intravessel variation**

# Catherine Klesner
# 2025

#R packages are required:  

library(Momocs)
library(rio)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(dplyr)
library(tidyr)

## Color conventions - based on brewer palette, Dark2: 
# "Capulí" = "#1B9E77", "Piartal" = "#7570B3", "Tuza" = "#D95F02", 
# "Tuza - Red Slip" = "#E7298A", "Piartal - 1" = "#E6AB02", 
# "Piartal - 2" = "#8D62C1"


#NOTE - All data has been digitised using tspUtil and tpsDig2 

# ------------------- Single cross sections (for illustration) -----------------

single_tpsdata <- Momocs::import_tps("./Scripts/single_cross_section.TPS")
single_database <- rio::import("./Scripts/single_cross_section.csv")

## **Data Cleaning - single cross-section**  

single_database$Vessel <- as.factor(single_database$Vessel)
single_database$ware <- as.factor(single_database$ware)

single_database$Vessel <- with(single_database, reorder(Vessel, as.numeric(ware)))

# check that order is correct

table(names(single_tpsdata$coo)==single_database$Vessel)

# creation of outlines

single_shape <- Out(single_tpsdata$coo, fac = single_database)

# visualisation of out shapes

single_database$ware <- factor(single_database$ware, levels = c("Capulí", "Tuza", 
                                                                "Tuza - Red slip", "Piartal - 1", "Piartal - 2"))


png(filename = "./Figures/Figure1a.png", width = 4800, height = 3200, res=300)
panel_single <- panel(single_shape, fac = single_database$ware, 
                      cex.names = 0.6,
                      palette = pal_manual(c("#1B9E77", "#D95F02", "#E7298A",
                                             "#E6AB02", "#8D62C1")), 
                      main = "Vessels under study")  
dev.off()


# -----------------Importing and cleaning GMM Data -----------------------------


## **Importing the GMM Data into R**

tpsdata <- Momocs::import_tps("./Scripts/intravessel30.TPS")
database <- rio::import("./Scripts/intravessel30_db.csv")

## **Data Cleaning - intravessel**  


database$Vessel <- as.factor(database$Vessel)
database$ware <- as.factor(database$ware)

database$`percent_complete` <- as.numeric(database$`percent_complete`)



# ----------------------Creation of the "Out" object----------------------------  

#WARNING - you must first ensure that the *tpsdata$coo* names and the database 
#IDs are in the same order 

table(names(tpsdata$coo)==database$ID)

# Subsample the database to only include vessel cross-sections that are over 
# 99% complete

over_99 <- database[database$`percent_complete`> 0.99,]
over_99_tps <- subset(tpsdata$coo, database$`percent_complete`> 0.99,)


shape_99 <- Out(over_99_tps, fac = over_99)



# ------------------------Outline Normalisation---------------------------------

shapenorm_99 <- coo_center(shape_99)
shapenorm_99 <- coo_scale(shapenorm_99)
shapenorm_99 <- coo_close(shapenorm_99)
shapenorm2_99 <- shapenorm_99 %>% coo_slidedirection("right") %>% coo_untiltx()

# ----------------------Elliptic Fourier Analysis-------------------------------

calibrate_harmonicpower_efourier(shapenorm2_99, plot = FALSE)
calibrate_reconstructions_efourier(shapenorm2_99, range = 1:40)
calibrate_deviations_efourier(shapenorm2_99, range = c(17,20,25,30,40))

#compute efa shape with number of harmonics that account for 99.9% of variation (n=17)

efashape_99 <- efourier(shapenorm2_99, nb.h = 17, smooth.it = 0, norm = TRUE)



#visual examination of harmonics 

coo_plot(shapenorm2_99[412], col = "transparent", centroid = TRUE,)

CA230378_24 <- shapenorm2_99[412]

efashape_CA230378_17 <- efourier(CA230378_24, nb.h = 17, smooth.it = 0, norm = TRUE)
efashape_CA230378_20 <- efourier(CA230378_24, nb.h = 20, smooth.it = 0, norm = TRUE)
efashape_CA230378_25 <- efourier(CA230378_24, nb.h = 25, smooth.it = 0, norm = TRUE)
efashape_CA230378_30 <- efourier(CA230378_24, nb.h = 30, smooth.it = 0, norm = TRUE)
efashape_CA230378_40 <- efourier(CA230378_24, nb.h = 40, smooth.it = 0, norm = TRUE)

CA230378_17_xy <- efourier_i(efashape_CA230378_15, nb.h = 17, nb.pts = 200)
CA230378_20_xy <- efourier_i(efashape_CA230378_20, nb.h = 20, nb.pts = 200)
CA230378_25_xy <- efourier_i(efashape_CA230378_25, nb.h = 25, nb.pts = 200)
CA230378_30_xy <- efourier_i(efashape_CA230378_30, nb.h = 30, nb.pts = 200)
CA230378_40_xy <- efourier_i(efashape_CA230378_40, nb.h = 40, nb.pts = 200)

CA230378_24 <- as.data.frame(CA230378_24)
CA230378_17_xy <- as.data.frame(CA230378_17_xy)
CA230378_20_xy <- as.data.frame(CA230378_20_xy)
CA230378_25_xy <- as.data.frame(CA230378_25_xy)
CA230378_30_xy <- as.data.frame(CA230378_30_xy)
CA230378_40_xy <- as.data.frame(CA230378_40_xy)


harmonics <- ggplot(CA230378_24, aes(x = V1, y = V2)) + 
  geom_path(color = "black") + 
  theme_minimal() +
  labs(title = "CA230378_24", x = "x", y = "y") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(-1.4, 1.4)) + 
  scale_y_continuous(limits = c(-1.4, 1.4))

# Add the Fourier shapes with different colors and legend
harmonics + 
  geom_path(data = CA230378_17_xy, aes(x = x, y = y, color = "17 harmonics")) +
  geom_path(data = CA230378_20_xy, aes(x = x, y = y, color = "20 harmonics")) +
  geom_path(data = CA230378_25_xy, aes(x = x, y = y, color = "25 harmonics")) +
  geom_path(data = CA230378_30_xy, aes(x = x, y = y, color = "30 harmonics")) +
  geom_path(data = CA230378_40_xy, aes(x = x, y = y, color = "40 harmonics")) +
  scale_color_manual(values = c("17 harmonics" = "red", 
                                "20 harmonics" = "orange", 
                                "25 harmonics" = "yellow", 
                                "30 harmonics" = "green", 
                                "40 harmonics" = "blue")) +
  labs(color = "Number of Harmonics")

png(filename = "./Figures/Figure4.png", width = 2400, height = 2400, res=300)
plot(harmonics + 
       geom_path(data = CA230378_17_xy, aes(x = x, y = y, color = "17 harmonics")) +
       geom_path(data = CA230378_20_xy, aes(x = x, y = y, color = "20 harmonics")) +
       geom_path(data = CA230378_25_xy, aes(x = x, y = y, color = "25 harmonics")) +
       geom_path(data = CA230378_30_xy, aes(x = x, y = y, color = "30 harmonics")) +
       geom_path(data = CA230378_40_xy, aes(x = x, y = y, color = "40 harmonics")) +
       scale_color_manual(values = c("17 harmonics" = "red", 
                                     "20 harmonics" = "orange", 
                                     "25 harmonics" = "yellow", 
                                     "30 harmonics" = "green", 
                                     "40 harmonics" = "blue")) +
       labs(color = "Number of Harmonics"))
dev.off()



#We will also run the analysis considering 30 harmonics.

efashape_99_30 <- efourier(shapenorm2_99, nb.h = 30, smooth.it = 0, norm = TRUE)

#figures

png(filename = "./Figures/Figure3a.png", width = 1600, height = 1600, res=300)
calibrate_harmonicpower_efourier(shapenorm2_99, nb.h = 9, plot = TRUE)
dev.off()

set.seed(121)

png(filename = "./Figures/Figure3b.png", width = 2400, height = 1600, res=300)
calibrate_reconstructions_efourier(shapenorm2_99, 408, range = 1:30)
dev.off()

# ---------------------Principal Component Analysis-----------------------------  

pcashape_99 <- PCA(efashape_99)
pcashape_99_30 <- PCA(efashape_99_30)

scree(pcashape_99)

scree(pcashape_99_30)


#plot our vessels within a morphospace representative of these PCs:  


pdf("./Figures/Figure5_1.pdf", width = 11, height = 7)
plot_PCA(pcashape_99_30, axes = c(1,2), over_99$ware, 
         morphospace_position = "range", zoom = 0.9, chull = FALSE, center_origin = FALSE,
         palette = pal_manual(c("#1B9E77", "#E6AB02", "#8D62C1","#D95F02", "#E7298A")), 
         title = "PC 1 vs PC 2 considering cross-sections >99% preserved") %>% 
  layer_points(cex = 0.5)
dev.off()


pdf("./Figures/Figure5_2.pdf", width = 11, height = 7)
plot_PCA(pcashape_99_30, axes = c(1,2), over_99$Vessel, 
         morphospace_position = "range_axes", zoom = 0.9, chull = FALSE, 
         points = FALSE,
         palette = pal_manual(c("#000000")), 
         center_origin = FALSE,
         title = "PC 1 vs PC 2 considering cross-sections >99% preserved")%>% layer_ellipses(conf = 0.9)
dev.off()


#examine differences in PC contributions by vessel or ware

boxplot_99 <- boxplot(pcashape_99_30, over_99$ware, nax = 1:7)
boxplot_99 + scale_fill_brewer(palette = "Dark2") + 
  ggtitle("PC contribution sections which are > 99% complete")


#Export the EFA as txt

file.create('./Scripts/over_99_EFA_30.txt')
Momocs::export(efashape_99_30, './Scripts/over_99_EFA_30.txt')
 
# -------------------------------Assess EFA variance ---------------------------


EFA <- import('./Scripts/over_99_EFA_30.txt')

EFA$Vessel <- as.factor(EFA$Vessel)
EFA$ware <- as.factor(EFA$ware)
EFA <- as.data.frame(EFA)

EFA_dist <- EFA %>%
  select(-name, -ID, -Vessel, -ware, -collection, -percent_complete)

# calculated multivariate homogeneity of groups dispersions (variances)

d_EFA <- dist(EFA_dist, method = "euclidean")
dist_EFA_c <- betadisper(d_EFA, EFA$Vessel, type = c("centroid"))


anova(dist_EFA_c)
plot(dist_EFA_c, ellipse = TRUE, hull = FALSE, conf = 0.90, label = FALSE) 
boxplot(dist_EFA_c, aes(x = Vessel, y = EFA$Vessel), fill = (EFA$ware))
permutest_EFA <- permutest(dist_EFA_c, pairwise = TRUE, permutations = 99)

permutest_EFA

dist_EFA_c.HSD <- TukeyHSD(dist_EFA_c)
plot(dist_EFA_c.HSD)

distances_EFA <- dist_EFA_c$distances


# Ensure EFA has the same number of rows as the distances vector
if(nrow(EFA) == length(distances_EFA)) {
  # Add the distances as a new column to the PC dataframe
  EFA$Distance <- distances_EFA
  
  # View the updated dataframe
  print(EFA)
} else {
  cat("Error: Mismatch in the number of rows in PC and the length of distances.\n")
}

EFA$Vessel <- with(EFA, reorder(Vessel, as.numeric(ware)))

write.csv(EFA, "./Scripts/EFA_output.csv", row.names = FALSE)


# Plotting distances

boxplot_EFA_30 <- ggplot(EFA, aes(x = Vessel, y = Distance, fill = ware)) +
  geom_boxplot(
    width = .6, 
    outlier.shape = NA
  ) +
  coord_cartesian(xlim = c(0, 30), clip = "off") + 
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, color = "black", fill = "white", stroke = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_fill_manual(values = c(
    "Capulí" = "#1B9E77", 
    "Tuza" = "#D95F02", 
    "Tuza - Red slip" = "#E7298A", 
    "Piartal - 1" = "#E6AB02", 
    "Piartal - 2" = "#8D62C1"
  ))


png(filename = "./Figures/Figure6.png", width = 2400, height = 1600, res=300)
plot(boxplot_EFA_30)
dev.off()



#---------------------Individual vessel morphological variation ----------------

# First we will subsample the database to only include cross-sections from specific vessels 
# Selected the following vessels based on their average group distance to the centroid
# values: CA230335 (lowest value), CA230365 (highest value), CA230377 (median Capuli value), 
# CA230417 (median Piartal - 1 value), CA230347 (median Piartal - 2 value), 
# and CA230715 (median Tuza value)

CA230335 <- database[database$`Vessel`== 'CA230335',]
CA230365 <- database[database$`Vessel`== 'CA230365',]
CA230377 <- database[database$`Vessel`== 'CA230377',]
CA230417 <- database[database$`Vessel`== 'CA230417',]
CA230347 <- database[database$`Vessel`== 'CA230347',]
CA230715 <- database[database$`Vessel`== 'CA230715',]

CA230335_tps <- subset(tpsdata$coo, database$`Vessel`== 'CA230335',)
CA230365_tps <- subset(tpsdata$coo, database$`Vessel`== 'CA230365',)
CA230377_tps <- subset(tpsdata$coo, database$`Vessel`== 'CA230377',)
CA230417_tps <- subset(tpsdata$coo, database$`Vessel`== 'CA230417',)
CA230347_tps <- subset(tpsdata$coo, database$`Vessel`== 'CA230347',)
CA230715_tps <- subset(tpsdata$coo, database$`Vessel`== 'CA230715',)

shape_CA230335 <- Out(CA230335_tps, fac = CA230335)
shape_CA230365 <- Out(CA230365_tps, fac = CA230365)
shape_CA230377 <- Out(CA230377_tps, fac = CA230377)
shape_CA230417 <- Out(CA230417_tps, fac = CA230417)
shape_CA230347 <- Out(CA230347_tps, fac = CA230347)
shape_CA230715 <- Out(CA230715_tps, fac = CA230715)


## **normalisation**

shapenorm_CA230335 <- coo_center(shape_CA230335)
shapenorm_CA230335 <- coo_scale(shapenorm_CA230335)
shapenorm_CA230335 <- coo_close(shapenorm_CA230335)
shapenorm_CA230335_2 <- shapenorm_CA230335 %>% coo_slidedirection("right") %>% 
  coo_untiltx()

shapenorm_CA230365 <- coo_center(shape_CA230365)
shapenorm_CA230365 <- coo_scale(shapenorm_CA230365)
shapenorm_CA230365 <- coo_close(shapenorm_CA230365)
shapenorm_CA230365_2 <- shapenorm_CA230365 %>% coo_slidedirection("right") %>% 
  coo_untiltx()

shapenorm_CA230377 <- coo_center(shape_CA230377)
shapenorm_CA230377 <- coo_scale(shapenorm_CA230377)
shapenorm_CA230377 <- coo_close(shapenorm_CA230377)
shapenorm_CA230377_2 <- shapenorm_CA230377 %>% coo_slidedirection("right") %>% 
  coo_untiltx()

shapenorm_CA230417 <- coo_center(shape_CA230417)
shapenorm_CA230417 <- coo_scale(shapenorm_CA230417)
shapenorm_CA230417 <- coo_close(shapenorm_CA230417)
shapenorm_CA230417_2 <- shapenorm_CA230417 %>% coo_slidedirection("right") %>% 
  coo_untiltx()

shapenorm_CA230347 <- coo_center(shape_CA230347)
shapenorm_CA230347 <- coo_scale(shapenorm_CA230347)
shapenorm_CA230347 <- coo_close(shapenorm_CA230347)
shapenorm_CA230347_2 <- shapenorm_CA230347 %>% coo_slidedirection("right") %>% 
  coo_untiltx()

shapenorm_CA230715 <- coo_center(shape_CA230715)
shapenorm_CA230715 <- coo_scale(shapenorm_CA230715)
shapenorm_CA230715 <- coo_close(shapenorm_CA230715)
shapenorm_CA230715_2 <- shapenorm_CA230715 %>% coo_slidedirection("right") %>% 
  coo_untiltx()

# EFA

efashape_CA230335_30 <- efourier(shapenorm_CA230335_2, nb.h = 30, smooth.it = 0, norm = TRUE)
efashape_CA230365_30 <- efourier(shapenorm_CA230365_2, nb.h = 30, smooth.it = 0, norm = TRUE)
efashape_CA230377_30 <- efourier(shapenorm_CA230377_2, nb.h = 30, smooth.it = 0, norm = TRUE)
efashape_CA230417_30 <- efourier(shapenorm_CA230417_2, nb.h = 30, smooth.it = 0, norm = TRUE)
efashape_CA230347_30 <- efourier(shapenorm_CA230347_2, nb.h = 30, smooth.it = 0, norm = TRUE)
efashape_CA230715_30 <- efourier(shapenorm_CA230715_2, nb.h = 30, smooth.it = 0, norm = TRUE)

#PCA

pcashape_CA230335 <- PCA(efashape_CA230335_30)
pcashape_CA230365 <- PCA(efashape_CA230365_30)
pcashape_CA230377 <- PCA(efashape_CA230377_30)
pcashape_CA230417 <- PCA(efashape_CA230417_30)
pcashape_CA230347 <- PCA(efashape_CA230347_30)
pcashape_CA230715 <- PCA(efashape_CA230715_30)

scree(pcashape_CA230335)
scree(pcashape_CA230365)
scree(pcashape_CA230377)
scree(pcashape_CA230417)
scree(pcashape_CA230347)
scree(pcashape_CA230715)

png(filename = "./Figures/Figure7a.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape_CA230365, axes = c(1,2), CA230365$Vessel, palette = pal_manual(c("#1B9E77")), 
         morphospace_position = "range", zoom = 0.8, chull = FALSE, center_origin = FALSE,
         title = "CA230365") %>% 
  layer_points(cex = 1)
dev.off()

png(filename = "./Figures/Figure7b.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape_CA230377, axes = c(1,2), CA230377$Vessel, palette = pal_manual(c("#1B9E77")), 
         morphospace_position = "range", zoom = 0.8, chull = FALSE, center_origin = FALSE,
         title = "CA230377") %>% 
  layer_points(cex = 1)
dev.off()

png(filename = "./Figures/Figure7c.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape_CA230417, axes = c(1,2), CA230417$Vessel, palette = pal_manual(c("#E6AB02")), 
         morphospace_position = "range", zoom = 0.8, chull = FALSE, center_origin = FALSE,
         title = "CA230417") %>% 
  layer_points(cex = 1)
dev.off()

png(filename = "./Figures/Figure7d.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape_CA230347, axes = c(1,2), CA230347$Vessel, palette = pal_manual(c("#8D62C1")), 
         morphospace_position = "range", zoom = 0.8, chull = FALSE, center_origin = FALSE,
         title = "CA230347") %>% 
  layer_points(cex = 1)
dev.off()

png(filename = "./Figures/Figure7e.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape_CA230715, axes = c(1,2), CA230715$Vessel, palette = pal_manual(c("#D95f02")), 
         morphospace_position = "range", zoom = 0.8, chull = FALSE, center_origin = FALSE,
         title = "CA230715") %>% 
  layer_points(cex = 1)
dev.off()

png(filename = "./Figures/Figure7f.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape_CA230335, axes = c(1,2), CA230335$Vessel, palette = pal_manual(c("#D95F02")), 
         morphospace_position = "range", zoom = 0.8, chull = FALSE, center_origin = FALSE,
         title = "CA230335") %>% 
  layer_points(cex = 1)
dev.off()

# ------------------------------Acknowledgement---------------------------------

# The GMM workflow was adapted from a workshop designed by Chris Hoggard 
# many thanks to him for providing it online!

citation()
citation("Momocs")
citation("vegan")
citation("ggplot2")
citation("dplyr")
citation("tidyr")




