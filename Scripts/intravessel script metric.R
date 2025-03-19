## **Metric Intravessel variation**

# Catherine Klesner and Jasmine Vieri
# 2025

#R packages are required:  

library(rio)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(ggdist)
library(gghalves)
library(dplyr)
library(tidyr)
library(car)

## Color conventions - based on brewer palette, Dark2: 
# "Capulí" = "#1B9E77", "Piartal" = "#7570B3", "Tuza" = "#D95F02", 
# "Tuza - Red Slip" = "#E7298A", "Piartal - 1" = "#E6AB02", 
# "Piartal - 2" = "#8D62C1"


# -----------------------------Wall thickness calcs-----------------------------

# Assessing metric measurements - horizontal wall thickness boxplot 
# Need to select 200 measurements

horizontal_wall_thickness <- read.csv("./Scripts/horizontal_wall_thickness.csv",
                                      na.strings="",header = TRUE)


horizontal_lists <- list()
horizontal_selected <- list()

for (i in 1:ncol(horizontal_wall_thickness)) {
  
  horizontal_lists[[i]] <- horizontal_wall_thickness[, i]
  
  #exclude NA observations
  horizontal_lists[[i]] <-
    horizontal_lists[[i]][!is.na(horizontal_lists[[i]])]
  
  #select 200 observations
  horizontal_selected[[i]] <-
    horizontal_lists[[i]][seq(
      from = 2,
      to = length(horizontal_lists[[i]]),
      length.out = 200
    )]
  horizontal_selected[[i]] <-
    cbind(horizontal_selected[[i]],
          names(horizontal_wall_thickness)[i],
          horizontal_wall_thickness[1, i])
  
}


horizontal_df <- as.data.frame(do.call(rbind, horizontal_selected))

horizontal_df[,1] <- as.numeric(horizontal_df[,1]) 

names(horizontal_df) <- c("Width_mm","Vessel","Ware")

horizontal_df$Ware <- as.factor(horizontal_df$Ware)
horizontal_df$Vessel <- as.factor(horizontal_df$Vessel)

horizontal_df <- horizontal_df[order(horizontal_df$Ware),]

# calculating average wall thickness

averages <- horizontal_df %>%
  group_by(Vessel) %>%
  summarise(mean_value = mean(Width_mm, na.rm = TRUE))

averages <- as.data.frame(averages)


# Get unique Vessel values
unique_vessels <- unique(horizontal_df$Vessel)

# Loop through unique vessels
for (i in seq_along(unique_vessels)) {
  vessel_name <- unique_vessels[i]  # Extract vessel name
  
  # Find the corresponding average value
  avg_value <- averages[averages$Vessel == vessel_name, 2]  # Assuming column 2 has the relevant average
  
  # Update the Width_mm_div column
  horizontal_df[horizontal_df$Vessel == vessel_name, "Width_mm_div"] <- 
    horizontal_df[horizontal_df$Vessel == vessel_name, "Width_mm"] / avg_value
}


# plotting wall thickness

boxplot_horizontal <- ggplot(horizontal_df, aes(x = reorder(Vessel,order(Ware)), y = Width_mm, fill = Ware)) +
  geom_boxplot(
    width = .6, 
    outlier.shape = NA
  ) +
  xlab("Vessel") +
  ylab("Width (mm) of horizontal section") +
  coord_cartesian(xlim = c(0, 30), clip = "off") + 
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



png(filename = "./Figures/Figure8a_1.png", width = 2400, height = 1600, res=300)
plot(boxplot_horizontal)
dev.off()




boxplot_horizontal_30 <- ggplot(horizontal_df, aes(x = reorder(Vessel,order(Ware)), y = Width_mm_div, fill = Ware)) +
  geom_boxplot(
    width = .6, 
    outlier.shape = NA
  ) +
  xlab("Vessel") +
  ylab("Horizontal width/Vessel's average horizontal width") +
  coord_cartesian(xlim = c(0, 30), clip = "off") + 
  coord_cartesian(ylim = c(0.5, 1.5)) +
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


png(filename = "./Figures/Figure8a_3.png", width = 2400, height = 1600, res=300)
plot(boxplot_horizontal_30)
dev.off()




#vertical wall thickness boxplot

vertical_wall_thickness <- read.csv("./Scripts/vertical_wall_thickness.csv", na.strings="",header = TRUE)


vertical_lists <- list()
vertical_selected <- list()

for (i in 1:ncol(vertical_wall_thickness)) {
  
  vertical_lists[[i]] <- vertical_wall_thickness[, i]
  
  #exclude NA observations
  vertical_lists[[i]] <-
    vertical_lists[[i]][!is.na(vertical_lists[[i]])]
  
  #select 200 observations
  vertical_selected[[i]] <-
    vertical_lists[[i]][seq(
      from = 2,
      to = length(vertical_lists[[i]]),
      length.out = 200
    )]
  vertical_selected[[i]] <-
    cbind(vertical_selected[[i]],
          names(vertical_wall_thickness)[i],
          vertical_wall_thickness[1, i])
  
}


vertical_df <- as.data.frame(do.call(rbind, vertical_selected))

vertical_df[,1] <- as.numeric(vertical_df[,1]) 

names(vertical_df) <- c("Width_mm","Vessel","Ware")

vertical_df$Ware <- as.factor(vertical_df$Ware)
vertical_df$Vessel <- as.factor(vertical_df$Vessel)

vertical_df <- vertical_df[order(vertical_df$Ware),]

# calculating average wall thickness

averages <- vertical_df %>%
  group_by(Vessel) %>%
  summarise(mean_value = mean(Width_mm, na.rm = TRUE))

averages <- as.data.frame(averages)


# Get unique Vessel values
unique_vessels <- unique(vertical_df$Vessel)

# Loop through unique vessels
for (i in seq_along(unique_vessels)) {
  vessel_name <- unique_vessels[i]  # Extract vessel name
  
  # Find the corresponding average value
  avg_value <- averages[averages$Vessel == vessel_name, 2]  # Assuming column 2 has the relevant average
  
  # Update the Width_mm_div column
  vertical_df[vertical_df$Vessel == vessel_name, "Width_mm_div"] <- 
    vertical_df[vertical_df$Vessel == vessel_name, "Width_mm"] / avg_value
}


#combine right and left 'arms' of the bowl together

vertical_df$Vessel <- as.character(vertical_df$Vessel)

vertical_df$Vessel_together <- substr(vertical_df$Vessel, 1, nchar(vertical_df$Vessel) - 2)

vertical_df$Vessel <- as.factor(vertical_df$Vessel)
vertical_df$Vessel_together <- as.factor(vertical_df$Vessel_together)

boxplot_vertical_comb <- ggplot(vertical_df, aes(x = reorder(Vessel_together,
                                                             order(Ware)), y = Width_mm, fill = Ware)) +
  geom_boxplot(
    width = .6, 
    outlier.shape = NA
  ) +
  xlab("Vessel") +
  ylab("Width (mm) of vertical section") +
  coord_cartesian(xlim = c(0, 30), clip = "off") + 
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


png(filename = "./Figures/Figure8a_2.png", width = 2400, height = 1600, res=300)
plot(boxplot_vertical_comb)
dev.off()



boxplot_vertical_30_comb <- ggplot(vertical_df, aes(x = reorder(Vessel_together,
                                                                order(Ware)), y = Width_mm_div, fill = Ware)) +
  geom_boxplot(
    width = .6, 
    outlier.shape = NA
  ) +
  xlab("Vessel") +
  ylab("Vertical width/Vessel's average vertical width") +
  coord_cartesian(xlim = c(0, 30), clip = "off") +
  coord_cartesian(ylim = c(0.5, 1.5)) +
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


png(filename = "./Figures/Figure8a_4.png", width = 2400, height = 1600, res=300)
plot(boxplot_vertical_30_comb)
dev.off()



------------------# Assessing roundness and circularity-------------------------

wall_thickness <- import('./Scripts/wall_thickness.csv')

wall_thickness <-  wall_thickness %>%
  mutate(ware = case_when(
    ware == "Piartal - 1" ~ "Piartal",
    ware == "Piartal - 2" ~ "Piartal",
    ware == "Tuza - Red slip" ~ "Tuza",
    TRUE ~ as.character(ware)  # Keep all other values as they are
  ))

wall_thickness$ware <- as.factor(wall_thickness$ware)
wall_thickness$subware <- as.factor(wall_thickness$subware)

wall_thickness$CV_horizontal_wall_thickness <- as.numeric(wall_thickness$CV_horizontal_wall_thickness)
wall_thickness$CV_vertical_wall_thickness_mm <- as.numeric(wall_thickness$CV_vertical_wall_thickness_mm)

manual_colors <- c("#1B9E77", "#7570B3", "#D95F02")

Roundness <- ggplot(wall_thickness, aes(x = ware, y = Round, 
                                        fill = ware, color = ware)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  labs(y = "Roundness") +
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") +
  theme_minimal()
plot(Roundness)

png(filename = "./Figures/Figure8b_1.png", width = 1350, height = 1000, res=300)
plot(Roundness)
dev.off()

Circ <- ggplot(wall_thickness, aes(x = ware, y = circ, 
                                   fill = ware, color = ware)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  labs(y = "Circularity") +
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") + 
  theme_minimal()

plot(Circ)

png(filename = "./Figures/Figure8b_2.png", width = 1350, height = 1000, res=300)
plot(Circ)
dev.off()






# -------------------Assessing wall thickness and circ parameters-------------------

anova_distance <- aov(Group_distance_to_centroid ~ ware, data = wall_thickness)
anova_height <- aov(Vessel_height_mm ~ ware, data = wall_thickness)
anova_circ <- aov(circ ~ ware, data = wall_thickness)
anova_round <- aov(Round ~ ware, data = wall_thickness)
anova_meanvwall <- aov(mean_vertical_wall_thickness_mm ~ ware, data = wall_thickness)
anova_meanhwall <- aov(mean_horizontal_wall_thickness ~ ware, data = wall_thickness)
anova_CVvwall <- aov(CV_vertical_wall_thickness_mm ~ ware, data = wall_thickness)
anova_CVhwall <- aov(CV_horizontal_wall_thickness ~ ware, data = wall_thickness)

summary(anova_distance)
summary(anova_height)
summary(anova_circ)
summary(anova_round)
summary(anova_meanvwall)
summary(anova_meanhwall)
summary(anova_CVvwall)
summary(anova_CVhwall)

leveneTest_distance <- leveneTest(Group_distance_to_centroid ~ ware, data = wall_thickness, center = mean)
leveneTest_height <- leveneTest(Vessel_height_mm ~ ware, data = wall_thickness, center = mean)
leveneTest_circ <- leveneTest(circ ~ ware, data = wall_thickness, center = mean)
leveneTest_round <- leveneTest(Round ~ ware, data = wall_thickness, center = mean)
leveneTest_meanvwall <- leveneTest(mean_vertical_wall_thickness_mm ~ ware, data = wall_thickness, center = mean)
leveneTest_meanhwall <- leveneTest(mean_horizontal_wall_thickness ~ ware, data = wall_thickness, center = mean)
leveneTest_CVvwall <- leveneTest(CV_vertical_wall_thickness_mm ~ ware, data = wall_thickness, center = mean)
leveneTest_CVhwall <- leveneTest(CV_horizontal_wall_thickness ~ ware, data = wall_thickness, center = mean)

print(leveneTest_distance)
print(leveneTest_height)
print(leveneTest_circ)
print(leveneTest_round)
print(leveneTest_meanvwall)
print(leveneTest_meanhwall)
print(leveneTest_CVvwall)
print(leveneTest_CVhwall)


# Pairwise t-tests for each dependent variable against ware
pairwise_distance <- pairwise.t.test(wall_thickness$Group_distance_to_centroid, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_height <- pairwise.t.test(wall_thickness$Vessel_height_mm, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_circ <- pairwise.t.test(wall_thickness$circ, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_round <- pairwise.t.test(wall_thickness$Round, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_meanvwall <- pairwise.t.test(wall_thickness$mean_vertical_wall_thickness_mm, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_meanhwall <- pairwise.t.test(wall_thickness$mean_horizontal_wall_thickness, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_CVvwall <- pairwise.t.test(wall_thickness$CV_vertical_wall_thickness_mm, wall_thickness$ware, p.adjust.method = "bonferroni")
pairwise_CVhwall <- pairwise.t.test(wall_thickness$CV_horizontal_wall_thickness, wall_thickness$ware, p.adjust.method = "bonferroni")

pairwise_distance
pairwise_height
pairwise_circ
pairwise_round
pairwise_meanvwall
pairwise_meanhwall
pairwise_CVvwall
pairwise_CVhwall

# Function to perform permutation test for variance
pairwise_permutation_variance <- function(data, response, group, num_permutations = 1000) {
  # Create all pairwise combinations of group levels
  combinations <- combn(levels(data[[group]]), 2, simplify = FALSE)
  
  # Initialize results dataframe
  results <- data.frame(Group1 = character(),
                        Group2 = character(),
                        Observed_F = numeric(),
                        Perm_p_value = numeric(),
                        stringsAsFactors = FALSE)
  
  # Perform pairwise permutation test for variance
  for (comb in combinations) {
    subset_data <- data %>% filter(.data[[group]] %in% comb)
    group1 <- subset_data %>% filter(.data[[group]] == comb[1]) %>% pull(.data[[response]])
    group2 <- subset_data %>% filter(.data[[group]] == comb[2]) %>% pull(.data[[response]])
    
    # Observed F-statistic for variance ratio
    observed_F <- var(group1) / var(group2)
    
    # Permutation test
    combined <- c(group1, group2)
    group_labels <- c(rep(comb[1], length(group1)), rep(comb[2], length(group2)))
    
    perm_F <- replicate(num_permutations, {
      permuted_labels <- sample(group_labels)
      perm_group1 <- combined[permuted_labels == comb[1]]
      perm_group2 <- combined[permuted_labels == comb[2]]
      var(perm_group1) / var(perm_group2)
    })
    
    # Calculate permutation p-value
    perm_p_value <- mean(abs(perm_F) >= abs(observed_F))
    
    # Save results
    results <- rbind(results, data.frame(Group1 = comb[1],
                                         Group2 = comb[2],
                                         Observed_F = observed_F,
                                         Perm_p_value = perm_p_value))
  }
  
  # Adjust p-values for multiple comparisons
  results$Adjusted_p_value <- p.adjust(results$Perm_p_value, method = "bonferroni")
  
  return(results)
}

# Variables to test from wall_thickness dataframe
wall_thickness_vars <- c("Group_distance_to_centroid", 
                         "circ", 
                         "Round", 
                         "mean_vertical_wall_thickness_mm", 
                         "mean_horizontal_wall_thickness", 
                         "CV_vertical_wall_thickness_mm", 
                         "CV_horizontal_wall_thickness")

# Perform tests for wall_thickness
wall_thickness_results <- lapply(wall_thickness_vars, function(var) {
  cat("\nTesting variance for:", var, "\n")
  pairwise_permutation_variance(wall_thickness, var, "ware")
})

names(wall_thickness_results) <- wall_thickness_vars
wall_thickness_results

#-------------------------------Density plots-----------------------------------

## All figures produced in this section are included in Table 4

# Load necessary library
library(ggplot2)
library(dplyr)

# Import DGroup values

DGroup <- import("EFA_output.csv")

head(DGroup)
DGroup$Vessel <- as.factor(DGroup$Vessel)
is.factor(DGroup$Vessel)
DGroup$ware <- as.factor(DGroup$ware)
is.factor(DGroup$ware)
DGroup <- as.data.frame(DGroup)
is.data.frame(DGroup)

DGroup_three <-  DGroup %>%
  mutate(ware = case_when(
    ware == "Piartal - 1" ~ "Piartal",
    ware == "Piartal - 2" ~ "Piartal",
    ware == "Tuza - Red slip" ~ "Tuza",
    TRUE ~ as.character(ware)  
  ))

##SHOWING VARIABILITY - Mean value plotted

group_distance_density_plot <- ggplot(DGroup_three, aes(x = Distance, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  geom_vline(data = wall_thickness %>% group_by(ware) %>% 
               summarise(mean_value = mean(Group_distance_to_centroid, na.rm = TRUE)),
             aes(xintercept = mean_value, color = ware), linetype = "dashed", size = 0.7) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Group Distance to Centroid",
    y = "Density"
  )

group_distance_density_plot

png(filename = "group_distance_density_plot.png", width = 1600, height = 1200, res=300)
plot(group_distance_density_plot)
dev.off()

Circ_density_plot <- ggplot(wall_thickness, aes(x = circ, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  geom_vline(data = wall_thickness %>% group_by(ware) %>% 
               summarise(mean_value = mean(circ, na.rm = TRUE)),
             aes(xintercept = mean_value, color = ware), linetype = "dashed", size = 0.7) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Circularity",
    y = "Density"
  )

png(filename = "Circ_density_plot.png", width = 1600, height = 1200, res=300)
plot(Circ_density_plot)
dev.off()

Round_density_plot <- ggplot(wall_thickness, aes(x = Round, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  geom_vline(data = wall_thickness %>% group_by(ware) %>% 
               summarise(mean_value = mean(Round, na.rm = TRUE)),
             aes(xintercept = mean_value, color = ware), linetype = "dashed", size = 0.7) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Roundness",
    y = "Density"
  )

png(filename = "Round_density_plot.png", width = 1600, height = 1200, res=300)
plot(Round_density_plot)
dev.off()


CVvthickness_density_plot <- ggplot(wall_thickness, aes(x = CV_vertical_wall_thickness_mm, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  geom_vline(data = wall_thickness %>% group_by(ware) %>% 
               summarise(mean_value = mean(CV_vertical_wall_thickness_mm, na.rm = TRUE)),
             aes(xintercept = mean_value, color = ware), linetype = "dashed", size = 0.7) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "CV Vertical wall thickness",
    y = "Density"
  )

png(filename = "CVvthickness_density_plot.png", width = 1600, height = 1200, res=300)
plot(CVvthickness_density_plot)
dev.off()

CVhthickness_density_plot <- ggplot(wall_thickness, aes(x = CV_horizontal_wall_thickness, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  geom_vline(data = wall_thickness %>% group_by(ware) %>% 
               summarise(mean_value = mean(CV_horizontal_wall_thickness, na.rm = TRUE)),
             aes(xintercept = mean_value, color = ware), linetype = "dashed", size = 0.7) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "CV Horizontal wall thickness",
    y = "Density"
  )

png(filename = "CVhthickness_density_plot.png", width = 1600, height = 1200, res=300)
plot(CVhthickness_density_plot)
dev.off()


CVhthickness_density_plot <- ggplot(wall_thickness, aes(x = CV_horizontal_wall_thickness, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  geom_vline(data = wall_thickness %>% group_by(ware) %>% 
               summarise(mean_value = mean(CV_horizontal_wall_thickness, na.rm = TRUE)),
             aes(xintercept = mean_value, color = ware), linetype = "dashed", size = 0.7) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "CV Horizontal wall thickness",
    y = "Density"
  )

png(filename = "CVhthickness_density_plot.png", width = 1600, height = 1200, res=300)
plot(CVhthickness_density_plot)
dev.off()


##SHOWING VARIABILITY - Range (90%) values plotted


# Group_distance
Group_distance_range_data <- wall_thickness %>%
  group_by(ware) %>%
  summarise(
    P5 = quantile(Group_distance_to_centroid, 0.05, na.rm = TRUE),
    P95 = quantile(Group_distance_to_centroid, 0.95, na.rm = TRUE)
  )

horizontal_range_data <- horizontal_df %>%
  group_by(Ware) %>%
  summarise(
    P5 = quantile(Width_mm_div, 0.05, na.rm = TRUE),
    P95 = quantile(Width_mm_div, 0.95, na.rm = TRUE)
  )



horizontal_df <-  horizontal_df %>%
  mutate(Ware = case_when(
    Ware == "Piartal - 1" ~ "Piartal",
    Ware == "Piartal - 2" ~ "Piartal",
    Ware == "Tuza - Red slip" ~ "Tuza",
    TRUE ~ as.character(Ware)  
  ))

hthickenss_density_plot <- ggplot(horizontal_df, aes(x = Width_mm_div, fill = Ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors)+
  theme_minimal() +
  # Add vertical lines for 5th and 95th percentiles
  geom_vline(data = horizontal_range_data, aes(xintercept = P5, color = Ware), linetype = "dotted", size = 1) +
  geom_vline(data = horizontal_range_data, aes(xintercept = P95, color = Ware), linetype = "dotted", size = 1) +
  # Add staggered tie lines between P5 and P95
  geom_segment(data = horizontal_range_data, 
               aes(x = P5, xend = P95, y = as.numeric(factor(Ware)) * 2, yend = as.numeric(factor(Ware)) * 2, color = Ware), 
               size = 1) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "relative wall thickness",
    y = "Density"
  )

hthickenss_density_plot 


distance_range_data <- DGroup_three %>%
  group_by(ware) %>%
  summarise(
    P5 = quantile(Distance, 0.023, na.rm = TRUE),
    P95 = quantile(Distance, 0.977, na.rm = TRUE)
  )


group_distance_range_density_plot <- ggplot(DGroup_three, aes(x = Distance, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  # Add vertical lines for 5th and 95th percentiles
  geom_vline(data = distance_range_data, aes(xintercept = P5, color = ware), linetype = "dotted", size = 1) +
  geom_vline(data = distance_range_data, aes(xintercept = P95, color = ware), linetype = "dotted", size = 1) +
  # Add staggered tie lines between P5 and P95
  geom_segment(data = distance_range_data, 
               aes(x = P5, xend = P95, y = as.numeric(factor(ware)) * 15, yend = as.numeric(factor(ware)) * 15, color = ware), 
               size = 1) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Group distance to centroid",
    y = "Density"
  )

group_distance_range_density_plot

png(filename = "distance_range_density_plot.png", width = 1600, height = 1200, res=300)
plot(group_distance_range_density_plot)
dev.off()

mvthickness_range_data <- wall_thickness %>%
  group_by(ware) %>%
  summarise(
    P5 = quantile(mean_vertical_wall_thickness_mm, 0.05, na.rm = TRUE),
    P95 = quantile(mean_vertical_wall_thickness_mm, 0.95, na.rm = TRUE)
  )

mvthickness_density_plot <- ggplot(wall_thickness, aes(x = mean_vertical_wall_thickness_mm, fill = ware)) +
  geom_density(alpha = 0.3) +
  # Add vertical lines for 5th and 95th percentiles
  geom_vline(data = mvthickness_range_data, aes(xintercept = P5, color = ware), linetype = "dotted", size = 1) +
  geom_vline(data = mvthickness_range_data, aes(xintercept = P95, color = ware), linetype = "dotted", size = 1) +
  # Add staggered tie lines between P5 and P95
  geom_segment(data = mvthickness_range_data, 
               aes(x = P5, xend = P95, y = as.numeric(factor(ware)) * 0.6, yend = as.numeric(factor(ware)) * 0.6, color = ware), 
               size = 1) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Vertical wall thickness (mm)",
    y = "Density"
  )

png(filename = "mvthickness_density_plot.png", width = 1600, height = 1200, res=300)
plot(mvthickness_density_plot)
dev.off()

mhthickness_range_data <- wall_thickness %>%
  group_by(ware) %>%
  summarise(
    P5 = quantile(mean_horizontal_wall_thickness, 0.05, na.rm = TRUE),
    P95 = quantile(mean_horizontal_wall_thickness, 0.95, na.rm = TRUE)
  )

mhthickness_density_plot <- ggplot(wall_thickness, aes(x = mean_horizontal_wall_thickness, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  # Add vertical lines for 5th and 95th percentiles
  geom_vline(data = mhthickness_range_data, aes(xintercept = P5, color = ware), linetype = "dotted", size = 1) +
  geom_vline(data = mhthickness_range_data, aes(xintercept = P95, color = ware), linetype = "dotted", size = 1) +
  # Add staggered tie lines between P5 and P95
  geom_segment(data = mhthickness_range_data, 
               aes(x = P5, xend = P95, y = as.numeric(factor(ware)) * 0.3, yend = as.numeric(factor(ware)) * 0.3, color = ware), 
               size = 1) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Horizontal wall thickness (mm)",
    y = "Density"
  )

png(filename = "mhthickness_density_plot.png", width = 1600, height = 1200, res=300)
plot(mhthickness_density_plot)
dev.off()

Circ_range_data <- wall_thickness %>%
  group_by(ware) %>%
  summarise(
    P5 = quantile(circ, 0.05, na.rm = TRUE),
    P95 = quantile(circ, 0.95, na.rm = TRUE)
  )

Circ_range_density_plot <- ggplot(wall_thickness, aes(x = circ, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  # Add vertical lines for 5th and 95th percentiles
  geom_vline(data = Circ_range_data, aes(xintercept = P5, color = ware), linetype = "dotted", size = 1) +
  geom_vline(data = Circ_range_data, aes(xintercept = P95, color = ware), linetype = "dotted", size = 1) +
  # Add staggered tie lines between P5 and P95
  geom_segment(data = Circ_range_data, 
               aes(x = P5, xend = P95, y = as.numeric(factor(ware)) * 80, yend = as.numeric(factor(ware)) * 80, color = ware), 
               size = 1) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Circularity",
    y = "Density"
  )

png(filename = "Circ_range_density_plot.png", width = 1600, height = 1200, res=300)
plot(Circ_range_density_plot)
dev.off()


Round_range_data <- wall_thickness %>%
  group_by(ware) %>%
  summarise(
    P5 = quantile(Round, 0.05, na.rm = TRUE),
    P95 = quantile(Round, 0.95, na.rm = TRUE)
  )

Round_range_density_plot <- ggplot(wall_thickness, aes(x = Round, fill = ware)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = manual_colors) + 
  theme_minimal() +
  # Add vertical lines for 5th and 95th percentiles
  geom_vline(data = Round_range_data, aes(xintercept = P5, color = ware), linetype = "dotted", size = 1) +
  geom_vline(data = Round_range_data, aes(xintercept = P95, color = ware), linetype = "dotted", size = 1) +
  # Add staggered tie lines between P5 and P95
  geom_segment(data = Round_range_data, 
               aes(x = P5, xend = P95, y = as.numeric(factor(ware)) * 15, yend = as.numeric(factor(ware)) * 15, color = ware), 
               size = 1) +
  scale_color_manual(values = manual_colors) + 
  labs(
    x = "Roundness",
    y = "Density"
  )

png(filename = "Round_range_density_plot.png", width = 1600, height = 1200, res=300)
plot(Round_range_density_plot)
dev.off()

# ------------------------------Acknowledgement---------------------------------

citation()
citation("car")
citation("ggplot2")
citation("ggbeeswarm")
citation("ggdist")
citation("gghalves")
citation("dplyr")
citation("tidyr")
