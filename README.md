# GMM_intravessel

This repository contains the data and scripts used in the article:

Catherine Klesner, Rosie Crawford, Jasmine Vieri, & Marcos Martinón-Torres (2026) **A method to assess inter- and intra-vessel shape variation in pottery using outline-based geometric morphometrics (GMM)** *Antiquity*

---

## Repository structure

The repository is organised into four main directories:

- **Scripts/** – R scripts used to process the data, run analyses, and generate outputs (including non-figure outputs).
- **Data/** – Original input datasets used in the study (e.g., `.TPS` and `.csv` files), plus any compiled/cleaned datasets required by the scripts.
- **vessel_cross_sections/** - contains all the cross-sections
- **Figures/** – All figures included in the manuscript (both those generated via the R scripts and those prepared by other methods).
- **Supplementary_Materials/** – Supplementary documents and figures referenced in the manuscript, as well as `3D_model_links.csv` which contains links to the 3D models analysed in this study (hosted on Sketchfab).

---

## Overview of analyses and scripts

This repository includes three main workflows:

1. **Calculation of completeness of individual cross-sections**
2. **GMM (outline-based) intravessel variation**
3. **Metric intravessel variation**


## Overview

This repository supports three linked workflows that move from (1) evaluating the completeness of digitised cross-sections, to (2) assessing metric intravessel variability, and (3) analysing outline shape variability using geometric morphometrics (Elliptic Fourier Analysis and PCA). The scripts are designed so that intermediate outputs (e.g., completeness estimates and EFA coefficient tables) can be reused across workflows, and the key figures in the manuscript are reproducible directly from the R code.

---

## 1: Calculation of completeness of individual cross-sections

The script **`Scripts/cross_section_completeness.R`** calculates how complete each digitised vessel cross-section is, based on illustration conventions for reconstructing missing portions of archaeological profiles. Prior to running the script, cross-section images must be processed in an image-editing program (we used Inkscape). Missing portions of each cross-section are infilled in **RED**, following standard archaeological illustration practice. The script then reads the processed image set from `./vessel_cross_sections/`.

Each image is converted to a raster and the script counts pixels classified as **red**, **black**, and **white**. Completeness is calculated, and the completeness table is saved as  `Scripts/color_counts_combined.csv`.

This output is later used in the GMM workflow to filter cross-sections by preservation quality (e.g., retaining only outlines >99% complete).

---

## 2: GMM intravessel variation (outline-based shape analysis)

The script **`Scripts/GMM_intravessel_variation.R`** performs the outline-based geometric morphometric analysis of vessel cross-sections, using the `Momocs` framework for Elliptic Fourier Analysis (EFA) and PCA.

All outlines were digitised using `tpsUtil` and `tpsDig2`, and imported into R as TPS coordinate files with associated CSV metadata. The workflow uses two outline datasets:

1. A single-cross-section dataset used for Figure 1a:
   - `Data/single_cross_section.TPS`
   - `Data/single_cross_section.csv`

2. The main intravessel dataset:
   - `Data/intravessel30.TPS`
   - `Data/intravessel30_db.csv`

Before analysis, the script checks that TPS coordinate names match the metadata order (critical for correct grouping), then constructs `Momocs::Out()` objects.

Completeness values (`percent_complete`) generated in `Scripts/cross_section_completeness.R` are used to subset the intravessel dataset. Cross-sections with `percent_complete > 0.99` are retained for the core analysis, producing an “over 99% complete” outline dataset used throughout downstream EFA and PCA.

Outlines are first normalised to produce a consistent outline set for EFA. The script calibrates harmonic power and reconstructions, then computes EFA coefficients. While harmonic calibration is explored (`Figures/Figure3a.png` and `Figures/Figure3b.png`), the main analysis is run using **30 harmonics** (`Figures/Figure4a.png` and `Figures/Figure4b.png`) 

PCA is applied to the EFA coefficient set, and variation is visualised in PCA morphospace (PC1 vs PC2) by ware (`Figures/Figure5_1.pdf`) and by vessel (`Figures/Figure5_2.pdf`). To quantify how much cross-sections vary within vessels (and compare dispersion among groups), the script computes Euclidean distances among EFA coefficients and uses `vegan::betadisper` to estimate each outline’s distance to its group centroid (**DGroup**). The EFA coefficients are exported for transparency and reuse (`Scripts/over_99_EFA_30.txt`), and the DGroup values are also saved (`Scripts/EFA_output.csv`) and prepared as a summary figure (`Figures/Figure6.png`). This distance output is reused in the metric variability workflow.

Finally, the script explores intravessel shape variation by selecting a set of representative vessels (chosen to reflect low, high, and median dispersion values across wares). PCA morphospace plots for these individual vessels are saved as:

- `Figures/Figure7a.png` – `CA230365`
- `Figures/Figure7b.png` – `CA230377`
- `Figures/Figure7c.png` – `CA230417`
- `Figures/Figure7d.png` – `CA230347`
- `Figures/Figure7e.png` – `CA230715`
- `Figures/Figure7f.png` – `CA230335`

---

## 3: Metric intravessel variation

The script **`Scripts/metric_intravessel_variation.R`** assesses intravessel variability using metric measurements and circularity/roundness descriptors derived from cross-sections. The script imports horizontal wall thickness measurements from `Data/horizontal_wall_thickness.csv`. For each vessel a subsample of **200** measurements is selected to standardise comparisons across vessels. It then calculates vessel-specific mean thickness values and creates a normalised thickness metric (`Width_mm_div`), defined as each measurement divided by the vessel mean. This produces two complementary views of thickness variation: absolute values (mm) and relative variation within each vessel. Figures are saved to `Figures/Figure8b_1.png` (boxplot of absolute horizontal wall thickness by vessel and ware) and `Figures/Figure8b_2.png` (boxplot of relative thickness (measurement / vessel mean)). 

Roundness and circularity metrics calculated using ImageJ are imported from `Data/wall_thickness.csv`. The script then visualises the distributions of roundness and circularity by ware, and figures are saved to `Figures/Figure8a_1.png` (roundness distributions by ware) and `Figures/Figure8a_2.png` (circularity distributions by ware). 

### Statistical testing and variability plots

The script evaluates whether variability differs among wares using a combination of Levene’s tests (variance homogeneity), ANOVA and Welch one-way tests, pairwise t-tests (Bonferroni-adjusted), and pairwise permutation tests for variance ratios. It also generates density plots to visualise variability across groups, which ara saved to `Figures/Table_4/` (multiple `.png` files). 

Some of these density plots rely on distance-to-centroid outputs produced by the GMM workflow, imported from `Scripts/EFA_output.csv`

---

## Data

The `Data/` directory contains the original data files used in the analyses, including:
- **TPS files** containing digitised cross-section outlines (generated using `tpsUtil` and `tpsDig2`)
- **CSV files** containing metadata and measurements (e.g., ware assignments, vessel IDs, wall thickness measurements, and completeness values)

---

## Figures

The `Figures/` directory contains all figures included in the manuscript:
- Figures generated directly by the R scripts (saved as `.png` and `.pdf`)
- Figures prepared by other methods but included in the manuscript

`Figures/Table_4/` contains the supplementary variability plots generated in R (e.g., density plots referenced in Table 4).

---

## Supplementary materials

The `Supplementary_Materials/` directory contains all supplementary documents and figures referenced in the article.

---

## R session info

The code was run using R version 4.3.1 

The following R packages are required to execute the code:

- magick_2.9.0 (Ooms, 2025)  
- raster_3.6-32 (Hijmans, 2025)  
- rio_1.2.4 (Chan et al., 2023)  
- ggplot2_4.0.2 (Wickham, 2016)  
- dplyr_1.1.4 (Wickham et al., 2023)  
- tidyr_1.3.1 (Wickham et al., 2024)  
- Momocs_1.5.0 (Bonhomme et al., 2014)  
- vegan_2.7-2 (Oksanen et al., 2025)  
- car_3.1-5 (Fox and Weisberg, 2019)  
- ggbeeswarm_0.7.2 (Clarke et al., 2025)  
- ggdist_3.3.3 (Kay, 2025)  
- gghalves_0.1.4 (Tiedemann, 2025)

Citations:

Bonhomme V, Picq S, Gaucherel C, Claude J (2014). _Momocs: Outline Analysis Using R_, volume 56 number 13. <https://www.jstatsoft.org/v56/i13/>.

Chan C, Leeper T, Becker J, Schoch D (2023). _rio: A Swiss-army knife for data file I/O_. <https://cran.r-project.org/package=rio>.

Clarke E, Sherrill-Mix S, Dawson C (2023). _ggbeeswarm: Categorical Scatter (Violin Point) Plots_. R package version 0.7.2, <https://CRAN.R-project.org/package=ggbeeswarm>.

Fox J, Weisberg S (2019). _An R Companion to Applied Regression_, Third edition. Sage, Thousand Oaks CA. <https://www.john-fox.ca/Companion/>.

Hijmans R (2025). _raster: Geographic Data Analysis and Modeling_. R package version 3.6-32, <https://CRAN.R-project.org/package=raster>.

Kay M (2025). _ggdist: Visualizations of Distributions and Uncertainty_. <https://doi.org/10.5281/zenodo.3879620>, R package version 3.3.3,
<https://mjskay.github.io/ggdist/>.

Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Borman T, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, Martino C, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2025). _vegan: Community Ecology Package_. R package version 2.7-2, <https://CRAN.R-project.org/package=vegan>.

Ooms J (2025). _magick: Advanced Graphics and Image-Processing in R_. R package version 2.9.0, <https://CRAN.R-project.org/package=magick>.

R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
  
Tiedemann F (2022). _gghalves: Compose Half-Half Plots Using Your Favourite Geoms_. R package version 0.1.4, <https://CRAN.R-project.org/package=gghalves>.

Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis.* Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.4, <https://CRAN.R-project.org/package=dplyr>.

Wickham H, Vaughan D, Girlich M (2024). _tidyr: Tidy Messy Data_. R package version 1.3.1, <https://CRAN.R-project.org/package=tidyr>.

---

## Funding

This project has received funding from the **European Research Council (ERC)** under the European Union’s **Horizon 2020** research and innovation programme (**Grant agreement No. 101021480**, granted to **Marcos Martinón-Torres**).

