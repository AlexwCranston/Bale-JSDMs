## Understanding Impacts from Domestic Livestock on Ungulate Communities in Bale Mountains National Park

### Overview

Authors: TBC

This repository contains all R scripts for a paper (in prep) using joint species distribution modelling using Hierarchical Modelling of Species Communities (HMSC) using transect data collected from the Bale Mountains National Park. The project aims to assess the impacts of livestock grazing on the ungulate community of Bale Mountains National Park,

Data associated with this paper to be used in conjunction with this code is available at: TBC

To cite this paper: TBC

To cite this repo: TBC

### Data Processing

The repository uses survey data on wild and domestic mammals collected in Bale Mountains National Park as its input. Data processing of the raw data is handled "Preliminary Data Processing.R" script. This script cleans the data and then bins all observations into a regular grid covering the survey area. It also extracts environmental variables (habitat, elevation, minimum distance from road) for this grid for use in later modelling.

### Spatial Analyses

To visualise the survey data, use "IDW Interpolation.R" to produce basic maps showing heatmaps of raw count data for all species in the survey area.

Co-occurrence analysis of the raw data  can be conducted in "EcoSimR.R". This script uses the EcoSim package to simulate a null distribution of checkerboard scores (measuring extent of species co-occurrences) to which the true data is then compared.

The main spatial analyses are handled by the following scripts. Run these in the following order:

1. "Define HMSC Models_Wet SeasonONLY.R" - This script defines the joint species distribution models using an HMSC framework and exports them as a list of unfitted models.
2. "Fit Models_Wet Season Only.R" - This script takes the list of unfitted models and fits them using MCMC sampling to estimate all parameters. It exports lists of fitted models.
3. "Evaluate Convergence HMSC_Wet SeasonOnly.R" - This script takes fitted models and evaluates convergence of the MCMC sampling.
4. "Evaluate Parameters HMSC_Wet SeasonOnly.R" - This script takes fitted models and evaluates the estimated parameters and produces some basic plots showing the model outputs.

### Population Estimates and Temporal Analyses

Population estimates for all species using survey data are conducted using "DistanceSampling_PopulationEstimates.R". This script outputs a .csv file of abundance estimates for all species over all years. "Time Series.R" uses this csv to run an additional HMSC that is temporally explicit (i.e. asking if years with higher abundances of one species have higher or lower abundances of another species?).

### Session Info

Please see below for the output of devtools::session_info() used to perform the analyses for this paper.

```
- Session info --------------------------------------------------------------------------------------------------------
 setting  value
 version  R version 4.4.2 (2024-10-31 ucrt)
 os       Windows 11 x64 (build 26100)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United Kingdom.1252
 ctype    English_United Kingdom.1252
 tz       Europe/London
 date     2025-04-22
 rstudio  2021.09.0+351 Ghost Orchid (desktop)
 pandoc   NA
 quarto   NA

- Packages ------------------------------------------------------------------------------------------------------------
 package       * version    date (UTC) lib source
 abind         * 1.4-8      2024-09-12 [1] CRAN (R 4.4.1)
 ape             5.8-1      2024-12-16 [1] CRAN (R 4.4.2)
 automap       * 1.1-12     2024-09-03 [1] CRAN (R 4.4.2)
 backports       1.5.0      2024-05-23 [1] CRAN (R 4.4.0)
 base64enc       0.1-3      2015-07-28 [1] CRAN (R 4.4.0)
 BayesLogit      2.1        2019-09-26 [1] CRAN (R 4.4.0)
 bitops          1.0-9      2024-10-03 [1] CRAN (R 4.4.1)
 broom         * 1.0.7      2024-09-26 [1] CRAN (R 4.4.1)
 cachem          1.1.0      2024-05-16 [1] CRAN (R 4.4.2)
 caTools         1.18.3     2024-09-04 [1] CRAN (R 4.4.1)
 class           7.3-22     2023-05-03 [2] CRAN (R 4.4.2)
 classInt        0.4-10     2023-09-05 [1] CRAN (R 4.4.2)
 cli             3.6.4      2025-02-13 [1] CRAN (R 4.4.3)
 coda          * 0.19-4.1   2024-01-31 [1] CRAN (R 4.4.2)
 codetools       0.2-20     2024-03-31 [2] CRAN (R 4.4.2)
 colorspace    * 2.1-1      2024-07-26 [1] CRAN (R 4.4.1)
 corrplot      * 0.95       2024-10-14 [1] CRAN (R 4.4.2)
 crosstalk       1.2.1      2023-11-23 [1] CRAN (R 4.4.2)
 DALEX           2.4.3      2023-01-15 [1] CRAN (R 4.4.3)
 data.table      1.16.4     2024-12-06 [1] CRAN (R 4.4.2)
 DBI             1.2.3      2024-06-02 [1] CRAN (R 4.4.2)
 devtools        2.4.5      2022-10-11 [1] CRAN (R 4.4.2)
 dials         * 1.4.0      2025-02-13 [1] CRAN (R 4.4.3)
 DiceDesign      1.10       2023-12-07 [1] CRAN (R 4.4.3)
 digest          0.6.37     2024-08-19 [1] CRAN (R 4.4.2)
 Distance      * 2.0.0      2024-10-24 [1] CRAN (R 4.4.3)
 doParallel    * 1.0.17     2022-02-07 [1] CRAN (R 4.4.2)
 dotCall64       1.2        2024-10-04 [1] CRAN (R 4.4.2)
 dplyr         * 1.1.4      2023-11-17 [1] CRAN (R 4.4.1)
 e1071           1.7-16     2024-09-16 [1] CRAN (R 4.4.2)
 EcoSimR       * 0.1.0      2015-04-03 [1] CRAN (R 4.4.2)
 ellipsis        0.3.2      2021-04-29 [1] CRAN (R 4.4.2)
 farver          2.1.2      2024-05-13 [1] CRAN (R 4.4.1)
 fastmap         1.2.0      2024-05-15 [1] CRAN (R 4.4.2)
 fields          16.3       2024-09-30 [1] CRAN (R 4.4.2)
 FNN             1.1.4.1    2024-09-22 [1] CRAN (R 4.4.2)
 forcats       * 1.0.0      2023-01-29 [1] CRAN (R 4.4.1)
 foreach       * 1.5.2      2022-02-02 [1] CRAN (R 4.4.2)
 fs              1.6.5      2024-10-30 [1] CRAN (R 4.4.2)
 furrr           0.3.1      2022-08-15 [1] CRAN (R 4.4.3)
 future          1.34.0     2024-07-29 [1] CRAN (R 4.4.2)
 future.apply    1.11.3     2024-10-27 [1] CRAN (R 4.4.2)
 generics        0.1.3      2022-07-05 [1] CRAN (R 4.4.1)
 ggplot2       * 3.5.1      2024-04-23 [1] CRAN (R 4.4.1)
 ggspatial     * 1.1.9      2023-08-17 [1] CRAN (R 4.4.3)
 globals         0.16.3     2024-03-08 [1] CRAN (R 4.4.0)
 glue            1.8.0      2024-09-30 [1] CRAN (R 4.4.1)
 gower           1.0.2      2024-12-17 [1] CRAN (R 4.4.2)
 GPfit           1.0-8      2019-02-08 [1] CRAN (R 4.4.3)
 gplots        * 3.2.0      2024-10-05 [1] CRAN (R 4.4.1)
 gstat         * 2.1-2      2024-09-05 [1] CRAN (R 4.4.2)
 gtable          0.3.6      2024-10-25 [1] CRAN (R 4.4.2)
 gtools          3.9.5      2023-11-20 [1] CRAN (R 4.4.1)
 hardhat         1.4.1      2025-01-31 [1] CRAN (R 4.4.3)
 hms             1.1.3      2023-03-21 [1] CRAN (R 4.4.2)
 Hmsc          * 3.0-13     2022-08-11 [1] CRAN (R 4.4.2)
 htmltools       0.5.8.1    2024-04-04 [1] CRAN (R 4.4.2)
 htmlwidgets     1.6.4      2023-12-06 [1] CRAN (R 4.4.2)
 httpuv          1.6.15     2024-03-26 [1] CRAN (R 4.4.2)
 ids           * 1.0.1      2017-05-31 [1] CRAN (R 4.4.2)
 igraph          2.1.4      2025-01-23 [1] CRAN (R 4.4.2)
 infer         * 1.0.7      2024-03-25 [1] CRAN (R 4.4.3)
 intervals       0.15.5     2024-08-23 [1] CRAN (R 4.4.1)
 ipred           0.9-15     2024-07-18 [1] CRAN (R 4.4.3)
 iterators     * 1.0.14     2022-02-05 [1] CRAN (R 4.4.2)
 janitor         2.2.1      2024-12-22 [1] CRAN (R 4.4.3)
 KernSmooth      2.23-24    2024-05-17 [2] CRAN (R 4.4.2)
 kknn          * 1.3.1      2016-03-26 [1] CRAN (R 4.4.3)
 later           1.4.1      2024-11-27 [1] CRAN (R 4.4.2)
 lattice         0.22-6     2024-03-20 [2] CRAN (R 4.4.2)
 lava            1.8.1      2025-01-12 [1] CRAN (R 4.4.3)
 leafem          0.2.3      2023-09-17 [1] CRAN (R 4.4.3)
 leaflet         2.2.2      2024-03-26 [1] CRAN (R 4.4.3)
 lhs             1.2.0      2024-06-30 [1] CRAN (R 4.4.3)
 lifecycle       1.0.4      2023-11-07 [1] CRAN (R 4.4.1)
 listenv         0.9.1      2024-01-29 [1] CRAN (R 4.4.2)
 lubridate     * 1.9.4      2024-12-08 [1] CRAN (R 4.4.2)
 magrittr        2.0.3      2022-03-30 [1] CRAN (R 4.4.1)
 maps            3.4.2.1    2024-11-10 [1] CRAN (R 4.4.2)
 mapview       * 2.11.2     2023-10-13 [1] CRAN (R 4.4.3)
 MASS          * 7.3-61     2024-06-13 [2] CRAN (R 4.4.2)
 Matrix          1.7-1      2024-10-18 [2] CRAN (R 4.4.2)
 MatrixModels    0.5-3      2023-11-06 [1] CRAN (R 4.4.2)
 matrixStats     1.5.0      2025-01-07 [1] CRAN (R 4.4.2)
 mcmc            0.9-8      2023-11-16 [1] CRAN (R 4.4.2)
 MCMCpack        1.7-1      2024-08-27 [1] CRAN (R 4.4.2)
 memoise         2.0.1      2021-11-26 [1] CRAN (R 4.4.2)
 mgcv            1.9-1      2023-12-21 [2] CRAN (R 4.4.2)
 mime            0.12       2021-09-28 [1] CRAN (R 4.4.0)
 miniUI          0.1.1.1    2018-05-18 [1] CRAN (R 4.4.2)
 modeldata     * 1.4.0      2024-06-19 [1] CRAN (R 4.4.3)
 mrds          * 3.0.0      2024-10-23 [1] CRAN (R 4.4.3)
 munsell         0.5.1      2024-04-01 [1] CRAN (R 4.4.1)
 nlme            3.1-166    2024-08-14 [2] CRAN (R 4.4.2)
 nloptr          2.1.1      2024-06-25 [1] CRAN (R 4.4.1)
 nnet            7.3-19     2023-05-03 [2] CRAN (R 4.4.2)
 numDeriv        2016.8-1.1 2019-06-06 [1] CRAN (R 4.4.0)
 optimx          2024-12.2  2024-12-10 [1] CRAN (R 4.4.3)
 parallelly      1.42.0     2025-01-30 [1] CRAN (R 4.4.1)
 parsnip       * 1.3.1      2025-03-12 [1] CRAN (R 4.4.3)
 patchwork     * 1.3.0      2024-09-16 [1] CRAN (R 4.4.2)
 pillar          1.10.1     2025-01-07 [1] CRAN (R 4.4.2)
 pkgbuild        1.4.6      2025-01-16 [1] CRAN (R 4.4.2)
 pkgconfig       2.0.3      2019-09-22 [1] CRAN (R 4.4.1)
 pkgload         1.4.0      2024-06-28 [1] CRAN (R 4.4.2)
 plyr            1.8.9      2023-10-02 [1] CRAN (R 4.4.2)
 png             0.1-8      2022-11-29 [1] CRAN (R 4.4.0)
 pracma          2.4.4      2023-11-10 [1] CRAN (R 4.4.2)
 pROC            1.18.5     2023-11-01 [1] CRAN (R 4.4.2)
 prodlim         2024.06.25 2024-06-24 [1] CRAN (R 4.4.3)
 profvis         0.4.0      2024-09-20 [1] CRAN (R 4.4.2)
 promises        1.3.2      2024-11-28 [1] CRAN (R 4.4.2)
 proxy           0.4-27     2022-06-09 [1] CRAN (R 4.4.2)
 purrr         * 1.0.4      2025-02-05 [1] CRAN (R 4.4.3)
 quantreg        6.00       2025-01-29 [1] CRAN (R 4.4.2)
 R6              2.6.1      2025-02-15 [1] CRAN (R 4.4.2)
 raster        * 3.6-30     2024-10-02 [1] CRAN (R 4.4.2)
 RColorBrewer  * 1.1-3      2022-04-03 [1] CRAN (R 4.4.0)
 Rcpp            1.0.13     2024-07-17 [1] CRAN (R 4.4.1)
 readr         * 2.1.5      2024-01-10 [1] CRAN (R 4.4.2)
 recipes       * 1.2.1      2025-03-25 [1] CRAN (R 4.4.3)
 remotes         2.5.0      2024-03-17 [1] CRAN (R 4.4.2)
 reshape         0.8.9      2022-04-12 [1] CRAN (R 4.4.2)
 reshape2      * 1.4.4      2020-04-09 [1] CRAN (R 4.4.3)
 rlang           1.1.5      2025-01-17 [1] CRAN (R 4.4.3)
 rpart           4.1.23     2023-12-05 [2] CRAN (R 4.4.2)
 rsample       * 1.2.1      2024-03-25 [1] CRAN (R 4.4.3)
 Rsolnp          1.16       2015-12-28 [1] CRAN (R 4.4.3)
 rstudioapi      0.17.1     2024-10-22 [1] CRAN (R 4.4.2)
 satellite       1.0.5      2024-02-10 [1] CRAN (R 4.4.3)
 scales        * 1.3.0      2023-11-28 [1] CRAN (R 4.4.1)
 sessioninfo     1.2.3      2025-02-05 [1] CRAN (R 4.4.2)
 sf            * 1.0-19     2024-11-05 [1] CRAN (R 4.4.2)
 shiny           1.10.0     2024-12-14 [1] CRAN (R 4.4.2)
 sm            * 2.2-6.0    2024-02-17 [1] CRAN (R 4.4.3)
 snakecase       0.11.1     2023-08-27 [1] CRAN (R 4.4.3)
 sp            * 2.1-4      2024-04-30 [1] CRAN (R 4.4.2)
 spacetime       1.3-2      2024-09-04 [1] CRAN (R 4.4.2)
 spam            2.11-1     2025-01-20 [1] CRAN (R 4.4.2)
 SparseM         1.84-2     2024-07-17 [1] CRAN (R 4.4.2)
 spatialsample * 0.6.0      2024-10-02 [1] CRAN (R 4.4.3)
 stars           0.6-7      2024-11-07 [1] CRAN (R 4.4.2)
 statmod         1.5.0      2023-01-06 [1] CRAN (R 4.4.2)
 stringi         1.8.4      2024-05-06 [1] CRAN (R 4.4.0)
 stringr       * 1.5.1      2023-11-14 [1] CRAN (R 4.4.1)
 survival        3.7-0      2024-06-05 [2] CRAN (R 4.4.2)
 terra           1.8-29     2025-02-26 [1] CRAN (R 4.4.3)
 tibble        * 3.2.1      2023-03-20 [1] CRAN (R 4.4.1)
 tidymodels    * 1.3.0      2025-02-21 [1] CRAN (R 4.4.3)
 tidyr         * 1.3.1      2024-01-24 [1] CRAN (R 4.4.1)
 tidysdm       * 1.0.0      2025-03-05 [1] CRAN (R 4.4.3)
 tidyselect      1.2.1      2024-03-11 [1] CRAN (R 4.4.1)
 tidyverse     * 2.0.0      2023-02-22 [1] CRAN (R 4.4.2)
 timechange      0.3.0      2024-01-18 [1] CRAN (R 4.4.2)
 timeDate        4041.110   2024-09-22 [1] CRAN (R 4.4.3)
 truncnorm       1.0-9      2023-03-20 [1] CRAN (R 4.4.2)
 tune          * 1.3.0      2025-02-21 [1] CRAN (R 4.4.3)
 tzdb            0.5.0      2025-03-15 [1] CRAN (R 4.4.3)
 units           0.8-5      2023-11-28 [1] CRAN (R 4.4.2)
 urlchecker      1.0.1      2021-11-30 [1] CRAN (R 4.4.2)
 usethis         3.1.0      2024-11-26 [1] CRAN (R 4.4.2)
 uuid            1.2-1      2024-07-29 [1] CRAN (R 4.4.1)
 vctrs           0.6.5      2023-12-01 [1] CRAN (R 4.4.1)
 vioplot       * 0.5.1      2025-02-23 [1] CRAN (R 4.4.3)
 viridisLite     0.4.2      2023-05-02 [1] CRAN (R 4.4.1)
 withr           3.0.2      2024-10-28 [1] CRAN (R 4.4.2)
 workflows     * 1.2.0      2025-02-19 [1] CRAN (R 4.4.3)
 workflowsets  * 1.1.0      2024-03-21 [1] CRAN (R 4.4.3)
 xtable          1.8-4      2019-04-21 [1] CRAN (R 4.4.2)
 xts             0.14.1     2024-10-15 [1] CRAN (R 4.4.2)
 yardstick     * 1.3.2      2025-01-22 [1] CRAN (R 4.4.3)
 zoo           * 1.8-12     2023-04-13 [1] CRAN (R 4.4.2)

 [1] C:/Users/alexw/AppData/Local/R/win-library/4.4
 [2] C:/Program Files/R/R-4.4.2/library
 * -- Packages attached to the search path.

-----------------------------------------------------------------------------------------------------------------------
```
