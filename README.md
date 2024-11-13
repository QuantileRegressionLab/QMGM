# README file

The scripts in this repository can be used to reproduce the results in *Mid-quantile mixed graphical models with an application to mass public shootings in the United States* by L. Merlo, M. Geraci and L. Petrella.

## Prerequisites
### Software requirements

-   [R](https://cran.r-project.org/) version 4.4.1 or higher
-   [RStudio](https://rstudio.com/) version 2024.09.0+375 or higher
-   [Qtools](https://github.com/marco-geraci/Qtools) version 1.5.9 or higher

## Dataset
The post-processed dataset is included in the excel file `Mass_Shootings.xlsx` in the `data` folder. The raw data (version 7 - Updated June 2023) can be freely downloaded at [www.theviolenceproject.org/mass-shooter-database](https://www.theviolenceproject.org/mass-shooter-database/).The codebook can be found in the main body of the text and more detailed information can be found at [www.theviolenceproject.org](https://www.theviolenceproject.org/).

## Script description
### WorkHorse.R
This script contains the functions to run the `R` scripts in this repository. The main functions are 
-    `QMGM` allows to fit the proposed mid-quantile mixed graphical model (QMGM).
-    `fit.cat` allows to estimate the regression parameters of a LASSO-penalized mid-quantile regression.
-    `midCDF_est` allows to estimate the conditional mid-CDF.


### Data_preparation_Mass_Shootings.R
Preliminary step to clean the data and provide a suitable format for the analysis. This script can be used to reproduce Table 1 and Figure 2 of the paper:
-  Table 1: Summary statistics of the variables in the sample.
-  Figure 2: Number of fatalities per year and number of shootings per year.


### Mass_Shootings_ave.R
This script can be used to reproduce Figures 3 and 4 of the paper.
-  Figure 3: Estimated network structures from QMGM7 and MGM.
-  Figure 4: Local centrality measures for the MGM and QMGM7 estimated networks.




### Running the Script
1. Open the `RealData.R` script in RStudio.
2. At line 32, upload the dataset `returns.RData` from the `data` folder.
3. At line 36, select the type of model (`expectile` or `quantile` for CEHMM and CQHMM respectively) to be used.
4. At line 48, set the number of cores to be used in the parallel computation.
5. From line 62 define the regressors, the dependent variable, the number of states, the number of restarts, the vector of quantiles or expectiles considered, and the copula to be used.


## Fig&Tab_out.R
This script contains the code to reproduce Table S14, Table S15, Figure S3, Figure S2 and Figure 1 of the paper. The script reads the RData files `returns.RData` and `price_xts.RData` from the `data` folder.

### Script Description
The script `Figures_out.R` generates the followings:

- Table S14: Descriptive statistics.
- Table S15: Empirical correlation matrix.
- Figure S3: QQ plots of standardized residuals.
- Figure S2: Cryptocurrencies daily normalized prices and log return series.
- Figure 1: Time series of returns for the five cryptocurrencies colored according to the two-state fitted models.
