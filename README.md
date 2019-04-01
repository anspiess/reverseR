# reverseR: Linear Regression Stability to Significance Reversal (`R`-package)

## Overview
              
In biomedical literature, the most widely employed statistical procedure to analyze and visualize the association between two variables is linear regression. Data points that exert influence on the fit and its parameters are routinely, but not as often as required, identified by established influence measures and their corresponding cut-off values. In this `R` package, we specifically address the presence of influential data points that directly impact the statistical inference of the models, which none of the established measures detect, such as

*leverage*  
*dffits*  
*dfbeta(s)*  
*covratio*  
*Cook's distance*  
*studentized residuals*  

We call these data points **"reversers"**. `reverseR` tests linear regressions for significance reversal through leave-one(multiple)-out and checking if $\alpha \in [p_{\beta1}, p_{\beta1(i)}]$, where $p_{\beta1}$ is the *p*-value of the regression's slope and $p_{\beta1(i)}$ is the *p*-value with the *i*-th point deleted. This paradigm is along the lines of the living-in-oblivion measure *dfstat* or *dfstud* (Belsley, Kuh & Welsch, Regression diagnostics: Identifying influential data and sources of collinearity, 2004) that checks the impact of each response value on statistical inference.

## Repo Contents
- [R](./R): `R` package code.
- [man](./man): package manuals for the different functions.
- [inst/reverseR](./tests): files for the Shiny-based GUI.
              
## Hardware Requirements
The `reverseR` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 4 GB of RAM. For optimal performance, we recommend a computer with the following specs: RAM: 16+ GB; CPU: 4+ Cores, 3.3+ GHz/Core.

## Software Requirements
`R` version greater than 3.5.0 for Linux Ubuntu 16, Windows 7, 8, 10 or Mac.  
Several CRAN packages may be needed.

## Runtimes
Runtimes vary for the different functions:  
*lmInfl* for "reverser" analysis: ~ 1-5 sec.  
*lmMult* for multiple leave-out analysis: ~ 5-30 sec.  
*simInfl* for Monte Carlo simulation: 5-600 sec.  

## Package dependencies
Users should install the following packages prior to installing `reverseR`, from an `R` terminal:
```
install.packages(c("shiny", "markdown", "knitr"))
```
## Installation Guide
From an `R` session, type:
```
if (!'devtools' %in% installed.packages()) install.packages(devtools)
devtools::install_github("anspiess/reverseR")
source("https://install-github.me/anspiess/reverseR")
```

## Running the GUI

The reverseR package has a fully functional shiny GUI which covers all functionality of the the package. To invoke the GUI within dedicated R IDEs (e.g., RStudio, RKWard) or a browser run the following line from the R command line.
```
library(reverseR)
shinyInfl()
```