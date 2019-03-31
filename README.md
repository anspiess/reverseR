## reverseR: Linear Regression Stability to Significance Reversal (R-package)

*reverseR* is an R package that tests linear regressions for significance reversal through leave-one(multiple)-out.

#### Installation
*reverseR* is not yet available [on CRAN](https://cran.r-project.org/). However, you can install the latest development version of the code using the [devtools](https://cran.r-project.org/package=devtools) R package.

```R
# Install devtools, if you haven't already.
if (!'devtools' %in% installed.packages()) install.packages(devtools)
# Fetech the latest source code from github
devtools::install_github("anspiess/reverseR")
source("https://install-github.me/anspiess/reverseR")
```

# Running the GUI

The reverseR package has a fully functional shiny GUI which covers all 
functionality of the the package. To invoke the GUI within dedicated R IDEs 
(e.g., RStudio, RKWard) or a browser run the following line from the R command
line.

```R
reverseR::shinyInfl()
```
