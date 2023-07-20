
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ipolygrowth

<!-- badges: start -->
<!-- badges: end -->

The goal of ipolygrowth is to calculate bacterial growth curve
parameters using fourth degree polynomial functions. Functions are
available for a single biological sample or multiple samples.

## Installation

The package can be installed from CRAN with:

``` r
install.packages("ipolygrowth")
```

or the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/kivanvan/ipolygrowth", upgrade = F, quiet = T)
```

## Example

The example data comes from the
[growthrates](https://CRAN.R-project.org/package=growthrates) package.  
This is a basic example which shows you how to use the single sample
function:

``` r
library(ipolygrowth)
library(dplyr)

# example data comes from the growthrates package (available on CRAN)
data <- growthrates::bactgrowth

# subset data to a single biological sample
df.singlesample <- data %>% dplyr::filter(strain == "D", conc == 0)

# calculate growth curve parameters using ipolygrowth function
out.singlesample <- ipg_singlesample(data = df.singlesample, time.name = "time", y.name = "value")
```

The output is a list, including a table of growth parameter estimates,
the polynomial model, a table of beta coefficients, and a table of
fitted values. Growth parameters include peak growth rate, peak growth
time, doubling time (at the peak growth), lag time, max y, and max y
time. View the results by calling each list element like:

``` r
out.singlesample$estimates
#>   peak growth rate peak growth time doubling time  lag time     max y
#> 1      0.005474298         3.636922      126.6185 0.1231345 0.1073791
#>   max y time
#> 1         30
```

For more instructions, please refer to the vignette.
