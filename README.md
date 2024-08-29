
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

Load the packages first.

``` r
library(ipolygrowth)
library(dplyr)
```

The example data comes from the
[growthrates](https://CRAN.R-project.org/package=growthrates) package.

``` r
# example data comes from the growthrates package (available on CRAN)
if (!"growthrates" %in% installed.packages()) {install.packages("growthrates")}
data <- growthrates::bactgrowth
```

Alternatively, download the `bactgrowth.txt` from
[here](https://github.com/tpetzoldt/growthrates/tree/master/data) to the
directory of your script. The data can then be read using the following
code.

``` r
data <- read.table("bactgrowth.txt", header = TRUE) %>%
  mutate(strain = factor(strain, levels = c("D", "R", "T")))
```

This is a basic example which shows you how to use the single sample
function:

``` r
# subset data to a single biological sample
df.singlesample <- data %>% dplyr::filter(strain == "D", conc == 0)

# calculate growth curve parameters using ipolygrowth function
out.singlesample <- ipg_singlesample(data = df.singlesample, time.name = "time", y.name = "value")
#> max y time is equal to the largest value of "time"
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

For more instructions and the expected output of vignette, please refer
to the [vignette on
CRAN](https://CRAN.R-project.org/package=ipolygrowth).
