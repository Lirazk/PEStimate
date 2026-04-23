# PEStimate: Predicting offspring disease risk after Polygenic Embryo Screening

PEStimate is a R package designed to estimate the relative and absolute risk reduction associated with Polygenic embryo screening (PES).

It utilizes the Liability threshold model to predict outcomes based on the disease prevalence, PRS accuracy and heritability.

* R/EmbryoSelection.R include the functions for the various risk reduction estimates.

* R/ui.R and R/server.R implements the Shiny application.

* for_graph.R is the script used to generate the figures for the publication (and a few which we didn't include).

## Installation

You can install the development version of PEStimate from GitHub using devtools:

```{r}
# If you don't have devtools installed:
# install.packages("devtools")

devtools::install_github("Lirazk/PEStimate")
```

## Quick start

1. Interactive Calculator

The easiest way to use PEStimate is through the built-in Shiny application:

```{r}
library(PEStimate)
PEStimate:::shiny_calculator()
```

2. Programmatic Usage

For batch processing or custom simulations, you can use the underlying risk functions directly. For example, to calculate risk reduction for the lowest-PRS strategy:

```{r}
risk_prediction_analytical(r2 = 0.1, K = 0.05, n = 5)
```

To estimate risk conditional on family history

```{r}
risk_prediction_exact(
  iter = 100000,
  n_embryos = 5, 
  r2 = 0.1, 
  h2 = 0.4, 
  K = 0.01,
  history = list(p1 = 1, sib_self = c(1, -1)) # Parent 1 affected, one sibling affected, one healthy
)
```

The same thing can be done with PRS. See ?risk_prediction_exact and ?risk_prediction_analytical for more information.