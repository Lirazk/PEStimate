# PEStimate

The source code for the PEStimate calculator, an updated version of SelectionCalc.
The different risk functions are implemented in EmbryoSelection.R, and can be used directly.

ui.R and server.R is where the shiny app is defined, and for_graph.R is used to create the graph in the supplementary material in the (upcoming) paper.

## Running

To run the calculator on a local machine, install the package and run:

```{r}
# Install the package by uncommenting next line:
# devtools::install_github("Lirazk/PEStimate")
# Run the calculator
PEStimate:::shiny_calculator()
```
