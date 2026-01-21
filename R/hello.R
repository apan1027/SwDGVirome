# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}



install.packages("pak")
pak::pak("rujinlong/projthis@aff52fb")
library("projthis")
projthis::proj_use_workflow("analyses")
projthis::use_qmd("01-tse-construction")

projthis::use_qmd("05-fig5-heatmap.qmd")
projthis::use_qmd("06-figs1-host-composition")
projthis::use_qmd("07-figs2-venn-diagram")
projthis::use_qmd("08-figs3-lifestyle")



projthis::use_qmd("09-figs4-abundant-patterns")






#
renv::snapshot()

renv::snapshot()

renv::restore()
