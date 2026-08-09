# ============================================================
# TCR_Explore - Global setup
# ============================================================

# ---- Core packages ----

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(DT)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)


# ---- Global options ----

options(
  shiny.maxRequestSize = 100 * 1024^2
)

# ---- Source shared functions ----

source("R/functions/helper_functions.R")
source("R/functions/data_processing.R")
source("R/functions/tcr_functions.R")
source("R/functions/plot_functions.R")
