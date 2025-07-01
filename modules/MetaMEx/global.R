options(shiny.sanitize.errors=F) # sanitize errors

# libraries
library(tidyverse)
library(feather)
library(shinyjs)
library(shinycssloaders)
library(rmarkdown)

library(forestplot)
library(metafor)

library(httr)
library(jsonlite)

# functions for server
source("server/data_loading.R")
source("server/functions_human.R")
source("server/functions_mouse.R")

# tabs for UI
source("ui/tab_home.R")
source("ui/tab_metaanalysis_human.R")
source("ui/tab_metaanalysis_mouse.R")
source("ui/tab_methods.R")
