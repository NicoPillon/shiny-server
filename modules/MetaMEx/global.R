options(shiny.sanitize.errors=F) # sanitize errors

# libraries
library(shinyjs)
library(shinycssloaders)
library(rmarkdown)

library(ggplot2)
library(gplots)
library(ggpubr)
library(ggprism)
library(ggfortify)

library(dplyr)
library(DT)
library(stringr)
library(scales)
library(rvest)
library(feather)

library(forestplot)
library(metafor)

# functions for server
source("server/data_loading.R")
source("server/functions_human.R")
source("server/functions_mouse.R")

# tabs for UI
source("ui/tab_home.R")
source("ui/tab_metaanalysis_human.R")
source("ui/tab_metaanalysis_mouse.R")
source("ui/tab_datasets.R")
source("ui/tab_download.R")
