library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

outlier_region <- 