###################################################
#Script to conduct model selection-------------------------------
###################################################

library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)

data_root<-file.path("data-raw")
results_root<-file.path("results")

######################################################################################################
################# Load prepared datasets ############
######################################################################################################

#Border outline for soap-film smoother
border.aut <- readRDS("Soap_film_boundaries.Rds")

temp_dataset <- readRDS("temperature_dataset.Rds")

#Add knots
knots_grid <- read.csv("knots_grid.csv")

######################################################################################################
################# Model run with just thin spline and cyclic cubic spline ############
######################################################################################################


model_00_thinspline_xy <- bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(35)),
                     data = temp_dataset, method="fREML", nthreads=3)
summary(model_00_thinspline_xy)


model_01_thinspline_xy_jd <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(35,5)),
                                 data = temp_dataset, method="fREML", nthreads=3)
#gam.check(model_01_thinspline_xy_jd)
summary(model_01_thinspline_xy_jd)

