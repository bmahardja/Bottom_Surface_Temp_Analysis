library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)

source("soap_checker/soap_check.R")
data_root<-file.path("data-raw")
results_root<-file.path("results")

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline_NoSSC_UTM10.shp"))

###################################################
################# Set up Boundary and Knots ############
###################################################

Delta.xy.aut <- tidy(Delta.aut)
head(Delta.xy.aut)

deltacoords <- Delta.xy.aut %>% dplyr::select(long,lat,piece)
names(deltacoords) <- c("x", "y", "piece")
borderlist <- split(deltacoords, deltacoords$piece)
names(borderlist)

Delta.xy.aut <- Delta.xy.aut %>% rename(x = long, y = lat)
#Delta.xy.aut$piece

border.aut <- lapply(borderlist, "[", c(1,2))
nr <- seq(1,length(borderlist))

border.aut <- lapply(nr, function(n) as.list.data.frame(border.aut[[n]]))

#Load knots
knots_grid <- read.csv("knots_grid.csv")


######################################################################################################
################# Read dataset and create temperature anomaly calculation ############
######################################################################################################

temp_dataset<-read.csv("temperature_dataset.csv")
######################
#Try to use the reduced  dataset for now
#temp_dataset<-read.csv("temperature_dataset_edited.csv")

str(temp_dataset)

# Calculate difference of bottom temperature from surface temperature as response variable
temp_dataset$Temperature_difference <- temp_dataset$Temperature_bottom-temp_dataset$Temperature
hist(temp_dataset$Temperature_difference)

#Temperature Anomaly Model---------------------------------------------
##Julian day and temperature model
temperature_anomaly_GAM<- gam(Temperature ~ s(Julian_day_s,bs="cc",k=5),data=temp_dataset)
summary(temperature_anomaly_GAM)
gam.check(temperature_anomaly_GAM)
plot(temperature_anomaly_GAM)

#Add term to the dataset
temp_dataset$Temperature_prediction <-predict(temperature_anomaly_GAM,temp_dataset)
temp_dataset$Temperature_anomaly <-temp_dataset$Temperature - temp_dataset$Temperature_prediction 
hist(temp_dataset$Temperature_anomaly)

############Summary of data
summary(temp_dataset$Temperature_difference)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-5.7600 -0.3000 -0.1000 -0.1669  0.0000  6.5600 
sd(temp_dataset$Temperature_difference)
#0.4546607

#Summarize by region for reference
temp_dataset$SubRegion=as.factor(temp_dataset$SubRegion)
temp_dataset$SampleSize=1

data_by_subregion<-temp_dataset %>% group_by(SubRegion) %>% summarize(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_subregion_month<-temp_dataset %>% group_by(SubRegion,Month) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_subregion_year<-temp_dataset %>% group_by(SubRegion,Year) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_year<-temp_dataset %>% group_by(Year) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_month<-temp_dataset %>% group_by(Month) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

#By inspecting the summarized data by subregion, it's apparent that temperature anomaly varies by region.
#South Delta is always warmer than North Delta as an example


##Construct model with longitude and latitude on top of julian day
temperature_anomaly_GAM_spatial<- gam(Temperature ~ te(x,y,Julian_day_s, d=c(2,1) ,bs=c("tp","cc"),k=c(10,5)),data=temp_dataset)
summary(temperature_anomaly_GAM_spatial)
#R-sq.(adj) =  0.907   Deviance explained = 90.7%

#K index is too low at 0.72, and edf is pretty close to k', but the goal was to remove collinearity, not fit
#k'  edf k-index p-value    
#59.0 57.8    0.51  <2e-16 ***
gam.check(temperature_anomaly_GAM_spatial)
plot(temperature_anomaly_GAM_spatial)

#Add term to the dataset
temp_dataset$Temperature_prediction_spatial <-predict(temperature_anomaly_GAM_spatial,temp_dataset)
temp_dataset$Temperature_anomaly_spatial <-temp_dataset$Temperature - temp_dataset$Temperature_prediction_spatial 
hist(temp_dataset$Temperature_anomaly_spatial)

#Make water year factor

temp_dataset$WaterYear<-as.factor(temp_dataset$WaterYear)
######################################################################################################
################# Model run with just thin spline and cyclic cubic spline ############
######################################################################################################


#Run models with Julian day and "temperature anomaly"

model_01_thinspline_xy_jd <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(25,5)),
                              data = temp_dataset, method="fREML", nthreads=3)
#gam.check(model_01_thinspline_xy_jd)
summary(model_01_thinspline_xy_jd)
#R-sq.(adj) =  0.352   Deviance explained = 35.7%
AICc(model_01_thinspline_xy_jd)
#12134.72

model_02_thinspline_xy_ta <- bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(25,7)),
                                 data = temp_dataset, method="fREML", nthreads=3)
#gam.check(model_02_thinspline_xy_ta)
summary(model_02_thinspline_xy_ta)
#R-sq.(adj) =  0.422   Deviance explained = 42.7%
AICc(model_02_thinspline_xy_ta)
#11076.24

model_03_thinspline_xy_jd_ta <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(25,5,7)),
                              data = temp_dataset, method="fREML", nthreads=3)
#gam.check(model_03_thinspline_xy_jd_ta)
summary(model_03_thinspline_xy_jd_ta)
#R-sq.(adj) =  0.541   Deviance explained = 55.5%
AICc(model_03_thinspline_xy_jd_ta)
#9083.439

##Add year interaction model
model_04_thinspline_xy_jd_yr <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(25,5), by=WaterYear)+WaterYear,
                                 data = temp_dataset, method="fREML", nthreads=3)
summary(model_04_thinspline_xy_jd_yr)
#R-sq.(adj) =  0.352   Deviance explained = 35.7%
AICc(model_04_thinspline_xy_jd_yr)
#12135.25

######################################################################################################
################# Model run with soap film smooth ############
######################################################################################################

model_05_soapfilm_xy_jd_ta <- bam(Temperature_difference ~  te(x, y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sf", "cc","tp"), k=c(25,5,7),xt = list(list(bnd = border.aut,nmax=1000),NULL,NULL))+
                           te(x, y, Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sw","cc","tp"), k=c(25,5,7),xt = list(list(bnd = border.aut,nmax=1000),NULL,NULL)),
                         data = temp_dataset, method="fREML", nthreads=3, knots =knots_grid)

summary(model_05_soapfilm_xy_jd_ta)
gam.check(model_05_soapfilm_xy_jd_ta)
#R-sq.(adj) =  0.366   Deviance explained = 38.7%
AICc(model_05_soapfilm_xy_jd_ta)
#7166.115


######################################################################################################
################# List of models to choose k  ############
######################################################################################################

model_k <-list( )
#Start from k=10 and increase by 5 until 40
model_k[[1]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(10,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[2]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(15,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[3]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(20,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[4]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(25,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[5]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(30,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[6]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(35,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[7]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5)),
                data = temp_dataset, method="fREML", nthreads=3)
model_k[[8]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(45,5)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k[[9]]<-bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(50,5)),
                  data = temp_dataset, method="fREML", nthreads=3)

############
#Evaluate adjusted R squared changes
#https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/summary.gam.html
k_test_results_Rsq <- vector("numeric", 9L)

for (i in 1:9){
  k_test_results_Rsq[i]<- summary(model_k[[i]])$r.sq
}
plot(k_test_results_Rsq)
#Looks like 35 might be the best number

#Evaluate AICc changes
k_test_results_AIC <- vector("numeric", 9L)

for (i in 1:9){
  k_test_results_AIC[i]<- AICc(model_k[[i]])
}

plot(k_test_results_AIC)

#k=35 for x and y spatial components seem to be most balanced between fit and complexity based on AICc and adjusted R^2







############### Fit test model  ######################
#m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut,nmax=1000)),
#              data = temp_dataset, method = "fREML", discrete=T, nthreads=3, knots = knots_grid)

#Test with 5000 nmax and the knots_custom didn't work and gave an error
knots_edit<-knots_grid[-7,]


m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut,nmax=1500)),
              data = temp_dataset, method = "fREML", discrete=T, nthreads=3, knots = knots_grid)
