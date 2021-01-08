library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)

# No need to specify absolute file paths in an Rstudio project (if you open it as a project).
data_root<-file.path("data-raw")

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline3_UTM10.shp"))

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

#Use custom knots
knots_custom <- read.csv("knots_custom.csv")
knots_grid <- read.csv("knots_grid.csv")


######################################################################################################
################# Read dataset and create temperature anomaly calculation ############
######################################################################################################

temp_dataset<-read.csv("temperature_dataset.csv")
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
#-5.7600 -0.3000 -0.1000 -0.1775  0.0000  6.5600 
sd(temp_dataset$Temperature_difference)
#0.4779963

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
temperature_anomaly_GAM_spatial<- gam(Temperature ~ te(x,y,Julian_day_s, d=c(2,1) ,bs=c("tp","cc"),k=c(15,5)),data=temp_dataset)
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

#This model has been causing crashes for whatever reason
#model_01_thinspline_xy <- bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(40)),
#                       data = temp_dataset, method="fREML", discrete=T, nthreads=3)
#gam.check(model_01_thinspline_xy)
#summary(model_01_thinspline_xy)
#R-sq.(adj) =   0.15   Deviance explained = 22.7%

model_02_thinspline_xy_jd <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5)),
                              data = temp_dataset, method="fREML", discrete=T, nthreads=3)
#gam.check(model_02_thinspline_xy_jd)
#summary(model_02_thinspline_xy_jd)
#R-sq.(adj) =   0.22   Deviance explained = 22.7%

model_03_thinspline_xy_ta <- bam(Temperature_difference ~  te(x,y,Temperature_anomaly, d=c(2,1), bs=c("tp","tp"), k=c(40,7)),
                                 data = temp_dataset, method="fREML", discrete=T, nthreads=3)
#gam.check(model_03_thinspline_xy_ta)
#summary(model_03_thinspline_xy_ta)
#R-sq.(adj) =  0.315   Deviance explained = 32.4%

model_04_thinspline_xy_jd_ta <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(40,5,7)),
                              data = temp_dataset, method="fREML", discrete=T, nthreads=3)
#gam.check(model_04_thinspline_xy_jd_ta)
#summary(model_04_thinspline_xy_jd_ta)
#R-sq.(adj) =   0.41   Deviance explained =   43%

model_05_thinspline_xy_jd_tas <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(40,5,7)),
                                    data = temp_dataset, method="fREML", discrete=T, nthreads=3)
summary(model_05_thinspline_xy_jd_tas)
#gam.check(model_05_thinspline_xy_jd_tas)
#R-sq.(adj) =   0.41   Deviance explained =   43.1%
AICc(model_05_thinspline_xy_jd_tas)
#8382.725

##Add year interaction model
model_06_thinspline_xy_jd_yr <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5), by=WaterYear)+WaterYear,
                                 data = temp_dataset, method="fREML", discrete=T, nthreads=3)
summary(model_06_thinspline_xy_jd_yr)
#R-sq.(adj) =  0.268   Deviance explained = 29.7%
AICc(model_06_thinspline_xy_jd_yr)
#10548.09
plot(model_06_thinspline_xy_jd_yr)

######################################################################################################
################# Model run with soap film smooth ############
######################################################################################################


model_06_soapfilm_xy_jd_ta <- bam(Temperature_difference ~  te(x, y,Julian_day_s,Temperature_anomaly, d=c(2,1,1), bs=c("sf", "cc","tp"), k=c(40,5,7),xt = list(list(bnd = border.aut,nmax=1500),NULL,NULL))+
                           te(x, y, Julian_day_s,Temperature_anomaly, d=c(2,1,1), bs=c("sw", "tp","cc"), k=c(40,5,7),xt = list(list(bnd = border.aut,nmax=1500),NULL,NULL)),
                         data = temp_dataset, method="fREML", discrete=T, nthreads=3, knots =knots_grid)
summary(model_06_soapfilm_xy_jd_ta)
gam.check(model_06_soapfilm_xy_jd_ta)
#R-sq.(adj) =  0.413   Deviance explained = 43.3%
AICc(model_06_soapfilm_xy_jd_ta)
#8345.612



######################################################################################################
################# Model run with just thin spline and cyclic cubic spline ############
######################################################################################################

model_list<- as.character(c(model_01_thinspline_xy,model_05_thinspline_xy_jd_tas))

























############### Fit test model  ######################
#m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut,nmax=1000)),
#              data = temp_dataset, method = "fREML", discrete=T, nthreads=3, knots = knots_grid)
#^ this one worked
#Test with 5000 nmax and the knots_custom didn't work and gave an error
knots_edit<-knots_grid[-7,]


m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut,nmax=1000)),
              data = temp_dataset, method = "fREML", discrete=T, nthreads=3, knots = knots_edit)

knots_edit<-knots_edit[-6,]
knots_edit<-knots_edit[-9,]
knots_edit<-knots_edit[-10,]
knots_edit<-knots_edit[-14,]
knots_edit<-knots_edit[-14,]
knots_edit<-knots_edit[-17,]



plot(Delta.aut, col="grey")
points(knots, pch=21, bg="red")
points(knots_edit, pch=21, bg="yellow")


#text(knots, labels=rownames(knots))