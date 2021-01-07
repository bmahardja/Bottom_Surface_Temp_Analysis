library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)
library(cvTools)

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


######################################################################################################
################# Build a set of candidate model list and conduct Cross-Validation ############
######################################################################################################



#Leave out the soap-film smooth for now since it takes longer to run
Cand.set.Temp.model <- list( )

Cand.set.Temp.model[[1]] <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5)),
                                 data = temp_dataset, method="fREML", discrete=T, nthreads=4)

Cand.set.Temp.model[[2]]<- bam(Temperature_difference ~  te(x,y,Temperature_anomaly, d=c(2,1), bs=c("tp","tp"), k=c(40,7)),
                                 data = temp_dataset, method="fREML", discrete=T, nthreads=4)

Cand.set.Temp.model[[3]] <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(40,5,7)),
                                    data = temp_dataset, method="fREML", discrete=T, nthreads=4)









#Cross validation test ---------------------------------------------------------------
#with K=5
k <- 5 #the number of folds

set.seed(2020)
folds <- cvFolds(NROW(Data_subset), K=k)

Data_subset$holdoutpred <- rep(0,nrow(Data_subset))

for(i in 1:k){
  train <- Data_subset[folds$subsets[folds$which != i], ] #Set the training set
  validation <- Data_subset[folds$subsets[folds$which == i], ] #Set the validation set
  
  newlm <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("tp", "tp","cc"), k=c(40, 25, 7, 7)) ,
               data = train, method="fREML", discrete=T, nthreads=4) #Get your new linear model (just fit on the train data)
  newpred <- predict(newlm,newdata=validation) #Get the predicitons for the validation set (from the model just fit on the train data)
  
  Data_subset[folds$subsets[folds$which == i], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
}

Data_subset$holdoutpred_CorrectSign<-ifelse(Data_subset$Temperature_difference==0,NA,(ifelse(Data_subset$Temperature_difference>0&Data_subset$holdoutpred>0|Data_subset$Temperature_difference<0&Data_subset$holdoutpred<0,1,0)))
