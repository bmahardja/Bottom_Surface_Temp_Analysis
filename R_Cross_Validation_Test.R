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

#Cross validation test ---------------------------------------------------------------

#Number of models to be tested
N_model<-5

#with K=10
k <- 10 #the number of folds
set.seed(2020)
folds <- cvFolds(NROW(temp_dataset), K=k)

#Set up data list
Cross_Validation_Results=list()

for(i in 1:N_model){
  Cross_Validation_Results[[i]]<-temp_dataset
  Cross_Validation_Results[[i]]$holdoutpred <- rep(0,nrow(temp_dataset))
}

########## MODEL 01
for(i in 1){
  for(j in 1:k){
  train <- Cross_Validation_Results[[i]][folds$subsets[folds$which != j], ] #Set the training set
  validation <-Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ] #Set the validation set
  
  new_model <- bam(Temperature_difference ~ te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5)),
               data = temp_dataset, method="fREML", discrete=T, nthreads=3) 
  newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
  
  Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
  message(paste0("Finished run ", j, "/",k))
  }
}


########## MODEL 02

for(i in 2){
  for(j in 1:k){
    train <- Cross_Validation_Results[[i]][folds$subsets[folds$which != j], ] #Set the training set
    validation <-Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ] #Set the validation set
    
    new_model <- bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(40,7)),
                     data = temp_dataset, method="fREML", discrete=T, nthreads=3)
    newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
    
    Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
    message(paste0("Finished run ", j, "/",k))
  }
}


########## MODEL 03

for(i in 3){
  for(j in 1:k){
    train <- Cross_Validation_Results[[i]][folds$subsets[folds$which != j], ] #Set the training set
    validation <-Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ] #Set the validation set
    
    new_model <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(40,5,7)),
                    data = temp_dataset, method="fREML", discrete=T, nthreads=3)
    newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
    
    Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
    message(paste0("Finished run ", j, "/",k))
  }
}

########## MODEL 04

for(i in 4){
  for(j in 1:k){
    train <- Cross_Validation_Results[[i]][folds$subsets[folds$which != j], ] #Set the training set
    validation <-Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ] #Set the validation set
    
    new_model <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5), by=WaterYear) + WaterYear,
                     data = temp_dataset, method="fREML", discrete=T, nthreads=3)
    newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
    
    Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
    message(paste0("Finished run ", j, "/",k))
  }
}

########## MODEL 05 - Soap-film smoother

for(i in 5){
  for(j in 1:k){
    train <- Cross_Validation_Results[[i]][folds$subsets[folds$which != j], ] #Set the training set
    validation <-Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ] #Set the validation set
    
    new_model <- bam(Temperature_difference ~  te(x, y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sf", "cc","tp"), k=c(40,5,7),xt = list(list(bnd = border.aut,nmax=1500),NULL,NULL))+
                       te(x, y, Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sw", "tp","cc"), k=c(40,5,7),xt = list(list(bnd = border.aut,nmax=1500),NULL,NULL)),
                     data = temp_dataset, method="fREML", discrete=T, nthreads=3, knots =knots_grid)
    newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
    
    Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
    message(paste0("Finished run ", j, "/",k))
  }
}

#Save cross validation results
saveRDS(Cross_Validation_Results, file="Cross_Validation_Results.Rds")

#Read saved cross validation results
Cross_Validation_Results<-readRDS("Cross_Validation_Results.Rds")

str(folds)

#Create a new column for kfold number
for(i in 1:N_model){
    Cross_Validation_Results[[i]]$kfold<-0
}

#Add kfold number to the kfold column
for(i in 1:N_model){
  for(j in 1:k){
   Cross_Validation_Results[[i]][folds$subsets[folds$which == j], "kfold"]<-j
  }
}

#Change column to factor
for(i in 1:N_model){
  Cross_Validation_Results[[i]]$kfold<-as.factor(Cross_Validation_Results[[i]]$kfold)
  Cross_Validation_Results[[i]]$resid_CV<-(Cross_Validation_Results[[i]]$Temperature_difference-Cross_Validation_Results[[i]]$holdoutpred)^2
}

###Create list for RMSE
RMSE_results=list()
for(i in 1:N_model){
  RMSE_results[[i]]<- Cross_Validation_Results[[i]] %>% group_by(kfold) %>% summarise(RMSE=sqrt(sum(resid_CV)/n()))
}


#Create calculation for predicting sign of bottom-surface temperature difference
for(i in 1:N_model){
  Cross_Validation_Results[[i]]$holdoutpred_CorrectSign<-ifelse(Cross_Validation_Results[[i]]$Temperature_difference==0,NA,(ifelse(Cross_Validation_Results[[i]]$Temperature_difference>0&Cross_Validation_Results[[i]]$holdoutpred>0|Cross_Validation_Results[[i]]$Temperature_difference<0&Cross_Validation_Results[[i]]$holdoutpred<0,1,0)))
}

#Calculate overall RMSE and Pearson correlation for each model
#https://stats.stackexchange.com/questions/85507/what-is-the-rmse-of-k-fold-cross-validation

Cross_Validation_Summary<-data.frame(ModelNumber=1:N_model,RMSE=0,Pearson=0)

for(i in 1:N_model){
  Cross_Validation_Summary$RMSE[i]<- sqrt(sum((RMSE_results[[i]]$RMSE)^2)/k)
  Cross_Validation_Summary$Pearson[i]<-cor(Cross_Validation_Results[[i]]$holdoutpred, Cross_Validation_Results[[i]]$Temperature_difference, method="pearson")
  Cross_Validation_Summary$CorrectSign[i]<-mean(Cross_Validation_Results[[i]]$holdoutpred_CorrectSign,na.rm=T)
}

Cross_Validation_Summary

write.csv(Cross_Validation_Summary,"Cross_Validation_Results.csv",row.names = F)



#####################################################Extra code

#Leave out the soap-film smooth for now since it takes longer to run

Cand.set.Temp.model <- list( )

Cand.set.Temp.model[[1]] <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(40,5)),
                                data = temp_dataset, method="fREML", discrete=T, nthreads=3)

Cand.set.Temp.model[[2]]<- bam(Temperature_difference ~  te(x,y,Temperature_anomaly, d=c(2,1), bs=c("tp","tp"), k=c(40,7)),
                               data = temp_dataset, method="fREML", discrete=T, nthreads=3)

Cand.set.Temp.model[[3]] <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(40,5,7)),
                                data = temp_dataset, method="fREML", discrete=T, nthreads=3)





