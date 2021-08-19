###################################################
#Script to conduct model construction and selection-------------------------------
###################################################

library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)
library(ggpubr)
library(cvTools)
library(parallel)
library(gridExtra)

data_root<-file.path("data-raw")
results_root<-file.path("results")

nc<-3 # number of cores
cl <- makeCluster(nc)

######################################################################################################
################# Load prepared datasets ############
######################################################################################################

#Border outline for soap-film smoother
border.aut <- readRDS("Soap_film_boundaries.Rds")

temp_dataset <- readRDS("temperature_dataset.Rds")

#Add knots
knots_grid <- read.csv("knots_grid.csv")

######################################################################################################
################# List of models to choose the proper k  for x and y coordinates  ############
######################################################################################################

#Create blank list to start with
model_k_spatial <-list( )

#Start from k=10 and increase by 5 until 50
model_k_spatial[[1]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(10)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[2]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(15)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[3]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(20)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[4]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(25)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[5]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(30)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[6]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(35)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[7]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(40)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[8]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(45)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_spatial[[9]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(50)),
                  data = temp_dataset, method="fREML", family=scat,cluster=cl)


#Extract adjusted R squared changes
#https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/summary.gam.html
k_test_results_Rsq <- vector("numeric", length(model_k_spatial))

for (i in 1:length(model_k_spatial)){
  k_test_results_Rsq[i]<- summary(model_k_spatial[[i]])$r.sq
}

#Extract AICc
k_test_results_AIC <- vector("numeric", length(model_k_spatial))

for (i in 1:length(model_k_spatial)){
  k_test_results_AIC[i]<- AICc(model_k_spatial[[i]])
}

#Store them in one data frame
model_k_spatial_results<-data.frame(k=c(seq(10,50,by=5)),AICc=k_test_results_AIC,R_squared=k_test_results_Rsq)

#Plot
plot_k_spatial_select_01<-ggplot() + 
  theme_bw() +
  geom_point(data=model_k_spatial_results, aes(x=k,y=R_squared),size=4) + 
  geom_line(data=model_k_spatial_results, aes(x=k,y=R_squared)) + 
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  scale_x_continuous(breaks=seq(10,50,5))+
  labs(x=NULL,title="GAM k evaluation for spatial component", y=expression(paste("Adj ",italic(R)^2,sep="")))
#plot_k_spatial_select_01
plot_k_spatial_select_02<-ggplot() + 
  theme_bw() +
  geom_point(data=model_k_spatial_results, aes(x=k,y=AICc),size=4) + 
  geom_line(data=model_k_spatial_results, aes(x=k,y=AICc)) + 
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  scale_x_continuous(breaks=seq(10,50,5))+
  labs(x="k", y="AICc")
#plot_k_spatial_select_02
########

#Save figure to results folder
tiff(filename=file.path(results_root,"Figure_spatial_component_k_selection.tiff"), 
     height=10, width=20, units="in", pointsize=14, res=300, type="cairo",compression="lzw")
ggarrange(plot_k_spatial_select_01, plot_k_spatial_select_02, ncol=1, nrow=2)
dev.off()

#Free up some space
remove(model_k_spatial)

#Looks like 20 might be the best number
#Set the k for spatial components
k_spatial_comp<-20

######################################################################################################
################# List of models to choose the proper k  for temperature anomaly covariate  ############
######################################################################################################
model_k_temp_anomaly <-list( )

model_k_temp_anomaly[[1]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,3)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[2]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,4)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[3]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,5)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[4]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,6)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[5]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,7)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[6]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,8)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[7]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,9)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)
model_k_temp_anomaly[[8]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(k_spatial_comp,10)),
                               data = temp_dataset, method="fREML", family=scat,cluster=cl)

#Extract adjusted R squared changes
k_test_results_Rsq <- vector("numeric", length(model_k_temp_anomaly))

for (i in 1:length(model_k_temp_anomaly)){
  k_test_results_Rsq[i]<- summary(model_k_temp_anomaly[[i]])$r.sq
}

#Extract AICc
k_test_results_AIC <- vector("numeric", length(model_k_temp_anomaly))

for (i in 1:length(model_k_temp_anomaly)){
  k_test_results_AIC[i]<- AICc(model_k_temp_anomaly[[i]])
}

#Store them in one data frame
model_k_temp_anomaly_results<-data.frame(k=c(3:10),AICc=k_test_results_AIC,R_squared=k_test_results_Rsq)

#Plot
plot_k_temp_anomaly_select_01<-ggplot() + 
  theme_bw() +
  geom_point(data=model_k_temp_anomaly_results, aes(x=k,y=R_squared),size=4) + 
  geom_line(data=model_k_temp_anomaly_results, aes(x=k,y=R_squared)) + 
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  scale_x_continuous(breaks=seq(3,10,1))+
  labs(x=NULL,title="GAM k evaluation for temperature anomaly covariate", y=expression(paste("Adj ",italic(R)^2,sep="")))
#plot_k_temp_anomaly_select_01
plot_k_temp_anomaly_select_02<-ggplot() + 
  theme_bw() +
  geom_point(data=model_k_temp_anomaly_results, aes(x=k,y=AICc),size=4) + 
  geom_line(data=model_k_temp_anomaly_results, aes(x=k,y=AICc)) + 
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  scale_x_continuous(breaks=seq(3,10,1))+
  labs(x="k", y="AICc")
#plot_k_temp_anomaly_select_02
########

#Save figure to results folder
tiff(filename=file.path(results_root,"Figure_temp_anomaly_k_selection.tiff"), 
     height=10, width=20, units="in", pointsize=14, res=300, type="cairo",compression="lzw")
ggarrange(plot_k_temp_anomaly_select_01, plot_k_temp_anomaly_select_02, ncol=1, nrow=2)
dev.off()

#Free up some space
remove(model_k_temp_anomaly)

#NOTE: k for julian day was already determined in Data setup file
#Looks like 3 might be the best number
#Set the k for temp anomaly
k_temp_anomaly<-3

######################################################################################################
################# Model runs for eventual extraction of AICc and R squared (and other metrics) ############
######################################################################################################

#01 Space only model
model_01_thinspline_xy <- bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(k_spatial_comp)),
                     data = temp_dataset, method="fREML", family=scat, cluster=cl)

#02 Space and season (julian day) model
model_02_thinspline_xy_jd <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(k_spatial_comp,5)),
                                 data = temp_dataset, method="fREML", family=scat, cluster=cl)

#03
#Space and season (julian day) and temperature anomaly (based on either interannual variability or time of day difference)
time1_conventional_model <- Sys.time() #Calculate how long it takes to run model
model_03_thinspline_xy_jd_ta <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(k_spatial_comp,5,k_temp_anomaly)),
                                    data = temp_dataset, method="fREML", family=scat, cluster=cl)

time2_conventional_model <- Sys.time() #Calculate how long it takes to run model
time2_conventional_model-time1_conventional_model
saveRDS(model_03_thinspline_xy_jd_ta, file.path(results_root,"Model_03_thinspline.Rds"))
#Roughly 6 minutes on 8/19/21

#04
#Space and season (julian day) and temperature anomaly + soap-film smoother to minimize bleedover from Suisun Marsh and Cache Slough Complex
time1_final_model <- Sys.time() #Calculate how long it takes to run model
model_04_soapfilm_xy_jd_ta <- bam(Temperature_difference ~  te(x, y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sf", "cc","tp"), k=c(k_spatial_comp,5,k_temp_anomaly),xt = list(list(bnd = border.aut,nmax=500),NULL,NULL))+
                                    te(x, y, Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sw","cc","tp"), k=c(k_spatial_comp,5,k_temp_anomaly),xt = list(list(bnd = border.aut,nmax=500),NULL,NULL)),
                                  data = temp_dataset, method="fREML", knots = knots_grid, family=scat, cluster=cl)
time2_final_model <- Sys.time()
time2_final_model-time1_final_model

#Generally ~13 mins at desktop, ~4 hours at work laptop (7/19/21)
#Model takes awhile to run, save it as Rds so it doesn't have to be repeated every time
saveRDS(model_04_soapfilm_xy_jd_ta, file.path(results_root,"Model_04_soapfilm.Rds"))

#Print out model diagnostics plot for best model
png(filename=file.path(results_root,"Model_diagnostics_scat.png"), units="in",type="cairo", bg="white", height=15, 
    width=15, res=300, pointsize=20)
par(mfrow=c(2,2)) 
gam.check(model_04_soapfilm_xy_jd_ta)
dev.off()


#######For skipping soap-film smoother model step
#Read saved model here
model_03_thinspline_xy_jd_ta<-readRDS(file.path(results_root,"Model_03_thinspline.Rds"))
model_04_soapfilm_xy_jd_ta<-readRDS(file.path(results_root,"Model_04_soapfilm.Rds"))


######################################################################################################
################# k-fold Cross Validation Runs ############
######################################################################################################

#Cross validation test ---------------------------------------------------------------

time1_cross_validation <- Sys.time() #To calculate how long it takes to run cross-validation

#Remove spatial component from dataset
temp_dataset$geometry<-NULL

#Number of models to be tested
N_model<-4

#with K=10
k <- 10 #the number of folds
set.seed(2021) #To replicate results if necessary
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
    
    new_model <- bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(k_spatial_comp)),
                     data = temp_dataset, method="fREML", family=scat, cluster=cl)
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
    
    new_model <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(k_spatial_comp,5)),
                     data = temp_dataset, method="fREML", family=scat, cluster=cl)
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
    
    new_model <- bam(Temperature_difference ~  te(x,y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("tp","cc", "tp"), k=c(k_spatial_comp,5,k_temp_anomaly)),
                     data = temp_dataset, method="fREML", family=scat, cluster=cl)
    newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
    
    Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
    message(paste0("Finished run ", j, "/",k))
  }
}

########## MODEL 04 - Soap-film smoother, final model

for(i in 4){
  for(j in 1:k){
    train <- Cross_Validation_Results[[i]][folds$subsets[folds$which != j], ] #Set the training set
    validation <-Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ] #Set the validation set
    
    new_model <- bam(Temperature_difference ~  te(x, y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sf", "cc","tp"), k=c(k_spatial_comp,5,k_temp_anomaly),xt = list(list(bnd = border.aut,nmax=500),NULL,NULL))+
                       te(x, y, Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sw","cc","tp"), k=c(k_spatial_comp,5,k_temp_anomaly),xt = list(list(bnd = border.aut,nmax=500),NULL,NULL)),
                     data = temp_dataset, method="fREML", knots = knots_grid)
    newpred <- predict(new_model,newdata=validation) #Get the predictions for the validation set (from the model just fit on the train data)
    
    Cross_Validation_Results[[i]][folds$subsets[folds$which == j], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
    message(paste0("Finished run ", j, "/",k))
  }
}


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

#To calculate how long it takes to run cross-validation
time2_cross_validation <- Sys.time() 
time2_cross_validation-time1_cross_validation

#Save cross validation results
saveRDS(Cross_Validation_Results, file="Cross_Validation_Results.Rds")

#Roughly 1.2 hours for cross validation on 8/19/21
###########Start here to skip long cross-validation process

#Read saved cross validation results
Cross_Validation_Results<-readRDS("Cross_Validation_Results.Rds")


###Create list for RMSE
RMSE_results=list()
for(i in 1:N_model){
  RMSE_results[[i]]<- Cross_Validation_Results[[i]] %>% group_by(kfold) %>% summarise(RMSE=sqrt(sum(resid_CV)/n()))
}


#Create calculation for predicting sign of bottom-surface temperature difference
for(i in 1:N_model){
  Cross_Validation_Results[[i]]$holdoutpred_CorrectSign<-ifelse(Cross_Validation_Results[[i]]$Temperature_difference==0,NA,(ifelse(Cross_Validation_Results[[i]]$Temperature_difference>0&Cross_Validation_Results[[i]]$holdoutpred>0|Cross_Validation_Results[[i]]$Temperature_difference<0&Cross_Validation_Results[[i]]$holdoutpred<0,1,0)))
}

#Create data frame for model selection summary table
Model_Selection_Summary<-data.frame(ModelNumber=1:N_model,RMSE=0,Pearson=0,Rsquared=0,AICc=0)

#Calculate overall RMSE and Pearson correlation for each model
#https://stats.stackexchange.com/questions/85507/what-is-the-rmse-of-k-fold-cross-validation

for(i in 1:N_model){
  Model_Selection_Summary$RMSE[i]<- sqrt(sum((RMSE_results[[i]]$RMSE)^2)/k)
  Model_Selection_Summary$Pearson[i]<-cor(Cross_Validation_Results[[i]]$holdoutpred, Cross_Validation_Results[[i]]$Temperature_difference, method="pearson")
  Model_Selection_Summary$CorrectSign[i]<-mean(Cross_Validation_Results[[i]]$holdoutpred_CorrectSign,na.rm=T)
}


#Add R squared to table
Model_Selection_Summary$Rsquared[1]<- summary(model_01_thinspline_xy)$r.sq
Model_Selection_Summary$Rsquared[2]<- summary(model_02_thinspline_xy_jd)$r.sq
Model_Selection_Summary$Rsquared[3]<- summary(model_03_thinspline_xy_jd_ta)$r.sq
Model_Selection_Summary$Rsquared[4]<- summary(model_04_soapfilm_xy_jd_ta)$r.sq

#Add R squared to table
Model_Selection_Summary$AICc[1]<- AICc(model_01_thinspline_xy)
Model_Selection_Summary$AICc[2]<- AICc(model_02_thinspline_xy_jd)
Model_Selection_Summary$AICc[3]<- AICc(model_03_thinspline_xy_jd_ta)
Model_Selection_Summary$AICc[4]<- AICc(model_04_soapfilm_xy_jd_ta)

write.csv(Model_Selection_Summary,file.path(results_root,paste("Model_Selection_Results",Sys.Date(),".csv",sep="")),row.names = F)

