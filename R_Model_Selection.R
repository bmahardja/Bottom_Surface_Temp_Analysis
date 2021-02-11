###################################################
#Script to conduct model selection-------------------------------
###################################################

library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)
library(ggpubr)

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
################# List of models to choose the proper k  for x and y coordinates  ############
######################################################################################################

#Create blank list to start with
model_k_spatial <-list( )

#Start from k=10 and increase by 5 until 50
model_k_spatial[[1]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(10)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[2]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(15)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[3]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(20)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[4]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(25)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[5]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(30)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[6]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(35)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[7]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(40)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[8]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(45)),
                  data = temp_dataset, method="fREML", nthreads=3)
model_k_spatial[[9]]<-bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(50)),
                  data = temp_dataset, method="fREML", nthreads=3)


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
#Looks like 35 might be the best number

#Save figure to results folder
tiff(filename=file.path(results_root,"Figure_spatial_component_k_selection.tiff"), 
     height=10, width=20, units="in", pointsize=14, res=300, type="cairo",compression="lzw")
ggarrange(plot_k_spatial_select_01, plot_k_spatial_select_02, ncol=1, nrow=2)
dev.off()

#Free up some space
remove(model_k_spatial)

######################################################################################################
################# List of models to choose the proper k  for temperature anomaly covariate  ############
######################################################################################################
model_k_temp_anomaly <-list( )

model_k_temp_anomaly[[1]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,3)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[2]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,4)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[3]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,5)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[4]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,6)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[5]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,7)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[6]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,8)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[7]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,9)),
                               data = temp_dataset, method="fREML", nthreads=3)
model_k_temp_anomaly[[8]]<-bam(Temperature_difference ~  te(x,y,Temperature_anomaly_spatial, d=c(2,1), bs=c("tp","tp"), k=c(35,10)),
                               data = temp_dataset, method="fREML", nthreads=3)

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
#Looks like 7 might be the best number

#Save figure to results folder
tiff(filename=file.path(results_root,"Figure_temp_anomaly_k_selection.tiff"), 
     height=10, width=20, units="in", pointsize=14, res=300, type="cairo",compression="lzw")
ggarrange(plot_k_temp_anomaly_select_01, plot_k_temp_anomaly_select_02, ncol=1, nrow=2)
dev.off()

#Free up some space
remove(model_k_temp_anomaly)

#NOTE: k for julian day was already determined in Data setup file

######################################################################################################
################# Model run with just thin spline and cyclic cubic spline ############
######################################################################################################

model_01_thinspline_xy <- bam(Temperature_difference ~  te(x,y, d=c(2), bs=c("tp"), k=c(35)),
                     data = temp_dataset, method="fREML", nthreads=3)
summary(model_01_thinspline_xy)


model_01_thinspline_xy_jd <- bam(Temperature_difference ~  te(x,y,Julian_day_s, d=c(2,1), bs=c("tp","cc"), k=c(35,5)),
                                 data = temp_dataset, method="fREML", nthreads=3)
#gam.check(model_01_thinspline_xy_jd)
summary(model_01_thinspline_xy_jd)

