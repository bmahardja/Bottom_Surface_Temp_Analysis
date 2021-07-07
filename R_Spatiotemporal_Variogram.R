###################################################
#Script to create spatiotemporal variogram per Sam Bashevkin's code to evaluate autocorrelation-
###################################################

library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)
library(sf)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(plotly)

require(spacetime)
require(sp)
require(gstat)
require(itsadug)
require(patchwork)

data_root<-file.path("data-raw")
results_root<-file.path("results")

######################################################################################################
################# Load prepared datasets ############
######################################################################################################

#To get UTM data
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline_NoSSC_UTM10.shp"))
#EPSG 4269

#Dataset
temp_dataset <- readRDS("temperature_dataset.Rds")
temp_dataset$Date<-as.Date(temp_dataset$Date)

#Read saved model results
model_05_soapfilm_xy_jd_ta<-readRDS(file.path(results_root,"Model_05_soapfilm.Rds"))

######################################################################################################
#Function to create variogram from Sam
######################################################################################################

Data_vario<-temp_dataset%>%
  mutate(Resid=resid_gam(model_05_soapfilm_xy_jd_ta))

Data_coords<-Data_vario%>%
  st_as_sf(coords=c("x", "y"), crs=4269)%>%
  st_coordinates()%>%
  as_tibble()#%>%
  #mutate(across(c(X,Y), ~(.x-mean(.x))/1000)) #From Sam. I believe this was to convert to kilometers

#Create data frame for the variogram analysis
Data_vario<-bind_cols(Data_vario%>%
                        select(Date, Resid), Data_coords)

sp<-SpatialPoints(coords=data.frame(X=Data_vario$X, Y=Data_vario$Y))

sp2<-STIDF(sp, time=Data_vario$Date, 
           data=data.frame(Residuals=Data_vario$Resid))

#Calculate gamma
vario_01<-variogramST(Residuals~1, data=sp2, tunit="days", cores=2, tlags=0:14)
vario_02<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=2, tlags=1:8)

#Adjust the 0.5
vario_plot_01<-vario_01%>%
  mutate(spacelag_km=spacelag/1000,timelag_day=timelag-0.5)

vario_plot_02<-vario_02%>%
  mutate(spacelag_km=spacelag/1000,timelag_week=timelag-0.5)

p<-plot_ly(data=vario_plot_01,x=vario_plot_01$spacelag_km,y=vario_plot_01$timelag_day,z=vario_plot_01$gamma,type="mesh3d",
           intensity=vario_plot_01$gamma, colorscale="Viridis")
p<-p %>% layout(
  scene = list(
    xaxis = list(title = "Distance (km)"),
    yaxis = list(title = "Time difference (days)"),
    zaxis = list(title = "Gamma")
  ))
p


p2<-plot_ly(data=vario_plot_02,x=vario_plot_02$spacelag_km,y=vario_plot_02$timelag_week,z=vario_plot_02$gamma,type="mesh3d",
           intensity=vario_plot_02$gamma, colorscale="Viridis")
p2<-p2 %>% layout(
  scene = list(
    xaxis = list(title = "Distance (km)"),
    yaxis = list(title = "Time difference (weeks)"),
    zaxis = list(title = "Gamma")
  ))
p2






###############Don't use
######################################################################################################

#If necessary to do what Sam did
p_time<-ggplot(vario_plot, aes(x=timelag, y=gamma, color=spacelag, group=spacelag))+
  geom_line()+
  geom_point()+
  scale_color_viridis_c(name="Distance (km)")+
  #scale_x_continuous(breaks=c(2,4,6,8,10))+
  xlab("Time difference (months)")+
  theme_bw()+
  theme(legend.justification = "left")

p_time

p_space<-ggplot(vario_plot, aes(x=spacelag, y=gamma, color=timelag, group=timelag))+
  geom_line()+
  geom_point()+
  scale_color_viridis_c(breaks=c(2,4,6,8,10), name="Time difference\n(months)")+
  xlab("Distance (km)")+
  theme_bw()+
  theme(legend.justification = "left")

p_space

#####
ST_variogram<-function(model, data, cores){
  norm_resids<-resid_gam(model)
  
  Data_vario<-data%>%
    mutate(Resid=norm_resids)
  
  Data_coords<-Data_vario%>%
    st_as_sf(coords=c("x", "y"), crs=4269)%>%
    st_coordinates()%>%
    as_tibble()#%>%
    #mutate(across(c(X,Y), ~(.x-mean(.x))/1000))
  
  Data_vario<-bind_cols(Data_vario%>%
                          select(Date, Resid), Data_coords)
  sp<-SpatialPoints(coords=data.frame(X=Data_vario$X, Y=Data_vario$Y))
  sp2<-STIDF(sp, time=Data_vario$Date, 
             data=data.frame(Residuals=Data_vario$Resid))
  vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=cores, tlags=1:8)
  
  vario_plot<-vario%>%
    mutate(monthlag=as.integer(as.factor(timelag))-0.5)
  
  p_time<-ggplot(vario_plot, aes(x=monthlag, y=gamma, color=spacelag, group=spacelag))+
    geom_line()+
    geom_point()+
    scale_color_viridis_c(name="Distance (km)")+
    scale_x_continuous(breaks=c(2,4,6,8,10))+
    xlab("Time difference (months)")+
    theme_bw()+
    theme(legend.justification = "left")
  
  p_space<-ggplot(vario_plot, aes(x=spacelag, y=gamma, color=monthlag, group=monthlag))+
    geom_line()+
    geom_point()+
    scale_color_viridis_c(breaks=c(2,4,6,8,10), name="Time difference\n(months)")+
    xlab("Distance (km)")+
    theme_bw()+
    theme(legend.justification = "left")
  
  p_variogram<-p_time/p_space+plot_annotation(tag_levels="A")
  
  return(p_variogram)
}

###############################################################################
#Evaluate autocorrelation

model_final_variogram<-ST_variogram(model_05_soapfilm_xy_jd_ta, temp_dataset, 3)

ggsave(model_final_variogram, file.path(results_root,"Model_SR_Variogram.png"),
       device="tiff", width=8, height=5, units="in", dpi=600,compression="lzw")


ggsave(file="Figure_BayesianModel_Diagonal_arrow.tiff", plot=diagonal_plot_arrows, device=tiff(),
       path=output_root, width=18, height=15, units="in", dpi=600,compression="lzw")

tiff(filename=file.path(output_root,"Figure_BayesianModel_Diagonal_arrow.tiff"),
     type="cairo",
     units="in", 
     width=6, #10*1, 
     height=5, #22*1, 
     pointsize=3, #12, 
     res=600,
     compression="lzw")
print(diagonal_plot_arrows)
dev.off()
