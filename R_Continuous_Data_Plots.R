
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
library(viridis)
require(spacetime)
require(sp)
require(gstat)
require(itsadug)
require(patchwork)

data_root<-file.path("data-raw")
results_root<-file.path("results")


contin_temp_dataset <- readRDS("depthdata_QC_20201005.Rds")
unique(contin_temp_dataset$Station)
#Only MRZ and RRI show stratification
str(contin_temp_dataset)

annual_summary <- contin_temp_dataset %>% filter(WaterCol=="Surface",Station %in% c("MRZ","RRI")) %>% 
  mutate(Date=as.Date(Datetime)) %>% group_by(Date,Station) %>% summarise(DailyTemp=mean(Temp)) %>% 
  mutate(Year=year(Date) )%>% group_by(Year,Station) %>% summarise(AnnualTemp=mean(DailyTemp)) 

annual_summary_MRZ <- annual_summary %>% filter(Station=="MRZ")
median(annual_summary_MRZ$AnnualTemp) #Median is between 2006 and 2001

annual_summary_RRI <- annual_summary %>% filter(Station=="RRI")
median(annual_summary_RRI$AnnualTemp) #Median is between 2008 and 2014

#Different median year for each station
#Just pick 2019 (most recent year) for now

contin_temp_dataset_subset<-contin_temp_dataset%>% filter(Station %in% c("MRZ","RRI")) %>% mutate(Date=as.Date(Datetime),Year=year(Datetime)) %>%
  filter(Year==2017) #Set year
contin_temp_dataset_subset$WaterCol<-as.factor(contin_temp_dataset_subset$WaterCol)

#Calculate percent of time temperature was suitable for Delta Smelt and Chinook Salmon
contin_temp_dataset_percent<- contin_temp_dataset_subset %>% mutate(PercentSuitable_DeltaSmelt=ifelse(Temp<25,1,0),PercentSuitable_ChinookSalmon=ifelse(Temp<20,1,0)) %>%
  group_by(Station,Date,Year,WaterCol) %>% summarise(PercentSuitable_ChinookSalmon=mean(PercentSuitable_ChinookSalmon),PercentSuitable_DeltaSmelt=mean(PercentSuitable_ChinookSalmon))

#Remove middle
contin_temp_dataset_percent<-contin_temp_dataset_percent %>% filter(WaterCol %in% c("Surface", "Bottom"))

#Create figure
plot_habitat_deltasmelt <-ggplot2::ggplot(contin_temp_dataset_percent,aes(y=WaterCol,x=Date,fill=PercentSuitable_DeltaSmelt))+
  ggplot2:: geom_tile() +scale_fill_viridis(name="Percent Suitable",option ="magma")+
  ggplot2::theme_bw()+
  ggplot2::facet_grid(Station~.)+
  ggplot2::guides(fill = guide_legend(title = "Proportion of time\n temperature was\n <25 C"))+
  ggplot2::labs(y = "Delta Smelt")+
  ggplot2::theme(plot.title=element_text(size=11), 
               axis.text.x=element_text(size=11, color="black"), 
               axis.text.y = element_text(size=10, color="black",angle=45), 
               axis.title.x = element_blank(), 
               axis.title.y = element_text(size=14, color="black"),
               strip.text = element_text(size = 11),
               legend.text=element_text(size = 11),
               strip.background = element_rect(size=0.3)) 
plot_habitat_deltasmelt

plot_habitat_chinooksalmon <-ggplot2::ggplot(contin_temp_dataset_percent,aes(y=WaterCol,x=Date,fill=PercentSuitable_ChinookSalmon))+
  ggplot2:: geom_tile() +scale_fill_viridis(name="Percent Suitable",option ="magma")+
  ggplot2::theme_bw()+
  ggplot2::facet_grid(Station~.)+
  ggplot2::guides(fill = guide_legend(title = "Proportion of time\n temperature was\n <20 C"))+
  ggplot2::labs(y = "Chinook Salmon")+
  ggplot2::theme(plot.title=element_text(size=9), 
                 axis.text.x=element_text(size=9, color="black"), 
                 axis.text.y = element_text(size=8, color="black",angle=45), 
                 axis.title.x = element_blank(), 
                 axis.title.y = element_text(size=12, color="black"),
                 strip.text = element_text(size = 9),
                 legend.text=element_text(size = 9),
                 strip.background = element_rect(size=0.3)) 
plot_habitat_chinooksalmon



#Print out smelt and salmon figure together
tiff(filename=file.path(results_root,"Figure_Continuous_Data_Suitability.tiff"), units="in",type="cairo", bg="white", height=5, 
     width=10, res=300, pointsize=10,compression="lzw")
ggarrange(plot_habitat_chinooksalmon, plot_habitat_deltasmelt, ncol=1, nrow=2,labels = c("A","B"),hjust=-0.1)
dev.off()