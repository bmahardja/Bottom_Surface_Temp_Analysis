library(tidyverse)
library(lubridate)
library(sf)

temp_dataset <- readRDS("temperature_dataset.Rds")
temp_dataset <- st_drop_geometry(temp_dataset)

min(temp_dataset$Temperature)
max(temp_dataset$Temperature)
mean(temp_dataset$Temperature)
median(temp_dataset$Temperature)

min(temp_dataset$Temperature_bottom)
max(temp_dataset$Temperature_bottom)
mean(temp_dataset$Temperature_bottom)
median(temp_dataset$Temperature_bottom)

min(temp_dataset$Temperature_difference)
max(temp_dataset$Temperature_difference)
mean(temp_dataset$Temperature_difference)
median(temp_dataset$Temperature_difference)
sd(temp_dataset$Temperature_difference)

#Proportion of data that is 0
mean(ifelse(temp_dataset$Temperature_difference==0,1,0))

mid_sac_river<-temp_dataset %>% filter(SubRegion == "Middle Sacramento River")

mid_sac_river_june_july<- mid_sac_river %>% filter(month(Datetime) %in% c(6,7))
unique(temp_dataset$SubRegion)
