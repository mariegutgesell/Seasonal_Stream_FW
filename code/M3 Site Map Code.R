##Map of Stream Sampling Sites - M3 Sites
##November 6, 2022

##Adapted from R code by K. Cazelles 
#### Packages ####
# read manipulate vector files tidyverse friendly
# introduce sf sfg and sfc objects
# rgeos
library(sf) # functions start with st_*
# read + manipulate raster files
library(raster)
# Visualize in a web browser
library(mapview)
# we already have a couple of nice tuto about it!
library(tidyverse)

library(leafem)
library(htmltools)
library(leaflet)
setwd("~/Desktop/Seasonal Stream Food Webs/Data/Land-use Data")
stream_gps <- read.csv("Stream Site Characteristics_2018.csv") %>%
  filter(sitecode == "AT" | sitecode == "HC" |sitecode=="P4")

stream_gps <- stream_gps%>%
  mutate(Site_Type = case_when(
    startsWith(sitecode, "HC") ~ "Low Impact",
    startsWith(sitecode, "P4") ~ "Mid Impact",
    startsWith(sitecode, "AT") ~ "High Impact"
  ) )

##stream maps -- all work, just look slightly different
streams_2018 <- st_as_sf(stream_gps, coords = c("Longitude", "Latitude"), crs = 4326)
mapview(streams_2018, zcol = "Site_Type") %>%
  addStaticLabels(label = stream_gps$Site.Code,
                  noHide = TRUE,
                  direction = 'right',
                  textOnly = TRUE,
                  textsize = "15px")

leaflet(stream_gps) %>%
  addTiles %>%
  addMarkers(~Longitude, ~Latitude, label = ~htmlEscape(Site_Type), labelOptions = labelOptions(noHide = TRUE))

leaflet(stream_gps) %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%
  addLabelOnlyMarkers(~Longitude, ~Latitude, label = ~as.character(Site_Type), 
                      labelOptions = labelOptions(noHide = T, direction = 'top', textOnly = T, textsize = "25px" )) %>%
  addScaleBar()

leaflet(stream_gps) %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%
  addCircleMarkers(~Longitude, ~Latitude, label = ~as.character(Site_Type), radius = 3, color = "black", stroke = FALSE, fillOpacity = 1.0) %>%
  addScaleBar() %>%
  addLabelOnlyMarkers(~Longitude, ~Latitude, label = ~htmlEscape(Site_Type), labelOptions = labelOptions(noHide =TRUE, textOnly = TRUE, textsize = "20px", direction = "right"))


