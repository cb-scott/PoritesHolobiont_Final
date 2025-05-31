#create a map
library(rgdal)
library(tidyverse)
library(sp)
library(sf)

meta <- read.csv("data/Porites_Meta.csv") #again issues with 
meta <- data.frame(sapply(meta, function(x) trimws(x)))
meta <- meta %>% filter(SampleID != "") 
reef_locs <- meta %>% select(Site, Latitude, Longitude) %>% distinct()
reef_locs <- reef_locs[-5,]
reef_locs$Latitude <- as.numeric(reef_locs$Latitude)
reef_locs$Longitude <- as.numeric(reef_locs$Longitude)
reef_locs <- SpatialPoints(reef_locs[,c("Longitude", "Latitude")], proj4string = CRS("+init=epsg:4326"))


spatialpolydf <- read_sf("data/TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features/TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features.shp")
orpheus = subset(spatialpolydf, X_COORD>= 146.45 & X_COORD <=146.55)
orpheus = subset(orpheus, Y_COORD >= -18.68 & Y_COORD <= -18.5)
#reef_locs <- spTransform(reef_locs, CRS("+init=epsg:GGRS87"))
#orpheus <- orpheus %>% st_as_sf()
png("Visualizations/orpheusmap.png", height = 8, width = 8, units = "in", res = 500)
plot(orpheus, bg = 'lightskyblue3', col = 'white')
points(reef_locs$Longitude, reef_locs$Latitude, col = 'firebrick4', cex = 2, pch = 19)
dev.off()

png("figures/BasicOrphMap.png", height = 8, width = 8, units = "in", res= 500)
plot(orpheus['CHART_NAME'], bg= "white", col = 'lightgray', main = "")
dev.off()
#also plot an australia map:

map_data_es <- map_data('world')[map_data('world')$region == "Spain",]

library(maps)
library(ggplot2)

## make a df with only the country to overlap
map_data_es <- map_data('world')[map_data('world')$region == "Australia",]

## The map (maps + ggplot2 )
orph_island <- ggplot() +
  ## First layer: worldwide map
  geom_polygon(data = map_data("world"),
               aes(x=long, y=lat, group = group),
               color = '#9c9c9c', fill = '#f3f3f3') +
  ## Second layer: Country map
  geom_polygon(data = map_data_es,
               aes(x=long, y=lat, group = group),
               color = 'black', fill = 'gray') +
  coord_map() +
  coord_fixed(1.3,
              xlim = c(114, 153),
              ylim = c(-44, -12)) +
  theme_classic() + 
  geom_point(aes(y=-18.61, x=146.48), size = 2, col = 'red3') + 
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel,
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text= element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
  
orph_island
#theme(panel.background =element_rect(fill = 'blue'))

ggsave("figures/Australia_Map.png", orph_island, bg = "transparent", width = 1, height = 1, units = "in", dpi=500)
