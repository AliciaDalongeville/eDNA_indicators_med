###############################################################################
# 
#  codes to reproduce analyses and figures of the Dalongeville et al. article accepted for publication in Journal of Applied Ecology: 
#  "Benchmarking eleven biodiversity indicators based on environmental DNA surveys: more diverse functional traits and evolutionary lineages inside marine reserves" 
#  
#   Code author : Alicia Dalongeville
#   Date : August 2022
###############################################################################


### Download libraries
library(ggplot2)
library(mapdata)
library(dplyr)
library(sf)
library(viridis)
library(colorBlindness)
library(spData)
library(cowplot)
library(ggsn)

## Load data
env <- read.csv("Data/environment.csv", header=T) %>%
  arrange(desc(protection))

env_reserve <- env %>%
  filter(protection == "reserve")

## Create a dataframe giving the position of each sites
sites_names <- c("Calvi", "Porquerolles", "Carry-le-Rouet", "Cap Roux", "Cerbere-Banyuls", "Riou ",
                 "Cerbicale", "Les Moines", "Lavezzi")
sites_coord <- env %>%
  group_by(Site, protection) %>%
  mutate(first_reserve = !duplicated(Site) & protection == "reserve")  %>%
  filter(first_reserve) %>%
  ungroup() %>%
  select(Site, longitude_start_DD,latitude_start_DD)

# Download the map for the Mediterranean Sea
wH <- map_data("worldHires",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea

### Map
x_title="Longitude"
y_title="Latitude"
hjust = c(1.2, -0.15,  1, -0.15, -0.1, 0.1, 0.8, 1.2, 1.3)
vjust = c(-0.9, 0, -1.1, 0, 0, 1.7, -2.3, 0, 0.3)

map_sampling <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim=c(2.3,9.3), ylim=c(41.2,44.2), ratio=1.2)+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+ #remove minor gridlines
  
  ## Add sampling sites
  geom_point(aes(x = longitude_start_DD, y = latitude_start_DD, 
                 color=protection, shape=protection), size=4,
             data=env_reserve)+
  geom_point(aes(x = longitude_start_DD, y = latitude_start_DD, 
                 color=protection, shape=protection), size = 2, 
             data=env)+
  geom_text(aes(x = longitude_start_DD, y = latitude_start_DD,label=sites_names),
            hjust=hjust, vjust=vjust, data=sites_coord, size=3.9)+
  xlab("Longitude") + ylab("Latitude") +
  scale_color_manual(values=c("cyan4","darkblue"), 
                     labels = c("fished area", "no-take marine reserve"),
                     name="Protection status") +
  scale_shape_manual(values=c(16,17), 
                     labels = c("fished area", "no-take marine reserve"),
                     name="Protection status") +
  theme(legend.position = "bottom") +
  north(env_reserve, location="bottomleft", scale=0.15, 
        anchor = c(x=4.3, y=41.2)) +
  scalebar(env_reserve, dist = 50, dist_unit = "km",
           transform = TRUE, model = "WGS84",
           st.bottom=T, location = "bottomright",
           anchor = c(x=6.3,y=41.3),
           st.dist = 0.05) +
  theme(legend.key=element_blank()) +
  theme(axis.text.x=element_text(colour="black",size=14))+
  theme(axis.text.y=element_text(colour="black", size=14))+
  theme(axis.title=element_text(colour="black",size=15)) +
  theme(legend.title=element_text(colour="black",size=16)) +
  theme(legend.text=element_text(colour="black",size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

###################################
## Inset map
# border of the area of interest
inset_map <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = "gray30") +
  coord_fixed(xlim=c(-8,37), ylim=c(29,47), ratio=1.2)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_blank(), axis.text=element_blank(),
        axis.ticks = element_blank())+
  geom_rect(aes(xmin = 2.3, xmax = 9.3, ymin = 41.2, ymax = 44.2), 
            fill = NA, colour = "black", size = 1.2) 

#################
## Combine sampling map and insert
fig1 = ggdraw() +
  draw_plot(map_sampling) +
  draw_plot(inset_map, x = 0.061, y = 0.29, width = 0.3, height = 0.2)

#fig1

### Save the map
ggsave(filename = "04_map/Figure1_Sampling.png", 
       plot = fig1,
       width = 10, 
       height = 10,
       dpi = 300)

