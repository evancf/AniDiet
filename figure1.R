# Working on a figure to show food web predictions...


load("m.mamm.pres.nat.RData")
rownames(m.mamm.pres.nat)
# Fix up rownames
rownames(m.mamm.pres.nat) <- gsub("/Users/efricke/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/Phylacine 1.2.1/Ranges/Present natural/",
                                  "", 
                                  rownames(m.mamm.pres.nat), fixed = T)
which(m.mamm.pres.nat["Thylacoleo carnifex",] == 1)


# extant <- web_current[[i]]
# extant_extinct <- web_pres_nat[[i]]
# 
# ?write.csv()
# setwd("~/Downloads")
# write.csv(extant, file = "extant.csv", row.names = F)
# write.csv(extant_extinct, file = "extant_extinct.csv", row.names = F)




# Diversion to get list of primarily carnivorous species 
# First, pull potential mammals from carnidiet
potential_mammals <- RCurl::getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/Supplementary%20data/Potential%20species%20list.csv") %>% 
  read.csv(text = .) %>% tibble()

mammal_predator <- c(gsub("_"," ",potential_mammals$Bin.),
                     impute_trait_data %>% 
                       filter(phylacine_binomial %in% 
                                te_dates$binomial & dphy_vertebrate > 50) %>% 
                       pull(phylacine_binomial))





# Some spatial data

library("raster")
library("sf")
continent_raster <- raster("~/Documents/git/AniDiet/Andermann/continent_shapes/raster_sum_all_regions.tif")

continent_proj <- continent_raster %>% 
  as("SpatialPixelsDataFrame") %>% 
  spTransform(CRS('+proj=longlat +datum=WGS84')) %>% 
  st_as_sf() %>% 
  mutate(cell = 1:nrow(.),
         terrestrial = raster_sum_all_regions > 0) %>% 
  dplyr::select(-raster_sum_all_regions)
cell_coords <- continent_proj %>% st_coordinates() %>% 
  as.data.frame() %>% 
  mutate(cell = continent_proj$cell)
colnames(cell_coords) <- c("lng", "lat", "cell")


cell_index_from_lat_lng <- function(coords){
  lat <- coords[1]
  lng <- coords[2]
  cell_coords$cell[which(abs(cell_coords$lng - lng) ==
                           min(abs(cell_coords$lng - lng)) &
                           abs(cell_coords$lat - lat) ==
                           min(abs(cell_coords$lat - lat)))]
}




# Use ggnetwork for network geometries so that nodes can be filtered out
# and they dont move location as with the D3 versions

library("ggnetwork")


make_ggnet <- function(x, pres_nat = T, panel_lab = "", ){
  
  
  intx_pres_nat <- web_pres_nat[[x]]
  intx_current <- web_current[[x]]
  intx_all <- bind_rows(intx_pres_nat, intx_current) %>% unique()
  
  # Get network geometry fitting all itneractions
  set.seed(4)
  ggnet_all <- ggnetwork(intx_all, 
                         arrow.gap = 0)
  
  # Get subset for pres_nat and 
  ggnet_pres_nat <- ggnet_all %>% filter(vertex.names %in% unlist(intx_pres_nat[,c("consumer_sp", "resource_sp")]))
  ggnet_pres_nat <- ggnet_pres_nat %>% filter(xend %in% ggnet_pres_nat$x)
  
  ggnet_current <- ggnet_all %>% filter(vertex.names %in% unlist(intx_current[,c("consumer_sp", "resource_sp")]))
  ggnet_current <- ggnet_current %>% filter(xend %in% ggnet_current$x)
  
  predator_color <- rgb(217,95,2, maxColorValue = 255)
  not_predator_color <- rgb(117,112,179, maxColorValue = 255)
  
  
  if(pres_nat) {
    ggplot() +
      geom_edges(data = ggnet_pres_nat,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color="grey50", curvature=0.1, size=0.2, alpha=1/2) +
      geom_nodes(data = ggnet_pres_nat,
                 aes(x=x, y=y),
                 size = 0.75,
                 color = ifelse(unique(ggnet_pres_nat$vertex.names) %in% mammal_predator, predator_color, not_predator_color)) +
      theme_void() +
      xlim(c(0,1)) +
      ylim(c(0,1.2)) + 
      theme(plot.background = element_rect(
        fill = NA,
        colour = "grey40",
        size = 0.75
      )) +
      coord_fixed() + 
      annotate(geom="text", x=0.5, y=1.15, label="Extant + Extinct", cex = 2)
  } else{
    ggplot() +
      geom_edges(data = ggnet_current,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color="grey50", curvature=0.1, size=0.2, alpha=1/2) +
      geom_nodes(data = ggnet_current,
                 aes(x=x, y=y),
                 size = 0.75,
                 color = ifelse(unique(ggnet_current$vertex.names) %in% mammal_predator, predator_color, not_predator_color)) +
      theme_void() +
      xlim(c(0,1)) +
      ylim(c(0,1.2)) + 
      theme(plot.background = element_rect(
        fill = NA,
        colour = "grey40",
        size = 0.75
      )) +
      coord_fixed() + 
      annotate(geom="text", x=0.5, y=1.15, label="Extant", cex = 2) + 
      annotate(geom="text", x=0.95, y=0.05, label=panel_lab, cex = 3, fontface = "bold")
  }

}




# Make a world map to show the locations on

library(rworldmap)
library("sf")
library(grid)
library(ggplotify)


worldmap <- rnaturalearth::ne_download(scale = 110,
                                       type = "countries",
                                       category = "cultural",
                                       destdir = tempdir(),
                                       load = TRUE,
                                       returnclass = "sf")

wm <- worldmap %>%
  select(CONTINENT) %>%
  filter(CONTINENT != "Antarctica") %>%
  st_combine() %>%
  st_transform(crs = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))


pdf("~/Downloads/map_examples.pdf", width = 7.2, height = 4.5)
par(mar = rep(8.5, 4))
plot(wm, col = "grey", border = "grey")

la_brea_coords <- c(34.063, -118.356)
la_brea_cell <- la_brea_coords %>% cell_index_from_lat_lng()

cerro_azul_coords <- c(1.2577, -72.612)
cerro_azul_cell <- cerro_azul_coords %>% cell_index_from_lat_lng()

lascaux_coords <- c(45.054, 1.168)
lascaux_coords <- c(49, 1)
lascaux_cell <- lascaux_coords %>% cell_index_from_lat_lng()

andriamamelo_coords <- c(-17.781896, 44.460088)
andriamamelo_cell <- andriamamelo_coords %>% cell_index_from_lat_lng()

#kimberley_coords <- c(-17.3, 127.3) # Not super exact coords
#kimberley_cell <- kimberley_coords %>% cell_index_from_lat_lng()

nsw_coords <- c(-34, 148.5)
nsw_cell <- nsw_coords %>% cell_index_from_lat_lng()

shanxi_coords <- c(38.5, 113)
shanxi_cell <- shanxi_coords %>% cell_index_from_lat_lng()





p1a <- make_ggnet(la_brea_cell, T)
p1b <- make_ggnet(la_brea_cell, F, panel_lab = "a")

p2a <- make_ggnet(cerro_azul_cell, T)
p2b <- make_ggnet(cerro_azul_cell, F, panel_lab = "b")

p3a <- make_ggnet(lascaux_cell, T)
p3b <- make_ggnet(lascaux_cell, F, panel_lab = "c")

p4a <- make_ggnet(andriamamelo_cell, T)
p4b <- make_ggnet(andriamamelo_cell, F, panel_lab = "d")

p5a <- make_ggnet(shanxi_cell, T)
p5b <- make_ggnet(shanxi_cell, F, panel_lab = "e")

p6a <- make_ggnet(nsw_cell, T)
p6b <- make_ggnet(nsw_cell, F, panel_lab = "f")




panel_dim <- 0.291

title_dim <- 0.4

x_diff <- 0.15
x_left <- 0.08
x_mid <- 0.43
x_right <- 0.92

x_title_left <- 0.005
x_title_mid <- 0.354
x_title_right <- 0.6935

y1 <- 0.2
y2 <- 0.3
y3 <- 0.7
y4 <- 0.8
y_title_diff <- 0.166

vp <- viewport(x_left, y3, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p1a))
upViewport()

vp <- viewport(x_title_left, y3 + y_title_diff, width = title_dim, height = title_dim/2)
pushViewport(vp)
grid.draw(as.grob(ggplot() + annotate(geom="text", x=0, y=0, label = "La Brea, California, USA", cex = 2.75, hjust = 0) + theme_void()))
upViewport()

vp <- viewport(x_left + x_diff, y3, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p1b))
upViewport()



vp <- viewport(x_left, y2, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p2a))
upViewport()

vp <- viewport(x_title_left, y2 - y_title_diff, width = title_dim, height = title_dim/2)
pushViewport(vp)
grid.draw(as.grob(ggplot() + annotate(geom="text", x=0, y=0, label = "Cerro Azul, Colombia", cex = 2.75, hjust = 0) + theme_void()))
upViewport()

vp <- viewport(x_left + x_diff, y2, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p2b))
upViewport()



vp <- viewport(x_mid, y4, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p3a))
upViewport()

vp <- viewport(x_title_mid, y4 + y_title_diff, width = title_dim, height = title_dim/2)
pushViewport(vp)
grid.draw(as.grob(ggplot() + annotate(geom="text", x=0, y=0, label = "Lascaux, France", cex = 2.75, hjust = 0) + theme_void()))
upViewport()

vp <- viewport(x_mid + x_diff, y4, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p3b))
upViewport()


vp <- viewport(x_mid, y1, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p4a))
upViewport()

vp <- viewport(x_title_mid, y1 - y_title_diff, width = title_dim, height = title_dim/2)
pushViewport(vp)
grid.draw(as.grob(ggplot() + annotate(geom="text", x=0, y=0, label = "Andriamamelo, Madagascar", cex = 2.75, hjust = 0) + theme_void()))
upViewport()

vp <- viewport(x_mid + x_diff, y1, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p4b))
upViewport()




vp <- viewport(x_right - x_diff, y3, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p5a))
upViewport()

vp <- viewport(x_title_right, y3 + y_title_diff, width = title_dim, height = title_dim/2)
pushViewport(vp)
grid.draw(as.grob(ggplot() + annotate(geom="text", x=0, y=0, label = "Shanxi, China", cex = 2.75, hjust = 0) + theme_void()))
upViewport()

vp <- viewport(x_right, y3, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p5b))
upViewport()



vp <- viewport(x_right - x_diff, y2, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p6a))
upViewport()

vp <- viewport(x_title_right, y2 - y_title_diff, width = title_dim, height = title_dim/2)
pushViewport(vp)
grid.draw(as.grob(ggplot() + annotate(geom="text", x=0, y=0, label = "New South Wales, Australia", cex = 2.75, hjust = 0) + theme_void()))
upViewport()

vp <- viewport(x_right, y2, width = panel_dim, height = panel_dim)
pushViewport(vp)
grid.draw(as.grob(p6b))
upViewport()






dev.off()







# In case you want to make a specific example food web (ie for talks)

make_ggnet2 <- function(x, pres_nat = T, panel_lab = ""){
  
  
  intx_pres_nat <- web_pres_nat[[x]]
  intx_current <- web_current[[x]]
  intx_all <- bind_rows(intx_pres_nat, intx_current) %>% unique()
  
  # Get network geometry fitting all itneractions
  set.seed(4)
  ggnet_all <- ggnetwork(intx_all, 
                         arrow.gap = 0)
  
  # Get subset for pres_nat and 
  ggnet_pres_nat <- ggnet_all %>% filter(vertex.names %in% unlist(intx_pres_nat[,c("consumer_sp", "resource_sp")]))
  ggnet_pres_nat <- ggnet_pres_nat %>% filter(xend %in% ggnet_pres_nat$x)
  
  ggnet_current <- ggnet_all %>% filter(vertex.names %in% unlist(intx_current[,c("consumer_sp", "resource_sp")]))
  ggnet_current <- ggnet_current %>% filter(xend %in% ggnet_current$x)
  
  predator_color <- rgb(217,95,2, maxColorValue = 255)
  not_predator_color <- rgb(117,112,179, maxColorValue = 255)
  
  
  if(pres_nat) {
    ggplot() +
      geom_edges(data = ggnet_pres_nat,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color="grey50", curvature=0.1, size=0.2, alpha=1/2) +
      geom_nodes(data = ggnet_pres_nat,
                 aes(x=x, y=y),
                 size = 0.75,
                 color = ifelse(unique(ggnet_pres_nat$vertex.names) %in% mammal_predator, predator_color, not_predator_color)) +
      theme_void() +
      xlim(c(0,1)) +
      ylim(c(0,1.2)) + 
      # theme(plot.background = element_rect(
      #   fill = NA,
      #   colour = "grey40",
      #   size = 0.75
      # )) +
      coord_fixed()
  } else{
    ggplot() +
      geom_edges(data = ggnet_current,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color="grey50", curvature=0.1, size=0.2, alpha=1/2) +
      geom_nodes(data = ggnet_current,
                 aes(x=x, y=y),
                 size = 0.75,
                 color = ifelse(unique(ggnet_current$vertex.names) %in% mammal_predator, predator_color, not_predator_color)) +
      theme_void() +
      xlim(c(0,1)) +
      ylim(c(0,1.2)) + 
      # theme(plot.background = element_rect(
      #   fill = NA,
      #   colour = "grey40",
      #   size = 0.75
      # )) +
      coord_fixed() + 
      annotate(geom="text", x=0.95, y=0.05, label=panel_lab, cex = 3, fontface = "bold")
  }
  
}
#pdf(file = "~/Downloads/brea.pdf", width = 1.2, height = 1.2)
png(file = "~/Downloads/brea.png", width = 1.2, height = 1.2, units = "in", res = 1000)
make_ggnet2(la_brea_cell, T)
dev.off()
png(file = "~/Downloads/brea2.png", width = 1.2, height = 1.2, units = "in", res = 1000)
make_ggnet2(la_brea_cell, F)
dev.off()

