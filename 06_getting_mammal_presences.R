# Want to get mammal species lists at each site.
# Will do this outside of this code, as it requires IUCN mammal range
# shapefiles that have access controlled by IUCN.


# This is the code that I'll execute on a cluster
library("raster")
library("rgdal")
library("sf")


# Load mammal range data
# Change to where your IUCN 
mamm.ranges <- shapefile("/nfs/efricke-data/species ranges/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp")


# Bring in study coordinates
library('tidyverse')
library("RCurl")

# Pull in the CarniDIET data
carnidiet_coords <- getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/CarniDIET%201.0.csv") %>% 
  read.csv(text = .) %>% tibble() %>% 
  dplyr::select("sourceTitle",
                "sourcePrimaryReference",
                "verbatimLocality", "protectedAreaHigher",
                "decimalLatitude", "decimalLongitude") %>% unique() %>% 
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>% 
  mutate(carni.id = paste(decimalLatitude, decimalLongitude))

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

carnidiet_coords <- st_as_sf(carnidiet_coords,
                             coords = c("decimalLongitude", "decimalLatitude"),
                             crs = projcrs)

# Want a 5km buffer (using appropriate projection for the study site)
utm_prj4 <- function(x) {
  
  coords <- cbind(x, st_coordinates(x))
  
  long <- coords$X
  lat <- coords$Y
  
  zone <- if(lat >= 56 && lat < 64 && long >= 3 && long < 12){x <- 32} else if(
    lat >= 72 && lat < 84 && long >= 0 && long < 9) {x <- 31} else if(
      lat >= 72 && lat < 84 && long >= 9 && long < 21) {x <- 33} else if(
        lat >= 72 && lat < 84 && long >= 21 && long < 33) {x <- 35} else if(
          lat >= 72 && lat < 84 && long >= 33 && long < 42) {x <- 37} else{
            x <- (floor((long + 180)/6) %% 60) + 1
          }
  prj <- purrr::map2_chr(zone, lat, function(y, z){
    if (z >= 0){
      paste0("+proj=utm +zone=", y, " +datum=WGS84 +units=m +no_defs")
    } else{
      paste0("+proj=utm +zone=", y, " +south", " +datum=WGS84 +units=m +no_defs")
    }})
  prj
}

carnidiet_coords <- map2(1:nrow(carnidiet_coords), utm_prj4(carnidiet_coords), function(x, y){
  st_transform(carnidiet_coords[x,], y)
})
carnidiet_coords <- map(carnidiet_coords, ~ st_buffer(., 5000))
carnidiet_coords <- map(carnidiet_coords, ~ st_transform(., projcrs)) # Go back to orig projection
carnidiet_coords <- do.call(rbind, carnidiet_coords) # Go from list back to df


# Spatial join and spread to make a big matrix of species presence
mamm.ranges.sf <- st_as_sf(mamm.ranges)

mamm.presence.by.carnicoords <- carnidiet_coords %>% st_join(filter(mamm.ranges.sf,
                                                                    presence == 1))

mamm.presence.by.carnicoords.wide <- mamm.presence.by.carnicoords[,c("carni.id", "binomial")]
st_geometry(mamm.presence.by.carnicoords.wide) <- NULL
mamm.presence.by.carnicoords.wide <- mamm.presence.by.carnicoords.wide[!duplicated(mamm.presence.by.carnicoords.wide),]
mamm.presence.by.carnicoords.wide$value <- T
mamm.presence.by.carnicoords.wide <- spread(data = mamm.presence.by.carnicoords.wide, key = carni.id, value = value, fill = F)

setwd("/nfs/efricke-data/species ranges/presence.by.carnicoords")
save(mamm.presence.by.carnicoords.wide, file = "mamm.presence.by.carnicoords.wide.RData")
