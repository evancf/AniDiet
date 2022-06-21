# Packages
devtools::source_gist("71d758f65261a72ab7dc") # gist for ipak


# Georeferenced records -------------------------------------------------------
# 
source("./02c_pull_CarniDIET_data.R")


ipak(c("rnaturalearth","rnaturalearthdata",
       "sf"))
world <- ne_countries(scale = "medium", 
                      returnclass = "sf")

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

world <- world %>% filter(continent != "Antarctica")

carnidiet_summary <- carnidiet_tidy %>% 
  filter(!is.na(decimalLatitude),
         !is.na(decimalLongitude)) %>% 
  group_by(round(decimalLatitude), round(decimalLongitude)) %>% 
  tally() %>% 
  st_as_sf(coords = c("round(decimalLongitude)", "round(decimalLatitude)"),
           crs = projcrs)


p1 <- ggplot(data = world) +
  geom_sf(lwd = 0) +
  theme_void() +
  geom_sf(data = carnidiet_summary, 
          aes(color = n),
          #color = log(n)),
          alpha = 1,
          #shape = 15,
          size = 0.9,
          stroke = 0) + 
  scale_color_gradient(low = "blue", high = "red", trans = "log",
                       breaks=c(1,20,400),
                       name = "Georeferenced\nrecords")# + theme(legend.position = "left")



# Predator coverage -----------------------------------------------

# Load a bunch of data, as how it's loaded in 08 and 09
intx_short <- read.csv("./data/intx_short.csv")[,-1] %>% tibble()
impute_trait_data <- read.csv("./data/impute_trait_data.csv")[,-1] %>% tibble()
potential_mammals <- RCurl::getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/Supplementary%20data/Potential%20species%20list.csv") %>% 
  read.csv(text = .) %>% tibble()
load("data/m.mamm.current.no.intro.RData")
rownames(m.mamm.current.no.intro) <- gsub("/Users/efricke/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/Phylacine 1.2.1/Ranges/Current/", 
                                          "", 
                                          rownames(m.mamm.current.no.intro),
                                          fixed = T)
# Need to also get this for extinct species
# (NOTE that this code is replicated in '09_web_change_metrics.R')
te_dates <- read.table("Andermann/global_pyrate_species_list.txt", header = T) %>% tibble() %>% 
  bind_cols(read.table("Andermann/global_ts_te_dates.txt", header = T) %>% tibble())

# These data are accessible here: https://www.science.org/doi/10.1126/sciadv.abb2313
africa_cells <- raster::raster("Andermann/continent_shapes/raster_africa.tif")
africa_cells <- tibble(cell = which(raster::values(africa_cells) == 1),
                       continent = "africa")
australia_cells <- raster::raster("Andermann/continent_shapes/raster_australia.tif")
australia_cells <- tibble(cell = which(raster::values(australia_cells) == 1),
                          continent = "australia")
caribbean_cells <- raster::raster("Andermann/continent_shapes/raster_carribean.tif")
caribbean_cells <- tibble(cell = which(raster::values(caribbean_cells) == 1),
                          continent = "caribbean")
eurasia_cells <- raster::raster("Andermann/continent_shapes/raster_eurasia.tif")
eurasia_cells <- tibble(cell = which(raster::values(eurasia_cells) == 1),
                        continent = "eurasia")
madagascar_cells <- raster::raster("Andermann/continent_shapes/raster_madagascar.tif")
madagascar_cells <- tibble(cell = which(raster::values(madagascar_cells) == 1),
                           continent = "madagascar")
northamerica_cells <- raster::raster("Andermann/continent_shapes/raster_northamerica.tif")
northamerica_cells <- tibble(cell = which(raster::values(northamerica_cells) == 1),
                             continent = "northamerica")
oceanic_main_cells <- raster::raster("Andermann/continent_shapes/raster_oceanic_main.tif")
oceanic_main_cells <- tibble(cell = which(raster::values(oceanic_main_cells) == 1),
                             continent = "oceanic_main")
southamerica_cells <- raster::raster("Andermann/continent_shapes/raster_southamerica.tif")
southamerica_cells <- tibble(cell = which(raster::values(southamerica_cells) == 1),
                             continent = "southamerica")

# Put these together to make. This associates cell number with region 
continent_cells <- bind_rows(africa_cells, australia_cells, caribbean_cells,
                             eurasia_cells, madagascar_cells, northamerica_cells,
                             oceanic_main_cells, southamerica_cells)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
continent_cells <- continent_cells %>% 
  mutate(continent = continent %>% firstup())



# Get the mammal consumers with data
predators_with_data <- intx_short %>% 
  filter(consumed == 1) %>% 
  pull(consumer_sp) %>% 
  unique() %>% 
  sort()


# Get all mammal consumers
te_dates <- te_dates %>% 
  mutate(binomial = paste(id, taxon, sep = " ")) %>% 
  dplyr::select(-starts_with("ts"))

te_dates <- te_dates %>% filter(status == "extinct")

# Need to reconcile names
te_dates$binomial <- plyr::revalue(te_dates$binomial, 
                                   c("Alces scotti" = "Cervalces scotti",
                                     "Candiacervus SpII" = "Candiacervus spII",
                                     "Dicroceros sp" = "Dicroceros spA",
                                     "Geocapromys SP_A" = "Geocapromys spA",
                                     "Homo denisovans" = "Homo spDenisova",
                                     "Megaoryzomys Sp_Now" = "Megaoryzomys spA",
                                     "Nesophontes SP_A" = "Nesophontes spA",
                                     "Nesoryzomys Sp_B" = "Nesoryzomys spB",
                                     "Nesoryzomys Sp_C" = "Nesoryzomys spC",
                                     "Nesoryzomys Sp_D" = "Nesoryzomys spD",
                                     "Nothrotheriops shastense" = "Nothrotheriops shastensis",
                                     "Pachyarmaterium brasiliense" = "Pachyarmatherium brasiliense",
                                     "Peroryctes SP_NOW" = "Peroryctes spA",
                                     "Tapirus copei" = "Tapirus merriami",
                                     "Valgipes deformis" = "Valgipes bucklandi"
                                   ))

mammal_predator <- c(gsub("_"," ",potential_mammals$Bin.),
                     impute_trait_data %>% 
                       filter(phylacine_binomial %in% 
                                te_dates$binomial & dphy_vertebrate > 50) %>% 
                       pull(phylacine_binomial))

impute_trait_data$mammal_predator <- impute_trait_data$phylacine_binomial %in% mammal_predator


# Get range rasters for all predators
m_all_predators <- m.mamm.current.no.intro[rownames(m.mamm.current.no.intro) %in% mammal_predator, ]
m_all_predators <- m_all_predators[rowSums(m_all_predators) > 0,]


# Get range rasters for only predators with data
good_predators_with_data <- predators_with_data[predators_with_data %in% rownames(m_all_predators)]

m_predators_with_data <- m.mamm.current.no.intro[good_predators_with_data, ]
m_predators_with_data_genus <- m_all_predators[
  word(rownames(m_all_predators),1) %in% word(good_predators_with_data, 1), ]

# Also want summaries by genus

# Now do this at the genus level
m_all_predators_gen <- m_all_predators %>% as.data.frame() %>%
  mutate(genus = word(rownames(m_all_predators), 1))

# This takes a few minutes
m_all_predators_gen <- m_all_predators_gen %>%
  group_by(genus) %>%
  dplyr::summarize(across(1:51120, sum))

m_all_predators_gen <- m_all_predators_gen %>%
  as.data.frame()

rownames(m_all_predators_gen) <- as.character(m_all_predators_gen$genus)

m_all_predators_gen <- m_all_predators_gen %>%
  dplyr::select(-genus)

m_all_predators_gen <- apply(m_all_predators_gen, 2, function(x) ifelse(x > 0, 1, 0))

m_predators_with_data_gen <- m_all_predators_gen[unique(word(good_predators_with_data,1)), ]

# Make a dataframe that we will then join to coverage_rel below
cell <- 1:(142*360)
coverage_rel_genus <- tibble(cell = cell,
                             sum_all_predators_gen = colSums(m_all_predators_gen))


# Create a tibble that can be used for plotting
coverage_rel <- tibble(cell = cell,
                       sum_all_predators = colSums(m_all_predators),
                       sum_predators_with_data = colSums(m_predators_with_data),
                       sum_predators_with_data_genus = colSums(m_predators_with_data_genus)) %>% 
  left_join(coverage_rel_genus) %>% 
  mutate(x = (cell %% 360),
         y = 360 - floor(cell/360)) %>% 
  left_join(continent_cells) %>% 
  filter(!is.na(continent)) %>%
  filter(!is.na(sum_all_predators)) %>%
  filter(sum_all_predators > 0) %>%
  mutate(prop_predators_with_data = sum_predators_with_data / sum_all_predators,
         prop_predators_with_data_genus = sum_predators_with_data_genus / sum_all_predators) %>% 
  dplyr::na_if("NaN")


pixel_size <- 0.1

p2 <- coverage_rel %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = sum_all_predators), size = pixel_size, shape = 15) +
  theme_void() +
  coord_equal() +
  scale_color_distiller(palette = "Spectral", name = "Mammal\npredator\nspecies\nrichness",
                        limits = c(1, 23))# + theme(legend.position = "left")

p3 <- coverage_rel %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = prop_predators_with_data), size = pixel_size, shape = 15) +
  theme_void() +
  coord_equal() +
  scale_color_gradient(low = "grey", high = "blue",
                       breaks=c(0,0.5,1),
                       name = "Proportion\npredator\nspecies\ncovered")# + theme(legend.position = "left")

p4 <- coverage_rel %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = sum_all_predators_gen), size = pixel_size, shape = 15) +
  theme_void() +
  coord_equal() +
  scale_color_distiller(palette = "Spectral", name = "Mammal\npredator\ngenus\nrichness",
                        limits = c(1, 23))# + theme(legend.position = "left")

p5 <- coverage_rel %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = prop_predators_with_data_genus), size = pixel_size, shape = 15) +
  theme_void() +
  coord_equal() +
  scale_color_gradient(low = "grey", high = "blue",
                       breaks=c(0,0.5,1),
                       name = "Proportion\npredator\ngenera\ncovered")# + theme(legend.position = "left")


png(filename = "./figures/coverage.png", units = "in",
    width = 10, height = 8, res = 440)
ggpubr::ggarrange(
  p1, 
  labels = "A",
  ggpubr::ggarrange(p2, p3, p4, p5,
                    ncol = 2,
                    nrow = 2,
                    labels = c("B", "C", "D", "E"),
                    align = "hv"),
  heights = c(0.9, 1),
  nrow = 2)
dev.off()


# Model Prediction Uncertainty -------------------------------------------------

rm(list = ls()) # Get a clean slate for the following analyses
library("keras")
library("tidymodels")
library("tidyverse")


# Read in data 
intx_tidy <- read.csv("data/intx_tidy.csv")[,-1] %>% tibble()
intx_tidy$consumed <- factor(intx_tidy$consumed)


# Load trait data (Need to run 07_prep_data_for_interaction_modelling to get this)
impute_trait_data <- read.csv("data/impute_trait_data.csv")[,-1] %>% tibble()


# Diversion to get list of primarily carnivorous species 
# First, pull potential mammals from carnidiet
potential_mammals <- RCurl::getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/Supplementary%20data/Potential%20species%20list.csv") %>% 
  read.csv(text = .) %>% tibble()

# Need to also get this for extinct species
# (NOTE that this code is replicated in '09_web_boot_change_metrics.R')
te_dates <- read.table("Andermann/global_pyrate_species_list.txt", header = T) %>% tibble() %>% 
  bind_cols(read.table("Andermann/global_ts_te_dates.txt", header = T) %>% tibble())

te_dates <- te_dates %>% 
  mutate(binomial = paste(id, taxon, sep = " ")) %>% 
  dplyr::select(-starts_with("ts"))

te_dates <- te_dates %>% filter(status == "extinct")

# Need to reconcile names
te_dates$binomial <- plyr::revalue(te_dates$binomial, 
                                   c("Alces scotti" = "Cervalces scotti",
                                     "Candiacervus SpII" = "Candiacervus spII",
                                     "Dicroceros sp" = "Dicroceros spA",
                                     "Geocapromys SP_A" = "Geocapromys spA",
                                     "Homo denisovans" = "Homo spDenisova",
                                     "Megaoryzomys Sp_Now" = "Megaoryzomys spA",
                                     "Nesophontes SP_A" = "Nesophontes spA",
                                     "Nesoryzomys Sp_B" = "Nesoryzomys spB",
                                     "Nesoryzomys Sp_C" = "Nesoryzomys spC",
                                     "Nesoryzomys Sp_D" = "Nesoryzomys spD",
                                     "Nothrotheriops shastense" = "Nothrotheriops shastensis",
                                     "Pachyarmaterium brasiliense" = "Pachyarmatherium brasiliense",
                                     "Peroryctes SP_NOW" = "Peroryctes spA",
                                     "Tapirus copei" = "Tapirus merriami",
                                     "Valgipes deformis" = "Valgipes bucklandi"
                                   ))

mammal_predator <- c(gsub("_"," ",potential_mammals$Bin.),
                     impute_trait_data %>% 
                       filter(phylacine_binomial %in% 
                                te_dates$binomial & dphy_vertebrate > 50) %>% 
                       pull(phylacine_binomial))

impute_trait_data$mammal_predator <- impute_trait_data$phylacine_binomial %in% mammal_predator


# Load matrix of mammal presence - these are PHYLACINE data, slightly reformatted
load("data/m.mamm.pres.nat.RData")
load("data/m.mamm.current.RData")

# Fix up rownames
rownames(m.mamm.pres.nat) <- gsub("/Users/efricke/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/Phylacine 1.2.1/Ranges/Present natural/",
                                  "", 
                                  rownames(m.mamm.pres.nat), fixed = T)


# Remove those that forage in marine environments
impute_trait_data <- impute_trait_data %>% 
  filter(foraging_stratum != -1)

# Also can skip these marine foraging species in the range matrices
m.mamm.pres.nat <- m.mamm.pres.nat[rownames(m.mamm.pres.nat) %in% impute_trait_data$phylacine_binomial,]
m.mamm.current <- m.mamm.current[rownames(m.mamm.current) %in% impute_trait_data$phylacine_binomial,]

# # Lastly, remove humans 
m.mamm.pres.nat <- m.mamm.pres.nat[!grepl("Homo ", rownames(m.mamm.pres.nat)),]
m.mamm.current <- m.mamm.current[!grepl("Homo ", rownames(m.mamm.current)),]




# Also want a no vulnerable and endangered species version
endangered_categories <- c("CR", "EN", "EP", "EW", "EX", "VU", "CR (PE)")
# Need to pull back in COMBINE for this to reconcile names...
temp <- tempfile()
download.file("https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3344&file=ecy3344-sup-0001-DataS1.zip",
              temp)
trait_data <- read.table(unz(temp, "COMBINE_archives/trait_data_imputed.csv"), sep = ",", header = T) %>% tibble()
unlink(temp)
trait_data <- trait_data %>% 
  mutate(iucn2020_binomial = word(iucn2020_binomial, 1, 2))
iucn_categories <- read.csv("data/mamm.iucn.categories.csv", header = T)[,-1] %>% tibble()

# Get vector of endangered species names while accounting for different names across iucn and phylacine
endangered_spp <- c(left_join(trait_data, iucn_categories, by = c("iucn2020_binomial" = "iucn.bin")) %>%  
                      filter(iucn.category %in% endangered_categories) %>% 
                      dplyr::select("iucn2020_binomial", "phylacine_binomial") %>%
                      unlist(),
                    left_join(trait_data, iucn_categories, by = c("phylacine_binomial" = "iucn.bin")) %>%  
                      filter(iucn.category %in% endangered_categories) %>% 
                      dplyr::select("iucn2020_binomial", "phylacine_binomial") %>%
                      unlist()) %>% 
  unique()

# Get a similar matrix for this
m.mamm.no.endangered <- m.mamm.current[!rownames(m.mamm.current) %in% endangered_spp,]




# Function to get all pairwise consumer / resource combination given a species list
get_combn_df <- function(x, sp_combn = T, gen_combn = T){
  
  if(length(x) < 2){
    
    dat <- as.data.frame(matrix(nrow=0, ncol = 2))
    colnames(dat) <- c("consumer_sp", "resource_sp")
    
  } else{
    dat <- cbind(combn(x, 2), combn(x, 2)[2:1,]) %>% t() %>% as.data.frame()
    colnames(dat) <- c("consumer_sp", "resource_sp")
    
    if(sp_combn){
      dat$sp_combn <- paste(dat$consumer_sp, dat$resource_sp)
    }
    
    if(gen_combn){
      dat$gen_combn <- paste(word(dat$consumer_sp, 1), word(dat$resource_sp, 1))
    }
  }
  
  dat
}


# Want to get a grid of points that we'll actually get values for (because this take a while...)
xx <- seq(1, 360, by = 1)
yy <- seq(1, 142, by = 1)


# Get some cells to actually sample (not every single one at this point)
cells_to_sample <- apply(expand.grid(yy, xx),
                         1,
                         function(x) ((x[1] - 1) * 360 + x[2])) %>% sort()

# Also want to skip ones in the ocean
cells_to_sample <- cells_to_sample[cells_to_sample %in% which(colSums(m.mamm.pres.nat)>1)]



# Bootstrapping spatial uncertainty
library("tidymodels")
sp_col <- c("consumer_sp", "resource_sp", "mammal_predator_c", "mammal_predator_r")

n_resamp <- 10

# Will save the networks in a list object
web_boot_pres_nat <- list()
web_boot_current <- list()

set.seed(4)
for(tt in 1:n_resamp){
  
  # Resample dataset
  intx_resampled <- intx_tidy[sample(1:nrow(intx_tidy),
                                     size = nrow(intx_tidy),
                                     replace = T),]
  train_data <- intx_resampled[order(intx_resampled$consumed, decreasing = T),]
  train_data <- train_data[!duplicated(select(train_data, c(consumer_sp, resource_sp))),] %>%
    dplyr::select(-all_of(sp_col))
  
  dim(train_data)
  
  
  nnet_rec <-
    recipe(consumed ~ ., data = train_data) %>%
    step_YeoJohnson(all_numeric())%>%
    step_normalize(all_numeric()) %>%
    prep(training = train_data, retain = TRUE)
  
  
  train_normalized <- bake(nnet_rec, new_data = train_data, all_predictors())
  
  tensorflow::set_random_seed(4)
  nnet_fit <-
    mlp(epochs = 557, hidden_units = 10, dropout = 0.778, activation = "elu") %>%
    set_mode("classification") %>%
    set_engine("keras", verbose = T) %>%
    fit(consumed ~ ., data = bake(nnet_rec, new_data = NULL))
  
  
  # For loop to make webs at every study location
  
  for(i in cells_to_sample){
    ind <- which(cells_to_sample == i)
    
    # Get all mammal species in the pres_nat scenario
    spp_pres_nat <- rownames(m.mamm.pres.nat)[m.mamm.pres.nat[,i]]
    if(length(spp_pres_nat) >= 2){
      # Get combination df
      dat_pres_nat <- get_combn_df(spp_pres_nat, sp_combn = F, gen_combn = F)
      # Join with traits
      dat_pres_nat <- dat_pres_nat %>%
        left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>%
        left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>%
        filter(mammal_predator_c)
      if(dim(dat_pres_nat)[1] == 0){ # Cases where we dont have any predators
        tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
        colnames(tmp) <- c("consumer_sp", "resource_sp")
        web_boot_pres_nat[[i]] <- tmp
      } else{
        # Make predictions
        dat_pres_nat_normalized <- bake(nnet_rec, new_data = dat_pres_nat, all_predictors())
        dat_pres_nat <- dat_pres_nat %>%
          bind_cols(predict(nnet_fit, new_data = dat_pres_nat_normalized),
                    predict(nnet_fit, new_data = dat_pres_nat_normalized, type = "prob"))
        # Save this to output list
        web_boot_pres_nat[[i]] <- dat_pres_nat %>% filter(.pred_class == 1) %>%
          dplyr::select("consumer_sp", "resource_sp") %>% unique()
      }
      
    } else{ # Cases where there's only one species
      tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
      colnames(tmp) <- c("consumer_sp", "resource_sp")
      web_boot_pres_nat[[i]] <- tmp
    }
    
    
    # Get all mammal species in the current scenario
    spp_current <- rownames(m.mamm.current)[m.mamm.current[,i]]
    if(length(spp_current) >= 2){
      # Get combination df
      dat_current <- get_combn_df(spp_current, sp_combn = F, gen_combn = F)
      # Join with traits
      dat_current <- dat_current %>%
        left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>%
        left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>%
        filter(mammal_predator_c)
      if(dim(dat_current)[1] == 0){ # Cases where we dont have any predators
        tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
        colnames(tmp) <- c("consumer_sp", "resource_sp")
        web_boot_current[[i]] <- tmp
      } else{
        # Make predictions
        dat_current_normalized <- bake(nnet_rec, new_data = dat_current, all_predictors())
        dat_current <- dat_current %>%
          bind_cols(predict(nnet_fit, new_data = dat_current_normalized),
                    predict(nnet_fit, new_data = dat_current_normalized, type = "prob"))
        # Save this to output list
        web_boot_current[[i]] <- dat_current %>% filter(.pred_class == 1) %>%
          dplyr::select("consumer_sp", "resource_sp") %>% unique()
      }
      
    } else{ # Cases where there's only one species
      tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
      colnames(tmp) <- c("consumer_sp", "resource_sp")
      web_boot_current[[i]] <- tmp
    }
    
    print(ind / length(cells_to_sample))
    
  }
  
  #
  boot_filename <- paste0("data/boot/hindcast_webs_", tt, ".RData")
  save(web_boot_pres_nat,
       web_boot_current,
       file = boot_filename)
  
  print(tt)
  
}

boot_files <- list.files("./data/boot/", full.names = T)

load(boot_files[1])
cells_with_webs <- which(!lapply(web_boot_pres_nat, is.null) %>% unlist())

n_links_fun <- function(x) dim(x)[1]
n_nodes_fun <- function(x) length(unique(unlist(x)))
predator_prey_ratio_fun <- function(x) length(unique(x[,1])) / length(unique(x[,2]))



grid_current_links <- tibble(cell = cells_with_webs,
                             "1" = NA,
                             "2" = NA,
                             "3" = NA,
                             "4" = NA,
                             "5" = NA,
                             "6" = NA,
                             "7" = NA,
                             "8" = NA,
                             "9" = NA,
                             "10" = NA
)

grid_current_nodes <- grid_pres_nat_links <- grid_pres_nat_nodes <- grid_current_links

for(i in boot_files){
  # Get a numeric version of i
  ni <- which(i == boot_files)
  
  # Load data
  load(i)
  
  # Get cell-wise summaries 
  grid_current_links[,as.character(ni)] <- lapply(web_boot_current[cells_with_webs], n_links_fun) %>% unlist
  grid_current_nodes[,as.character(ni)] <- lapply(web_boot_current[cells_with_webs], n_nodes_fun) %>% unlist
  
  grid_pres_nat_links[,as.character(ni)] <- lapply(web_boot_pres_nat[cells_with_webs], n_links_fun) %>% unlist
  grid_pres_nat_nodes[,as.character(ni)] <- lapply(web_boot_pres_nat[cells_with_webs], n_nodes_fun) %>% unlist
  
  print(i)
}

cv_fun <- function(x){
  m <- mean(x[2:11], na.rm = T)
  sd <- sd(x[2:11], na.rm = T)
  
  sd / m
}



cv_current_links <- apply(grid_current_links, 1, cv_fun)
cv_current_nodes <- apply(grid_current_nodes, 1, cv_fun)

cv_pres_nat_links <- apply(grid_pres_nat_links, 1, cv_fun)
cv_pres_nat_nodes <- apply(grid_pres_nat_nodes, 1, cv_fun)

hist(cv_current_nodes)


cell <- 1:(142*360) %>% rev()
cv_rel <- left_join(tibble(cell = cell), 
                    tibble(cell = cells_with_webs,
                           cv_current_links = cv_current_links)) %>%
  left_join(tibble(cell = cells_with_webs,
                   cv_current_nodes = cv_current_nodes)) %>%
  left_join(tibble(cell = cells_with_webs,
                   cv_pres_nat_links = cv_pres_nat_links)) %>%
  left_join(tibble(cell = cells_with_webs,
                   cv_pres_nat_nodes = cv_pres_nat_nodes)) %>%
  mutate(x = (cell %% 360),
         y = 360 - floor(cell/360))


threshold <- 0.8
pixel_size <-  0.1

pu1 <- cv_rel %>% 
  mutate(cv_current_links = ifelse(cv_current_links > threshold, threshold, cv_current_links)) %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = cv_current_links), size = pixel_size, shape = 14) +
  theme_void() +
  coord_equal() +
  ggtitle("Number of links") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "grey", high = "blue",
                       breaks=c(0,threshold/2,threshold),
                       labels = c(0, threshold/2, paste0(">",0.8)),
                       na.value = "white",
                       name = "Coef. var.\ncurrent")# + theme(legend.position = "left")

pu2 <- cv_rel %>% 
  mutate(cv_current_nodes = ifelse(cv_current_nodes > threshold, threshold, cv_current_nodes)) %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = cv_current_nodes), size = pixel_size, shape = 14) +
  theme_void() +
  coord_equal() +
  ggtitle("Number of species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "grey", high = "blue",
                       breaks=c(0,threshold/2,threshold),
                       labels = c(0, threshold/2, paste0(">",0.8)),
                       na.value = "white",
                       name = "Coef. var.\ncurrent")# + theme(legend.position = "left")

pu3 <- cv_rel %>% 
  mutate(cv_pres_nat_links = ifelse(cv_pres_nat_links > threshold, threshold, cv_pres_nat_links)) %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = cv_pres_nat_links), size = pixel_size, shape = 14) +
  theme_void() +
  coord_equal() +
  scale_color_gradient(low = "grey", high = "blue",
                       breaks=c(0,threshold/2,threshold),
                       labels = c(0, threshold/2, paste0(">",0.8)),
                       na.value = "white",
                       name = "Coef. var.\npresent\nnatural")# + theme(legend.position = "left")

pu4 <- cv_rel %>% 
  mutate(cv_pres_nat_nodes = ifelse(cv_pres_nat_nodes > threshold, threshold, cv_pres_nat_nodes)) %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = cv_pres_nat_nodes), size = pixel_size, shape = 14) +
  theme_void() +
  coord_equal() +
  scale_color_gradient(low = "grey", high = "blue",
                       breaks=c(0,threshold/2,threshold),
                       labels = c(0, threshold/2, paste0(">",0.8)),
                       na.value = "white",
                       name = "Coef. var.\npresent\nnatural")# + theme(legend.position = "left")


png(filename = "./figures/uncertainty.png", units = "in",
    width = 10, height = 5, res = 440)
ggpubr::ggarrange(pu2, pu1, pu4, pu3, 
                  nrow = 2,
                  ncol = 2,
                  labels = c("A", "B", "C", "D"),
                  align = "hv")
dev.off()

