library("tidyverse")

# Load raster data for region-level analyses ------------------------------

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

# Load hindcast webs ------------------------------

load("data/hindcast_webs.RData")


# Load trait data ------------------------------

impute_trait_data <- read.csv("data/impute_trait_data.csv")[,-1] %>% tibble()


# Load matrix of mammal presence ------------------------------

load("data/m.mamm.pres.nat.RData")
rownames(m.mamm.pres.nat) <- rownames(m.mamm.pres.nat) %>% basename(.)


# Load extinction dates from Andermann ------------------------------

te_dates <- read.table("Andermann/global_pyrate_species_list.txt", header = T) %>% tibble() %>% 
  bind_cols(read.table("Andermann/global_ts_te_dates.txt", header = T) %>% tibble())

te_dates <- te_dates %>% 
  mutate(binomial = paste(id, taxon, sep = " ")) %>% 
  dplyr::select(-starts_with("ts"))

te_dates <- te_dates %>% filter(status == "extinct")

# Need to reconcile names
te_dates$binomial[!te_dates$binomial %in% impute_trait_data$phylacine_binomial]

filter(impute_trait_data, word(phylacine_binomial,1) == "Tapirus") %>% select(phylacine_binomial, iucn2020_binomial)
filter(impute_trait_data, word(phylacine_binomial,1) == "copei") %>% select(phylacine_binomial, iucn2020_binomial)
filter(impute_trait_data, word(phylacine_binomial,1) == "Cervalces") %>% select(phylacine_binomial, iucn2020_binomial)

filter(te_dates, word(binomial,1) == "Tapirus") %>% select(binomial)

# A few naming things to reconcile
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


# Decide on some relevant focal years (rather than every year)

focal_years <- c(0, 400, 3000, 10000, 23000, 46000, 80000, 126000)
focal_years <- (seq(0, (126000)^(1/3), length.out = 40))^3 %>% round()


# Get metrics for the cells for which we were actually able to get networks.
# In some it wasn't possible because there aren't terrestrial mammalian carnivore
# interactions (e.g., oceanic islands)

cells_with_webs <- which(!lapply(web_pres_nat, is.null) %>% unlist())

n_links_fun <- function(x) dim(x)[1]
n_nodes_fun <- function(x) length(unique(unlist(x)))
predator_prey_ratio_fun <- function(x) length(unique(x[,1])) / length(unique(x[,2]))


grid_web_metrics <- tibble(cell = cells_with_webs,
                           n_links_pres_nat = lapply(web_pres_nat[cells_with_webs], n_links_fun) %>% unlist,
                           n_nodes_pres_nat = lapply(web_pres_nat[cells_with_webs], n_nodes_fun) %>% unlist,
                           predator_prey_ratio_pres_nat = lapply(web_pres_nat[cells_with_webs], predator_prey_ratio_fun) %>% unlist,
                           #mean_chain_length_pres_nat = NA,
                           n_links_current = lapply(web_current[cells_with_webs], n_links_fun) %>% unlist,
                           n_nodes_current = lapply(web_current[cells_with_webs], n_nodes_fun) %>% unlist,
                           predator_prey_ratio_current = lapply(web_current[cells_with_webs], predator_prey_ratio_fun) %>% unlist,
                           #mean_chain_length_current = NA
                           n_links_no_endangered = lapply(web_no_endangered[cells_with_webs], n_links_fun) %>% unlist,
                           n_nodes_no_endangered = lapply(web_no_endangered[cells_with_webs], n_nodes_fun) %>% unlist,
                           predator_prey_ratio_no_endangered = lapply(web_no_endangered[cells_with_webs], predator_prey_ratio_fun) %>% unlist,
                           #mean_chain_length_no_endangered = NA
) %>% as.data.frame()



# Decide how many of 100 extinction years to sample
te_samp <- 1 #paste("te", 1, sep = ".")

# Get matrix of link and node metrics

past_metrics <- expand_grid(cells_with_webs, focal_years, te_samp) %>% 
  rename(cell = cells_with_webs) %>% 
  mutate(n_links = NA, n_nodes = NA,
         n_links_null = NA, n_nodes_null = NA,
         predator_prey_ratio = NA, predator_prey_ratio_null = NA)

# set.seed(4)
# for(i in cells_with_webs){
#   # Get a vector of all species that naturally occur at this cell
#   all_sp <- rownames(m.mamm.pres.nat)[m.mamm.pres.nat[,i]]
# 
#   for(k in te_samp){ # Switching the indexing a little
# 
#     te_cell_samp <- paste("te", sample(1:99, 1), sep = ".")
# 
#     # Create vectors to be populated with a record of the mammal species
#     # that went extinct in a given cell, or simulated to go extinct
#     # under a random simulation
#     true_extinct <- c()
#     null_extinct <- c()
# 
#     for(j in rev(focal_years)){
# 
#       ind <- which(past_metrics$cell == i & past_metrics$focal_years == j & past_metrics$te_samp == k)
# 
#       # Globally extinct species by this year
#       global_extinct_j <- te_dates$binomial[te_dates[,te_cell_samp] >= j]
# 
#       # Locally extant and extinct species by this year
#       local_extant_j <- all_sp[!all_sp %in% global_extinct_j]
#       local_extinct_j <- all_sp[all_sp %in% global_extinct_j]
# 
#       # Null extant species at this year (random extinctions)
#       # Sample additional species from among those extant in the null scenario
#       new_null_extinct <- all_sp[!all_sp %in% null_extinct] %>%
#         sample(length(local_extinct_j) - length(null_extinct), replace = F)
# 
#       # Add these to the null_extinct vector, which will be updated at the next year
#       null_extinct <- c(null_extinct, new_null_extinct)
#       #local_extant_j_null <- all_sp %>% sample(length(local_extant_j), replace = F)
# 
#       # Determine if the consumer and resource are absent at the given date
#       c_bool <-  web_pres_nat[[i]][,1] %in% local_extant_j
#       r_bool <-  web_pres_nat[[i]][,2] %in% local_extant_j
# 
#       past_web <- web_pres_nat[[i]][c_bool & r_bool,]
# 
#       past_metrics$n_links[ind] <- past_web %>% n_links_fun()
#       past_metrics$n_nodes[ind] <- past_web %>% n_nodes_fun()
#       past_metrics$predator_prey_ratio[ind] <- past_web %>% predator_prey_ratio_fun()
# 
# 
# 
# 
# 
#       # Perform similar calculations for web with randomized extinctions
#       c_bool_null <-  !web_pres_nat[[i]][,1] %in% null_extinct
#       r_bool_null <-  !web_pres_nat[[i]][,2] %in% null_extinct
# 
#       past_web_null <- web_pres_nat[[i]][c_bool_null & r_bool_null,]
# 
#       past_metrics$n_links_null[ind] <- past_web_null %>% n_links_fun()
#       past_metrics$n_nodes_null[ind] <- past_web_null %>% n_nodes_fun()
#       past_metrics$predator_prey_ratio_null[ind] <- past_web_null %>% predator_prey_ratio_fun()
# 
# 
#     }
#   }
#   print(round(i / tail(cells_with_webs, 1), digits = 3))
# }
# 
# 
# write.csv(file = "data/past_metrics.csv", past_metrics)

past_metrics <- read.csv(file = "data/past_metrics.csv")[,-1]

# Left join to pres_nat values

prop_past_metrics <- past_metrics %>% 
  left_join(grid_web_metrics) %>% 
  mutate(n_links_change = n_links / n_links_pres_nat,
         n_nodes_change = n_nodes / n_nodes_pres_nat,
         n_links_change_null = n_links_null / n_links_pres_nat,
         n_nodes_change_null = n_nodes_null / n_nodes_pres_nat) %>% 
  dplyr::na_if("NaN") %>% 
  left_join(continent_cells)

prop_past_metrics <- prop_past_metrics %>%
  group_by(cell, focal_years) %>% # only relevant if more than one cell/focal year combination is sampled - in other words, te_samp > 1
  summarise(continent = first(continent),
            n_links_change = mean(n_links_change, na.rm = T),
            n_nodes_change = mean(n_nodes_change, na.rm = T),
            n_links_change_null = mean(n_links_change_null, na.rm = T),
            n_nodes_change_null = mean(n_nodes_change_null, na.rm = T),
            n_nodes_change_current = mean(n_nodes_current/n_nodes_pres_nat, na.rm = T),
            n_links_change_current = mean(n_links_current/n_links_pres_nat, na.rm = T)) %>% 
  dplyr::na_if("NaN")

prop_past_metrics$continent <- factor(prop_past_metrics$continent, 
                                      levels = c("Africa",
                                                 "Australia",
                                                 "Caribbean",
                                                 "Eurasia",
                                                 "Madagascar",
                                                 "Northamerica",
                                                 "Oceanic_main",
                                                 "Southamerica")#, 
                                      # labels = c("Africa",
                                      #            "Australia",
                                      #            #"Caribbean",
                                      #            "Eurasia",
                                      #            "Madagascar",
                                      #            "North America",
                                      #            "Oceanic areas",
                                      #            "South America")
                                      )

# Human arrival dates from here: https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01566
arrival_kya <-  tibble(continent = levels(prop_past_metrics$continent),
                       min = c(NA,
                               60,
                               7, # Carribbean
                               64,
                               8,
                               20,
                               NA, # Oceanic
                               16),
                       max = c(NA,
                               56,
                               4, # Caribbean
                               60,
                               4,
                               16,
                               NA, # Oceanic
                               12
                       ))




# Decline over time figure ---------------------------

link_col <- rgb(217,95,2, maxColorValue = 255)
link_col_poly <- rgb(217,95,2, maxColorValue = 255, alpha = 200)
node_col <- rgb(117,112,179, maxColorValue = 255)
node_col_poly <- rgb(117,112,179, maxColorValue = 255, alpha = 200)
null_col <- rgb(27,158,119, maxColorValue = 255)
null_col_poly <- rgb(27,158,119, maxColorValue = 255, alpha = 200)

pdf(file = "figures/change over time.pdf", width = 3.5, height = 6)

xrep <- 4
layout(matrix(data = c(rep(19, xrep*2+2),
                       c(rep(c(20, rep(1, xrep), rep(2, xrep),3),3)),
                       c(rep(c(20, rep(4, xrep), rep(5, xrep),6),3)),
                       c(rep(c(20, rep(7, xrep), rep(8, xrep),9),3)),
                       c(rep(c(20, rep(10, xrep), rep(11, xrep),12),3)),
                       c(rep(c(20, rep(13, xrep), rep(14, xrep),15),3)),
                       c(rep(c(20, rep(16, xrep), rep(17, xrep),18),3)),
                       rep(19, xrep*2+2)),
              ncol = xrep*2+2,
              byrow = T))



#par(mfrow=c(6,2))


counter <- 1
bottom_row_inds <- 1#5:6
left_col_inds <- c(1,3,5)
focal_year_labels <- c(0, 400, 10000, 46000, 126000)

samp_mean <- function(x,i){mean(x[i])}

se_fun <- function(x, dir = "hi"){
  if(dir == "hi"){
    return(mean(x, na.rm = T) + sd(x, na.rm = T)/sqrt(length(x[!is.na(x)])) * 1.96)
  } 
  if(dir == "lo"){
    return(mean(x, na.rm = T) - sd(x, na.rm = T)/sqrt(length(x[!is.na(x)])) * 1.96)
  } 
  
} 
continent_levels <- c("Africa", "Eurasia", "Australia", "Northamerica", "Southamerica", "Madagascar")
for(i in continent_levels){
  
  dat <- prop_past_metrics %>% filter(continent == i) %>% 
    mutate(y1 = n_links_change, #log(n_links_change),
           y2 = n_nodes_change, #log(n_nodes_change),
           y3 = n_links_change_null,
           y4 = n_nodes_change_null,
           x = c(focal_years ^ (1/3)))
  
  if(i %in% c("Oceanic_main", "Caribbean")) next()
  
  par(mar = c(0,2,1,0))
  
  ylo <- 0.1
  
  plot(NA, 
       #data = dat,
       xlim = (c(126000, 0) ^(1/3)),
       ylim = c(ylo, 1), #log(c(0.05, 1)),
       pch = 16,
       col = rgb(0,0,1,0.1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       frame = F,
       las = 1)
  
  rect(xleft = (filter(arrival_kya, continent == i)$min * 1000)^(1/3),
       xright = (filter(arrival_kya, continent == i)$max * 1000)^(1/3),
       ybottom = ylo + 0.1, #ifelse(i == "Madagascar", ylo + 0.2, ylo + 0.1),
       ytop = 1.01,
       border = F,
       col = "lightgrey")
  
  if(i == "Madagascar"){
    mtext("Species", side = 1,
          col = node_col,
          line = 0.5, font = 1,
          cex = 0.90, xpd = T)
  }
  
  if(counter == 2){
    
    
    text(x = 38000^(1/3),
         y = 0.35,
         pos = 4,
         font = 3,
         labels = "H. sapiens\narrival")
    arrows(x0 = 35000^(1/3),
           x1 = 58000^(1/3),
           y0 = 0.35 + 0.03,
           y1 = 0.47,
           length = 0.06)
  }
  
  
  
  if(counter == 1) mtext("Thousand years before present", side = 3, line = 2.15, cex = 0.9, adj = -0.2, xpd = T)
  if(counter == 3) mtext("Percent change due to extinction", side = 2, line = 3.5, cex = 0.9, adj = 0.8, xpd = T)
  
  if(counter %in% bottom_row_inds){
    axis(3, 
         at = rev(focal_year_labels^(1/3)), 
         labels = rev(focal_year_labels/1000),
         mgp = c(3,0.75,0))
  } else{
    axis(3, 
         at = rev(focal_year_labels^(1/3)), 
         labels = NA)
  }
  
  #if(counter %in% left_col_inds){
  axis(2, at = c(1, 0.8, 0.6, 0.4, 0.2), #log(c(1, 0.5, 0.2, 0.05)), 
       labels = c("0%", "", "-40%", "", "-80%"), las = 1)
  # } else{
  #   axis(2, at = c(1, 0.8, 0.6, 0.4, 0.2), #log(c(1, 0.5, 0.2, 0.05)), 
  #        labels = NA, las = 1)
  # }
  
  xx <- c(focal_years ^ (1/3))
  
  
  yy1_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y1, na.rm = T)) %>% ungroup() %>% select(mean) %>% unlist()
  yy1_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y1, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy1_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y1, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  yy2_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y2, na.rm = T)) %>% ungroup() %>% select(mean) %>% unlist()
  yy2_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y2, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy2_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y2, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  yy3_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y3, na.rm = T)) %>% ungroup() %>% select(mean) %>% unlist()
  yy3_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y3, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy3_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y3, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  yy4_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y4, na.rm = T)) %>% ungroup() %>% select(mean) %>% unlist()
  yy4_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y4, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy4_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y4, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  
  
  polygon(x = c(xx, rev(xx)), y = c(yy4_sehi, rev(yy4_selo)),
          border = F,
          col = "lightgrey")
  lines(x = xx, y = yy4_mean, col = "darkgrey")
  
  polygon(x = c(xx, rev(xx)), y = c(yy2_sehi, rev(yy2_selo)),
          border = F,
          col = node_col_poly)
  lines(x = xx, y = yy2_mean, col = node_col)
  
  
  plot(NA, 
       #data = dat,
       xlim = (c(126000, 0) ^(1/3)),
       ylim = c(ylo, 1), #log(c(0.05, 1)),
       pch = 16,
       col = rgb(0,0,1,0.1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       frame = F,
       xpd = T,
       las = 1)
  
  
  axis(2, at = c(1, 0.8, 0.6, 0.4, 0.2), #log(c(1, 0.5, 0.2, 0.05)), 
       labels = NA, las = 1)
  
  rect(xleft = (filter(arrival_kya, continent == i)$min * 1000)^(1/3),
       xright = (filter(arrival_kya, continent == i)$max * 1000)^(1/3),
       ybottom = ylo + 0.1,
       ytop = 1.01,
       border = F,
       col = "lightgrey")
  
  if(i == "Madagascar"){
    mtext("Links", side = 1, 
          line = 0.5, font = 1,
          col = link_col,
          cex = 0.90, xpd = T)
  }
  
  polygon(x = c(xx, rev(xx)), y = c(yy3_sehi, rev(yy3_selo)),
          border = F,
          col = "lightgrey")#null_col_poly)
  lines(x = xx, y = yy3_mean, col = "darkgrey")#null_col)
  
  polygon(x = c(xx, rev(xx)), y = c(yy1_sehi, rev(yy1_selo)),
          border = F,
          col = link_col_poly)
  lines(x = xx, y = yy1_mean, col = link_col)
  
  
  
  
  
  if(counter %in% bottom_row_inds){
    axis(3, 
         at = rev(focal_year_labels^(1/3)), 
         labels = rev(focal_year_labels/1000),
         mgp = c(3,0.75,0))
  } else{
    axis(3, 
         at = rev(focal_year_labels^(1/3)), 
         labels = NA)
  }
  
  
  
  
  
  
  if(counter == 2){
    
    xxo <- 33
    yyo <- 0.68
    text(x = xxo,
         y = yyo,
         pos = 4,
         font = 3,
         labels = "Observed")
    arrows(x0 = xxo - 27,
           x1 = xxo - 31,
           y0 = yyo + 0.02,
           y1 = yyo + 0.09,
           length = 0.06)
    
    xxn <- 20
    yyn <- 0.95
    text(x = xxn,
         y = yyn,
         pos = 4,
         font = 3,
         labels = "Null")
    arrows(x0 = xxn - 13,
           x1 = xxn - 18,
           y0 = yyn - 0.01,
           y1 = yyn - 0.08,
           length = 0.06)
  }
  
  counter <- counter + 1
  
  par(mar = c(0,0,0,0))
  plot(NA,
       xlim = c(0,1),
       ylim = c(0,1),
       frame = F,
       xaxt = "n",
       yaxt = "n")
  region_labels <- c("Africa", "Australia", NA, "Eurasia", "Madagascar",
                     "N. America", NA, "S. America")
  text(0.6, 0.4, region_labels[which(levels(prop_past_metrics$continent) ==i)], 
       srt = 90,
       cex = 1.15, xpd = T, pos = 3)
  
}

dev.off()



# Want to calculate how much less link loss would occur had species declined randomly
# "Comparing randomized extinction to the global average effect of all observed 
# extinctions, randomized extinction resulted in X% fewer links lost and X% 
# fewer species lost from food webs"
dd <- 1 - prop_past_metrics %>% filter(focal_years == 0) %>% 
  select(n_links_change:n_nodes_change_null) %>% 
  colMeans(na.rm=T)

dd[2] / dd[4]
dd[3] / dd[5]

# Could also frame this in the inverse sense (how many less there would be if
# it were random)
1 - dd[4] / dd[2] # links
1 - dd[5] / dd[3]




# Percent decline attribution ---------------------------

# For the following section, the scenarios are coded this way

# pn: present natural - the counterfactual scenario of species composition as
# it would occur today had defaunation / introductions not occurred

# ex: extinction only scenario. Had extinctions occured, but extant species
# still filled their natural range

# c: current. Species composition as it occurs today, with ranges affected
# by defaunation and species introductions

# e: endangered species extinction scenario. This represents the current
# scenario, but with all vulnerable and endangered species extinct

cell <- 1:(142*360) %>% rev()
grid_rel <- left_join(tibble(cell = cell), grid_web_metrics) %>% 
  mutate(x = (cell %% 360),
         y = 360 - floor(cell/360)) %>% 
  left_join(continent_cells)

grid_rel_extinction_only <- past_metrics %>% 
  tibble() %>% 
  filter(focal_years == 0) %>% 
  rename(n_nodes_extinct_only = n_nodes,
         n_links_extinct_only = n_links) %>% 
  select(cell, n_nodes_extinct_only, n_links_extinct_only)


grid_rel <- grid_rel %>% 
  left_join(grid_rel_extinction_only) %>% 
  # Get absolute declines percent declines relative to the present natural
  mutate(abs_expn_nodes = (n_nodes_extinct_only - n_nodes_pres_nat) / n_nodes_pres_nat,
         abs_expn_links = (n_links_extinct_only - n_links_pres_nat) / n_links_pres_nat,
         abs_expn_ld = ((n_links_extinct_only/n_nodes_extinct_only) - (n_links_pres_nat/n_nodes_pres_nat)) / (n_links_pres_nat/n_nodes_pres_nat),
         
         abs_cex_nodes = (n_nodes_current - n_nodes_extinct_only) / n_nodes_extinct_only,
         abs_cex_links = (n_links_current - n_links_extinct_only) / n_links_extinct_only,
         abs_cex_ld = ((n_links_current/n_nodes_current) - (n_links_extinct_only/n_nodes_extinct_only)) / (n_links_extinct_only/n_nodes_extinct_only),
         
         
         abs_ec_nodes = (n_nodes_no_endangered - n_nodes_current) / n_nodes_current,
         abs_ec_links = (n_links_no_endangered - n_links_current) / n_links_current,
         abs_ec_ld = ((n_links_no_endangered/n_nodes_no_endangered) - (n_links_current/n_nodes_current)) / (n_links_current/n_nodes_current),
  )


# Reshape to get three panel
grid_rel_longer <- grid_rel %>% 
  select(cell, continent,
         abs_expn_nodes, abs_expn_links, abs_expn_ld,
         abs_cex_nodes, abs_cex_links, abs_cex_ld,
         abs_ec_nodes, abs_ec_links, abs_ec_ld) %>% 
  pivot_longer(c(abs_expn_nodes, abs_expn_links, abs_expn_ld,
                 abs_cex_nodes, abs_cex_links, abs_cex_ld,
                 abs_ec_nodes, abs_ec_links, abs_ec_ld), 
               names_to = "metric", 
               values_to = "value") %>% 
  dplyr::na_if("NaN")  %>% 
  mutate(x = (cell %% 360),
         y = 360 - floor(cell/360))

grid_rel_longer$metric <- factor(grid_rel_longer$metric, 
                                 levels = c("abs_expn_nodes", "abs_expn_links", "abs_expn_ld",
                                            "abs_cex_nodes", "abs_cex_links", "abs_cex_ld",
                                            "abs_ec_nodes", "abs_ec_links", "abs_ec_ld"), 
                                 labels = c("Extinction: Species", "Extinction: Links", "Extinction: Linkage Density",
                                            "Range Change: Species", "Range Change: Links", "Range Change: Linkage Density",
                                            "Endangerment: Species", "Endangerment: Links", "Endangerment: Linkage Density"))


pdf(file = "figures/relative change maps.pdf", width = 4.25, height = 2.75)
grid_rel_longer %>% 
  mutate(value = ifelse(value > 1, 1, value)) %>%
  mutate(value = ifelse(value < -0.75, -0.75, value)) %>%
  mutate(value = value + 1) %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(colour = log(value)), size = 0.005, shape = 15) +
  facet_wrap(vars(metric), nrow=3) +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0,
                         na.value = "white",
                         name = "Percent\nchange\n                  \n\n",
                         breaks=c(log(2), log(1), log(1- 0.5), log(1 - 0.75)), 
                         labels=c(">100%\ngain", "0%\n(no change)", "50%\ndecline", ">75%\ndecline"),
                         limits = c(log(1 - 0.75), log(2))) +
  theme_void() +
  coord_equal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.3, 'cm'),
    legend.key.height = unit(0.15, 'cm'),
    panel.spacing = unit(-0.2, "lines")
  )
dev.off()



continent_rel_df <- grid_rel %>%
  group_by(continent) %>% 
  summarize(expn_nodes_mean = mean(abs_expn_nodes, na.rm = T),
            expn_nodes_hi = se_fun(abs_expn_nodes, dir = "hi"),
            expn_nodes_lo = se_fun(abs_expn_nodes, dir = "lo"),
            expn_links_mean = mean(abs_expn_links, na.rm = T),
            expn_links_hi = se_fun(abs_expn_links, dir = "hi"),
            expn_links_lo = se_fun(abs_expn_links, dir = "lo"),
            
            
            cex_nodes_mean = mean(abs_cex_nodes, na.rm = T),
            cex_nodes_hi = se_fun(abs_cex_nodes, dir = "hi"),
            cex_nodes_lo = se_fun(abs_cex_nodes, dir = "lo"),
            cex_links_mean = mean(abs_cex_links, na.rm = T),
            cex_links_hi = se_fun(abs_cex_links, dir = "hi"),
            cex_links_lo = se_fun(abs_cex_links, dir = "lo"),
            
            
            ec_nodes_mean = mean(abs_ec_nodes, na.rm = T),
            ec_nodes_hi = se_fun(abs_ec_nodes, dir = "hi"),
            ec_nodes_lo = se_fun(abs_ec_nodes, dir = "lo"),
            ec_links_mean = mean(abs_ec_links, na.rm = T),
            ec_links_hi = se_fun(abs_ec_links, dir = "hi"),
            ec_links_lo = se_fun(abs_ec_links, dir = "lo")) %>% 
  filter(!continent %in% c("Oceanic_main", "Caribbean") & !is.na(continent))

range(continent_rel_df$cex_nodes_mean)
range(continent_rel_df$cex_links_mean)

range(continent_rel_df$ec_nodes_mean)
range(continent_rel_df$ec_links_mean)


grid_rel %>%
  dplyr::na_if("NaN") %>% 
  filter(n_nodes_pres_nat > 0) %>% 
  summarize(expn_nodes_mean = mean(abs_expn_nodes, na.rm = T),
            expn_nodes_hi = se_fun(abs_expn_nodes, dir = "hi"),
            expn_nodes_lo = se_fun(abs_expn_nodes, dir = "lo"),
            expn_links_mean = mean(abs_expn_links, na.rm = T),
            expn_links_hi = se_fun(abs_expn_links, dir = "hi"),
            expn_links_lo = se_fun(abs_expn_links, dir = "lo"))


grid_rel %>%   
  dplyr::na_if("NaN") %>% 
  filter(n_nodes_extinct_only > 0) %>% 
  summarize(cex_nodes_mean = mean(abs_cex_nodes, na.rm = T),
            cex_nodes_hi = se_fun(abs_cex_nodes, dir = "hi"),
            cex_nodes_lo = se_fun(abs_cex_nodes, dir = "lo"),
            cex_links_mean = mean(abs_cex_links, na.rm = T),
            cex_links_hi = se_fun(abs_cex_links, dir = "hi"),
            cex_links_lo = se_fun(abs_cex_links, dir = "lo"))
            
grid_rel %>%     
  summarize(ec_nodes_mean = mean(abs_ec_nodes, na.rm = T),
            ec_nodes_hi = se_fun(abs_ec_nodes, dir = "hi"),
            ec_nodes_lo = se_fun(abs_ec_nodes, dir = "lo"),
            ec_links_mean = mean(abs_ec_links, na.rm = T),
            ec_links_hi = se_fun(abs_ec_links, dir = "hi"),
            ec_links_lo = se_fun(abs_ec_links, dir = "lo"))


# Percent decline attribution ---------------------------

cell <- 1:(142*360) %>% rev()
grid_long <- left_join(tibble(cell = cell), grid_web_metrics) %>% 
  mutate(x = (cell %% 360),
         y = 360 - floor(cell/360)) %>% 
  left_join(continent_cells)

grid_long_extinction_only <- past_metrics %>% 
  tibble() %>% 
  filter(focal_years == 0) %>% 
  rename(n_nodes_extinct_only = n_nodes,
         n_links_extinct_only = n_links) %>% 
  select(cell, n_nodes_extinct_only, n_links_extinct_only)


grid_long <- grid_long %>% 
  left_join(grid_long_extinction_only) %>% 
  # Get absolute declines percent declines relative to the present natural
  mutate(abs_expn_nodes = (n_nodes_extinct_only - n_nodes_pres_nat) / n_nodes_pres_nat,
         abs_expn_links = (n_links_extinct_only - n_links_pres_nat) / n_links_pres_nat,
         abs_expn_ld = ((n_links_extinct_only/n_nodes_extinct_only) - (n_links_pres_nat/n_nodes_pres_nat)) / (n_links_pres_nat/n_nodes_pres_nat),
         
         abs_cpn_nodes = (n_nodes_current - n_nodes_pres_nat) / n_nodes_pres_nat,
         abs_cpn_links = (n_links_current - n_links_pres_nat) / n_links_pres_nat,
         abs_cpn_ld = ((n_links_current/n_nodes_current) - (n_links_pres_nat/n_nodes_pres_nat)) / (n_links_pres_nat/n_nodes_pres_nat),
         
         abs_epn_nodes = (n_nodes_no_endangered - n_nodes_pres_nat) / n_nodes_pres_nat,
         abs_epn_links = (n_links_no_endangered - n_links_pres_nat) / n_links_pres_nat,
         abs_epn_ld = ((n_links_no_endangered/n_nodes_no_endangered) - (n_links_pres_nat/n_nodes_pres_nat)) / (n_links_pres_nat/n_nodes_pres_nat),
  )




# Get continent level summaries

continent_df <- grid_long %>%
  group_by(continent) %>% 
  summarize(expn_nodes_mean = mean(abs_expn_nodes, na.rm = T),
            expn_nodes_hi = se_fun(abs_expn_nodes, dir = "hi"),
            expn_nodes_lo = se_fun(abs_expn_nodes, dir = "lo"),
            expn_links_mean = mean(abs_expn_links, na.rm = T),
            expn_links_hi = se_fun(abs_expn_links, dir = "hi"),
            expn_links_lo = se_fun(abs_expn_links, dir = "lo"),
            
            
            cpn_nodes_mean = mean(abs_cpn_nodes, na.rm = T),
            cpn_nodes_hi = se_fun(abs_cpn_nodes, dir = "hi"),
            cpn_nodes_lo = se_fun(abs_cpn_nodes, dir = "lo"),
            cpn_links_mean = mean(abs_cpn_links, na.rm = T),
            cpn_links_hi = se_fun(abs_cpn_links, dir = "hi"),
            cpn_links_lo = se_fun(abs_cpn_links, dir = "lo"),
            
            
            epn_nodes_mean = mean(abs_epn_nodes, na.rm = T),
            epn_nodes_hi = se_fun(abs_epn_nodes, dir = "hi"),
            epn_nodes_lo = se_fun(abs_epn_nodes, dir = "lo"),
            epn_links_mean = mean(abs_epn_links, na.rm = T),
            epn_links_hi = se_fun(abs_epn_links, dir = "hi"),
            epn_links_lo = se_fun(abs_epn_links, dir = "lo")) %>% 
  filter(!continent %in% c("Oceanic_main", "Caribbean") & !is.na(continent))

apply(continent_df, 2, range)

range(continent_df$cpn_nodes_mean)
range(continent_df$cpn_links_mean)


global_df <- grid_long %>%
  dplyr::na_if("NaN") %>% 
  filter(!is.na(n_nodes_pres_nat)) %>% 
  filter(n_nodes_pres_nat > 0) %>% 
  summarize(expn_nodes_mean = mean(abs_expn_nodes, na.rm = T),
            expn_nodes_hi = se_fun(abs_expn_nodes, dir = "hi"),
            expn_nodes_lo = se_fun(abs_expn_nodes, dir = "lo"),
            expn_links_mean = mean(abs_expn_links, na.rm = T),
            expn_links_hi = se_fun(abs_expn_links, dir = "hi"),
            expn_links_lo = se_fun(abs_expn_links, dir = "lo"),
            
            
            cpn_nodes_mean = mean(abs_cpn_nodes, na.rm = T),
            cpn_nodes_hi = se_fun(abs_cpn_nodes, dir = "hi"),
            cpn_nodes_lo = se_fun(abs_cpn_nodes, dir = "lo"),
            cpn_links_mean = mean(abs_cpn_links, na.rm = T),
            cpn_links_hi = se_fun(abs_cpn_links, dir = "hi"),
            cpn_links_lo = se_fun(abs_cpn_links, dir = "lo"),
            
            
            epn_nodes_mean = mean(abs_epn_nodes, na.rm = T),
            epn_nodes_hi = se_fun(abs_epn_nodes, dir = "hi"),
            epn_nodes_lo = se_fun(abs_epn_nodes, dir = "lo"),
            epn_links_mean = mean(abs_epn_links, na.rm = T),
            epn_links_hi = se_fun(abs_epn_links, dir = "hi"),
            epn_links_lo = se_fun(abs_epn_links, dir = "lo"))

global_df$expn_nodes_mean
global_df$expn_links_mean

global_df$cpn_nodes_mean
global_df$cpn_links_mean



pdf(file = "figures/extinction vs range loss.pdf", width = 3.15, height = 2.75)

par(mar = c(1, 4, 3, 3))

link_col <- rgb(217,95,2, maxColorValue = 255)
link_col_poly <- rgb(217,95,2, maxColorValue = 255, alpha = 175)
link_col_poly2 <- rgb(217,95,2, maxColorValue = 255, alpha = 100)
node_col <- rgb(117,112,179, maxColorValue = 255)
node_col_poly <- rgb(117,112,179, maxColorValue = 255, alpha = 175)
node_col_poly2 <- rgb(117,112,179, maxColorValue = 255, alpha = 100)
null_col <- rgb(27,158,119, maxColorValue = 255)
null_col_poly <- rgb(27,158,119, maxColorValue = 255, alpha = 150)

xx <- seq(1, 16, length.out = 6)
plot(NA,
     xlim = c(min(xx) - 1, max(xx) + 1),
     ylim = c(-1.1, 0),
     frame = F,
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "")
mtext("Percent change", side = 2, line = 3, cex = 1)

#axis(3, at = xx, labels = NA)
region_labels <- c("Africa", "Australia","Eurasia", "Madagascar",
                   "N. America", "S. America")
text(region_labels, x = xx-1.5, y = 0.04, pos = 4, srt = 35, xpd = T,
     cex = 0.85)
axis(2, at = c(0, -0.2, -0.4, -0.6, -0.8, -1), #log(c(1, 0.5, 0.2, 0.05)), 
     labels = c("0%", "-20%", "-40%", "-60%", "-80%", "-100%"), las = 1,
     cex.axis = 0.75)
seg_lwd <- 6.5

x_diff <- 0.4
segments(x0 = xx - x_diff,
         y0 = 0,
         y1 = continent_df$expn_nodes_mean,
         lend = "butt",
         lwd = seg_lwd,
         col = node_col)

segments(x0 = xx - x_diff,
         y0 = continent_df$expn_nodes_mean,
         y1 = continent_df$cpn_nodes_mean,
         lend = "butt",
         lwd = seg_lwd,
         col = node_col_poly)

segments(x0 = xx - x_diff,
         y0 = continent_df$cpn_nodes_mean,
         y1 = continent_df$epn_nodes_mean,
         lend = "butt",
         lwd = seg_lwd,
         col = node_col_poly2)

arrow_length <- 0.025

arrows(x0 = xx - x_diff,
       y0 = continent_df$expn_nodes_lo,
       y1 = continent_df$expn_nodes_hi,
       angle = 90,
       code = 3,
       length = arrow_length,
       col = "grey")

arrows(x0 = xx - x_diff,
       y0 = continent_df$cpn_nodes_lo,
       y1 = continent_df$cpn_nodes_hi,
       angle = 90,
       code = 3,
       length = arrow_length,
       col = "grey")

arrows(x0 = xx - x_diff,
       y0 = continent_df$epn_nodes_lo,
       y1 = continent_df$epn_nodes_hi,
       angle = 90,
       code = 3,
       length = arrow_length,
       col = "grey")


segments(x0 = xx + x_diff,
         y0 = 0,
         y1 = continent_df$expn_links_mean,
         lend = "butt",
         lwd = seg_lwd,
         col = link_col)

segments(x0 = xx + x_diff,
         y0 = continent_df$expn_links_mean,
         y1 = continent_df$cpn_links_mean,
         lend = "butt",
         lwd = seg_lwd,
         col = link_col_poly)

segments(x0 = xx + x_diff,
         y0 = continent_df$cpn_links_mean,
         y1 = continent_df$epn_links_mean,
         lend = "butt",
         lwd = seg_lwd,
         col = link_col_poly2)

arrows(x0 = xx + x_diff,
       y0 = continent_df$expn_links_lo,
       y1 = continent_df$expn_links_hi,
       angle = 90,
       code = 3,
       length = arrow_length,
       col = "grey")

arrows(x0 = xx + x_diff,
       y0 = continent_df$cpn_links_lo,
       y1 = continent_df$cpn_links_hi,
       angle = 90,
       code = 3,
       length = arrow_length,
       col = "grey")

arrows(x0 = xx + x_diff,
       y0 = continent_df$epn_links_lo,
       y1 = continent_df$epn_links_hi,
       angle = 90,
       code = 3,
       length = arrow_length,
       col = "grey")

segments(x0 = min(xx)-2, x1 = max(xx)+2, y0 = 0, lty = 2)

xx <- 13.5
xxd <- 1.8
yy <- -1-0.1
yyd <- 0.1

points(x = c(xx, xx, xx, xx + xxd, xx + xxd, xx + xxd),
       y = c(yy + yyd, yy, yy - yyd, yy + yyd, yy, yy - yyd),
       pch = 15,
       col = c(node_col,
               node_col_poly,
               node_col_poly2,
               link_col,
               link_col_poly,
               link_col_poly2),
       xpd = T,
       cex = 1.6)
text(c(xx, xx + xxd) - 1.1,
     yy + yyd + 0.06,
     pos = 4, 
     srt = 45,
     labels = c("Species", "Links"),
     xpd = T,
     cex = 0.8)
text(xx,
     c(yy, yy + yyd, yy - yyd),
     pos = 2,
     labels = c("Range change", "Extinction", "Endangered spp. extinction"),
     cex = 0.8, xpd = T)


dev.off()



#



# Supplementary figures

# Figure out the average degree for all species

degree_fun <- function(x) x %>% unlist() %>% table() %>% as.data.frame()

degree_df <- lapply(web_pres_nat[cells_with_webs], degree_fun) %>% 
  do.call(rbind, .)

colnames(degree_df) <- c("binomial", "degree")

head(degree_df)

degree_mean <- degree_df %>% 
  tibble() %>% 
  group_by(binomial) %>% 
  summarize(degree_mean = mean(degree))

# Get estimates of range change
load("data/m.mamm.current.RData")
range_change <- rowSums(m.mamm.current) / rowSums(m.mamm.pres.nat) %>% as.data.frame()
range_change$binomial <- rownames(range_change)
colnames(range_change) <- c("range_change", "binomial")


# Lastly get 

# Also want a no vulnerable and endangered species version
endangered_categories <- c("CR", "EN", "EP", "EW", "EX", "VU", "CR (PE)")
# Need to pull back in COMBINE for this to reconcile names...
temp <- tempfile()
download.file("https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3344&file=ecy3344-sup-0001-DataS1.zip",temp)
trait_data <- read.table(unz(temp, "COMBINE_archives/trait_data_imputed.csv"), sep = ",", header = T) %>% tibble()
unlink(temp)
trait_data <- trait_data %>% 
  mutate(iucn2020_binomial = word(iucn2020_binomial, 1, 2))
iucn_categories <- read.csv("data/mamm.iucn.categories.csv", header = T)[,-1] %>% tibble()

# Get vector of endangered species names while accounting for different names across iucn and phylacine
endangered_spp <- c(left_join(trait_data, iucn_categories, by = c("iucn2020_binomial" = "iucn.bin")) %>%  
                      filter(iucn.category %in% endangered_categories) %>% 
                      select("iucn2020_binomial", "phylacine_binomial") %>%
                      unlist(),
                    left_join(trait_data, iucn_categories, by = c("phylacine_binomial" = "iucn.bin")) %>%  
                      filter(iucn.category %in% endangered_categories) %>% 
                      select("iucn2020_binomial", "phylacine_binomial") %>%
                      unlist()) %>% 
  unique()



#

degree_mean <- degree_mean %>% 
  left_join(select(te_dates, binomial, te)) %>% 
  left_join(range_change) %>% 
  mutate(extinct = ifelse(is.na(te), "Extant", "Extinct")) %>% 
  mutate(endangered = ifelse(binomial %in% endangered_spp, "Endangered", "Not endangered"))

degree_mean$extinct <- factor(degree_mean$extinct, levels = c("Extinct", "Extant"))

# 
# p1 <- degree_mean %>% 
#   ggplot(aes(degree_mean, fill = extinct, color = extinct)) +
#   geom_density(alpha = 0.1, adjust = 1.5) +
#   scale_x_log10() + 
#   theme_classic() +
#   xlab("Species degree") +
#   ylab("Density") +
#   theme(legend.position = c(0.7, 0.6)) +
#   labs(fill='', color = '') +
#   theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")) +
#   annotate(geom="text", x=1, y=3.4, label="a", fontface = 2)
# 
# p2 <- degree_mean %>% 
#   filter(extinct != "Extinct") %>% 
#   ggplot(aes(degree_mean, fill = endangered, color = endangered)) +
#   geom_density(alpha = 0.1, adjust = 1.5) +
#   scale_x_log10() +
#   theme_classic() +
#   xlab("Species degree") +
#   ylab("Density") +
#   theme(legend.position = c(0.7, 0.6)) +
#   labs(fill='', color = '')  +
#   theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")) +
#   annotate(geom="text", x=1, y=4.1, label="c", fontface = 2)
# 
# set.seed(44)
# p3 <- degree_mean %>% 
#   filter(extinct != "Extinct") %>% 
#   ggplot(aes((range_change - 1) * -100, degree_mean)) +
#   scale_y_log10() +
#   #geom_point() +
#   geom_jitter(width = 2, height = 0.05, color = "grey", size = 0.5) +
#   theme_classic() +
#   geom_smooth(method = lm) +
#   xlim(0, 100) +
#   xlab("Percent range loss") +
#   ylab("Species degree") +
#   theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")) +
#   annotate(geom="text", x=0, y=120, label="b", fontface = 2)
# 
# pdf(file = "degree plots.pdf", width = 7.2, height = 2.5)
# gridExtra::grid.arrange(p1, p3, p2, nrow = 1)
# dev.off()

xx <- 0.5

p1 <- degree_mean %>% 
  ggplot(aes(y = degree_mean, x = extinct, fill = extinct, color = extinct)) +
  scale_y_log10() + 
  geom_violin(alpha = 0.1) +
  theme_classic() +
  ylab("Species degree") +
  xlab("") +
  xlim("Extant", "Extinct") +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")) +
  annotate(geom="text", x = xx, y=100, label="a", fontface = 2)


p2 <- degree_mean %>% 
  filter(extinct != "Extinct") %>% 
  ggplot(aes(y = degree_mean, x = endangered, fill = endangered, color = endangered)) +
  scale_y_log10() + 
  geom_violin(alpha = 0.1) +  
  # scale_x_discrete("endangered",
  #                  labels = c("Not endangered" = "Not\nendangered",
  #                             "Endangered" = "Endangered")) +
  theme_classic() +
  ylab("") +
  xlab("") +
  xlim("Not endangered", "Endangered") +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")) +
  annotate(geom="text", x = xx, y=100, label="c", fontface = 2)

set.seed(44)
p3 <- degree_mean %>% 
  filter(extinct != "Extinct") %>% 
  ggplot(aes((range_change - 1) * -100, degree_mean)) +
  scale_y_log10() +
  #geom_point() +
  geom_jitter(width = 2, height = 0.05, color = rgb(0,0,0, 0.1), size = 0.5) +
  theme_classic() +
  geom_smooth(method = lm, se = F) +
  xlim(0, 100) +
  xlab("Percent range loss") +
  ylab("") +
  theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm")) +
  annotate(geom="text", x=0, y=120, label="b", fontface = 2)

pdf(file = "figures/degree plots.pdf", width = 7.2, height = 2.5)
gridExtra::grid.arrange(p1, p3, p2, nrow = 1)
dev.off()


# Statistical models for differences in species degree

extinction_degree_mod <- lm(log(degree_mean) ~ extinct, data = degree_mean)
summary(extinction_degree_mod)
exp(coef(extinction_degree_mod)[1]) / exp(coef(extinction_degree_mod)[1] + coef(extinction_degree_mod)[2])

range_change_degree_mod <- lm(log(degree_mean) ~ range_change, data = filter(degree_mean, extinct != "Extinct"))
summary(range_change_degree_mod)
exp(coef(range_change_degree_mod)[1]) / exp(coef(range_change_degree_mod)[1] + coef(range_change_degree_mod)[2])

endangered_degree_mod <- lm(log(degree_mean) ~ endangered, data = filter(degree_mean, extinct != "Extinct"))
summary(endangered_degree_mod)
exp(coef(endangered_degree_mod)[1]) / exp(coef(endangered_degree_mod)[1] + coef(endangered_degree_mod)[2])






#


















# # I THINK THIS IS ALL OBSOLETE
# 
# # Want to get a tibble with change (delta) between 4 categories
# # c = current # pn = present natural # ex = extinctions alone # e = endangered species extinction
# cell <- 1:(142*360) %>% rev()
# grid_long <- left_join(tibble(cell = cell), grid_web_metrics) %>% 
#   mutate(delta_cpn_nodes = n_nodes_current / n_nodes_pres_nat,
#          delta_cpn_links = n_links_current / n_links_pres_nat,
#          delta_epn_nodes = n_nodes_no_endangered / n_nodes_pres_nat,
#          delta_epn_links = n_links_no_endangered / n_links_pres_nat,
#          delta_ec_nodes = n_nodes_no_endangered / n_nodes_current,
#          delta_ec_links = n_links_no_endangered / n_links_current) %>% 
#   filter(!is.na(n_links_pres_nat)) 
# 
# 
# # Want to also have an estimate for just the change attributable to extinction alone
# grid_long_extinction_only <- prop_past_metrics %>% 
#   filter(focal_years == 0) %>% 
#   rename(delta_expn_nodes = n_nodes_change,
#          delta_expn_links = n_links_change,
#          delta_cex_nodes = n_nodes_change_current,
#          delta_cex_links = n_links_change_current) %>% # Want to rename to follow the
#   select(cell, 
#          delta_expn_nodes,
#          delta_expn_links,
#          delta_cex_nodes,
#          delta_cex_links)
# 
# # Now join this to grid_long
# grid_long <- grid_long %>% 
#   left_join(grid_long_extinction_only) %>% 
#   # Also want a few absolute (abs) rather than relative (delta) differences
#   mutate(abs_expn_nodes = 1 - delta_expn_nodes,
#          abs_expn_links = 1 - delta_expn_links,
#          
#          abs_cex_nodes = (1 - delta_cpn_nodes) - (1 - delta_expn_nodes),
#          abs_cex_links = (1 - delta_cpn_links) - (1 - delta_expn_links),
#          
#          abs_ec_nodes = (1 - delta_epn_nodes) - (1 - delta_cpn_nodes),
#          abs_ec_links = (1 - delta_epn_links) - (1 - delta_cpn_links)) %>% 
#   # And linkage density
#   mutate(abs_expn_ld = abs_expn_links / abs_expn_nodes,
#          abs_cex_ld = abs_cex_links / abs_cex_nodes,
#          abs_ec_ld = abs_ec_links / abs_ec_nodes)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# grid_longer$metric <- factor(grid_longer$metric, levels = c("delta_cpn_nodes",
#                                                             "delta_cpn_links",
#                                                             "delta_cpn_linkage_density",
#                                                             "delta_cpn_connectance",
#                                                             "delta_cpn_predator_prey_ratio"), 
#                              labels = c("Number of species",
#                                         "Number of food web links",
#                                         "Linkage density (mean links per species)",
#                                         "Connectance",
#                                         "Predator prey ratio"))
# 
# pdf(file = "change maps.pdf", width = 7.2, height = 4.5)
# grid_longer %>% 
#   mutate(value = ifelse(value > 2, 2, value)) %>%
#   mutate(value = ifelse(value < 0.25, 0.25, value)) %>%
#   ggplot(aes(x, y)) +
#   geom_point(aes(colour = log(value)), size = 0.005, shape = 15) +
#   facet_wrap(vars(metric), nrow=3) +
#   scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0,
#                          na.value = "white",
#                          name = "Percent change\n(current / natural)\n",
#                          breaks=c(log(2), log(1), log(1- 0.5), log(1 - 0.75)), 
#                          labels=c(">200% gain", "0% (no change)", "50% decline", ">75% decline"),
#                          limits = c(log(1 - 0.75), log(2))) +
#   theme_void() +
#   coord_equal() +
#   theme(legend.position = c(0.86, 0.04),
#         legend.justification = c(1, 0),
#         legend.key.height = unit(0.4, 'cm'),
#         legend.title = element_text(size=10))
# dev.off()
# 
# 
# 
# 
# 
# 
# link_col <- rgb(217,95,2, maxColorValue = 255)
# link_col_poly <- rgb(217,95,2, maxColorValue = 255, alpha = 150)
# link_col_poly2 <- rgb(217,95,2, maxColorValue = 255, alpha = 100)
# node_col <- rgb(117,112,179, maxColorValue = 255)
# node_col_poly <- rgb(117,112,179, maxColorValue = 255, alpha = 150)
# node_col_poly2 <- rgb(117,112,179, maxColorValue = 255, alpha = 100)
# null_col <- rgb(27,158,119, maxColorValue = 255)
# null_col_poly <- rgb(27,158,119, maxColorValue = 255, alpha = 150)
# 
# pdf(file = "extinction vs range loss.pdf", width = 4.5, height = 4.5)
# 
# par(mar = c(1, 5, 4, 2))
# 
# continent_df <- prop_past_metrics %>% filter(focal_years == 0) %>%
#   group_by(continent) %>% 
#   summarize(extinction_nodes_mean = mean(n_nodes_change, na.rm = T),
#             extinction_nodes_hi = se_fun(n_nodes_change, dir = "hi"),
#             extinction_nodes_lo = se_fun(n_nodes_change, dir = "lo"),
#             extinction_links_mean = mean(n_links_change, na.rm = T),
#             extinction_links_hi = se_fun(n_links_change, dir = "hi"),
#             extinction_links_lo = se_fun(n_links_change, dir = "lo"),
#             
#             current_nodes_mean = mean(n_nodes_change_current, na.rm = T),
#             current_nodes_hi = se_fun(n_nodes_change_current, dir = "hi"),
#             current_nodes_lo = se_fun(n_nodes_change_current, dir = "lo"),
#             current_links_mean = mean(n_links_change_current, na.rm = T),
#             current_links_hi = se_fun(n_links_change_current, dir = "hi"),
#             current_links_lo = se_fun(n_links_change_current, dir = "lo")) %>% 
#   filter(!continent %in% c("Oceanic_main", "Caribbean") & !is.na(continent))
# 
# 
# xx <- seq(1, 16, length.out = 6)
# plot(NA,
#      xlim = c(min(xx) - 1, max(xx) + 1),
#      ylim = c(0.0, 1),
#      frame = F,
#      xaxt = "n",
#      yaxt = "n",
#      xlab = "",
#      ylab = "")
# mtext("Percent change", side = 2, line = 3.5, cex = 1.2)
# 
# #axis(3, at = xx, labels = NA)
# region_labels <- c("Africa", "Australia","Eurasia", "Madagascar",
#                    "N. America", "S. America")
# text(region_labels, x = xx-1, y = 1.02, pos = 4, srt = 45, xpd = T)
# axis(2, at = c(1, 0.8, 0.6, 0.4, 0.2, 0), #log(c(1, 0.5, 0.2, 0.05)), 
#      labels = c("0%", "-20%", "-40%", "-60%", "-80%", "-100%"), las = 1)
# seg_lwd <- 11
# 
# x_diff <- 0.4
# segments(x0 = xx - x_diff,
#          y0 = 1,
#          y1 = continent_df$extinction_nodes_mean,
#          lend = "butt",
#          lwd = seg_lwd,
#          col = node_col)
# 
# segments(x0 = xx - x_diff,
#          y0 = continent_df$extinction_nodes_mean,
#          y1 = continent_df$current_nodes_mean,
#          lend = "butt",
#          lwd = seg_lwd,
#          col = node_col_poly)
# 
# arrows(x0 = xx - x_diff,
#        y0 = continent_df$extinction_nodes_lo,
#        y1 = continent_df$extinction_nodes_hi,
#        angle = 90,
#        code = 3,
#        length = 0.04,
#        col = "grey")
# 
# arrows(x0 = xx - x_diff,
#        y0 = continent_df$current_nodes_lo,
#        y1 = continent_df$current_nodes_hi,
#        angle = 90,
#        code = 3,
#        length = 0.04,
#        col = "grey")
# 
# segments(x0 = xx + x_diff,
#          y0 = 1,
#          y1 = continent_df$extinction_links_mean,
#          lend = "butt",
#          lwd = seg_lwd,
#          col = link_col)
# 
# segments(x0 = xx + x_diff,
#          y0 = continent_df$extinction_links_mean,
#          y1 = continent_df$current_links_mean,
#          lend = "butt",
#          lwd = seg_lwd,
#          col = link_col_poly)
# 
# arrows(x0 = xx + x_diff,
#        y0 = continent_df$extinction_links_lo,
#        y1 = continent_df$extinction_links_hi,
#        angle = 90,
#        code = 3,
#        length = 0.04,
#        col = "grey")
# 
# arrows(x0 = xx + x_diff,
#        y0 = continent_df$current_links_lo,
#        y1 = continent_df$current_links_hi,
#        angle = 90,
#        code = 3,
#        length = 0.04,
#        col = "grey")
# 
# segments(x0 = min(xx)-2, x1 = max(xx)+2, y0 = 1, lty = 2)
# 
# xx <- 14.5
# xxd <- 1.2
# yy <- 0.1
# yyd <- 0.06
# 
# points(x = c(xx, xx + xxd, xx + xxd, xx),
#        y = c(yy, yy + yyd, yy, yy + yyd),
#        pch = 15,
#        col = c(node_col_poly,
#                link_col,
#                link_col_poly,
#                node_col),
#        cex = 2)
# text(c(xx, xx + xxd) - 0.7,
#      yy + yyd + 0.04,
#      pos = 4, 
#      srt = 45,
#      labels = c("Species", "Links"),
#      cex = 0.8)
# text(xx,
#      c(yy, yy + yyd),
#      pos = 2,
#      labels = c("Range change", "Extinction"),
#      cex = 0.8)
# 
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# pdf(file = "endangerment maps.pdf", width = 3.5, height = 3.5)
# 
# 
# endangered_df <- grid_web_metrics %>% 
#   dplyr::na_if("NaN") %>% 
#   left_join(continent_cells) %>% 
#   mutate(nodes_endangered = n_nodes_no_endangered / n_nodes_current,
#          links_endangered = n_links_no_endangered / n_links_current) %>% 
#   tibble() %>% 
#   mutate(x = (cell %% 360),
#          y = 360 - floor(cell/360))
# 
# 
# endangered_longer <- endangered_df %>% 
#   pivot_longer(c(nodes_endangered, links_endangered), names_to = "metric", values_to = "value")
# 
# endangered_longer$metric <- factor(endangered_longer$metric, levels = c("nodes_endangered",
#                                                                         "links_endangered"), 
#                                    labels = c("Number of species",
#                                               "Number of food web links"))
# 
# endangered_longer %>% 
#   mutate(value = ifelse(value < 0.25, 0.25, value)) %>% 
#   ggplot(aes(x, y)) +
#   geom_point(aes(colour = log(value)), size = 0.1, shape = 15) +
#   facet_wrap(vars(metric), nrow=3) +
#   scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0,
#                          na.value = "white",
#                          name = "Food web\nendangerment\n",
#                          breaks=c(log(1), log(1- 0.5), log(1 - 0.75)), 
#                          labels=c("0%", "50%", ">75%"),
#                          limits = c(log(1 - 0.75), log(1))) +
#   theme_void() +
#   coord_equal() +
#   theme(legend.position="bottom")
# 
# 
# dev.off()
# 
# 
# pdf(file = "endangerment bars.pdf", width = 3.5, height = 3.5)
# 
# par(mar = c(4, 5, 0, 0.5))
# 
# endangered_continent_df <- endangered_df %>%
#   group_by(continent) %>% 
#   summarize(extinction_nodes_mean = mean(nodes_endangered, na.rm = T),
#             extinction_nodes_hi = se_fun(nodes_endangered, dir = "hi"),
#             extinction_nodes_lo = se_fun(nodes_endangered, dir = "lo"),
#             extinction_links_mean = mean(links_endangered, na.rm = T),
#             extinction_links_hi = se_fun(links_endangered, dir = "hi"),
#             extinction_links_lo = se_fun(links_endangered, dir = "lo")) %>% 
#   filter(!continent %in% c("Oceanic_main", "Caribbean") & !is.na(continent))
# 
# 
# xx <- seq(1, 16, length.out = 6)
# plot(NA,
#      xlim = c(min(xx) - 1, max(xx) + 1),
#      ylim = c(1,0.1),
#      frame = F,
#      xaxt = "n",
#      yaxt = "n",
#      xlab = "",
#      ylab = "")
# mtext("Food web endangerment", side = 2, line = 3.5, cex = 1.2)
# 
# #axis(3, at = xx, labels = NA)
# region_labels <- c("Africa", "Australia","Eurasia", "Madagascar",
#                    "N. America", "S. America")
# text(region_labels, x = xx+0.5, y = 1.02, pos = 2, srt = 45, xpd = T)
# axis(2, at = c(1, 0.8, 0.6, 0.4, 0.2), #log(c(1, 0.5, 0.2, 0.05)), 
#      labels = c("0%", "20%", "40%", "60%", "80%"), las = 1)
# seg_lwd <- 9
# 
# x_diff <- 0.4
# segments(x0 = xx - x_diff,
#          y0 = 1,
#          y1 = endangered_continent_df$extinction_nodes_mean,
#          lend = "butt",
#          lwd = seg_lwd,
#          col = node_col)
# 
# 
# arrows(x0 = xx - x_diff,
#        y0 = endangered_continent_df$extinction_nodes_lo,
#        y1 = endangered_continent_df$extinction_nodes_hi,
#        angle = 90,
#        code = 3,
#        length = 0.03,
#        col = "grey")
# 
# segments(x0 = xx + x_diff,
#          y0 = 1,
#          y1 = endangered_continent_df$extinction_links_mean,
#          lend = "butt",
#          lwd = seg_lwd,
#          col = link_col)
# 
# 
# arrows(x0 = xx + x_diff,
#        y0 = endangered_continent_df$extinction_links_lo,
#        y1 = endangered_continent_df$extinction_links_hi,
#        angle = 90,
#        code = 3,
#        length = 0.03,
#        col = "grey")
# 
# 
# segments(x0 = min(xx)-2, x1 = max(xx)+2, y0 = 1, lty = 2)
# xx <- 6
# yy <- 0.3
# yyd <- 0.07
# points(x = c(xx, xx),
#        y = c(yy, yy + yyd),
#        pch = 15,
#        col = c(node_col,
#                link_col),
#        cex = 2)
# text(xx,
#      c(yy, yy + yyd),
#      pos = 2,
#      labels = c("Species", "Links"),
#      cex = 0.8)
# 
# dev.off()
# 
# 








