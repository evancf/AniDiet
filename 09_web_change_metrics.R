library("tidyverse")


# 
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

# Put these together
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

load("hindcast_webs.RData")


# Load trait data ------------------------------

impute_trait_data <- read.csv("impute_trait_data.csv")[,-1] %>% tibble()


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

te_dates$binomial <- plyr::revalue(te_dates$binomial, 
                                   c("Alces scotti" = "Cervalces scotti",
                                     "Candiacervus SpII" = "Candiacervus spII",
                                     #"Caprini indet" = "",
                                     "Dicroceros sp" = "Dicroceros spA",
                                     "Geocapromys SP_A" = "Geocapromys spA",
                                     #"HexolobodontinaeGen_NOW Sp_NOW" = "",
                                     "Homo denisovans" = "Homo spDenisova",
                                     #"Hydrodamalis gigas" = "",
                                     "Megaoryzomys Sp_Now" = "Megaoryzomys spA",
                                     #"Neomonachus tropicalis" = "",
                                     "Nesophontes SP_A" = "Nesophontes spA",
                                     "Nesoryzomys Sp_B" = "Nesoryzomys spB",
                                     "Nesoryzomys Sp_C" = "Nesoryzomys spC",
                                     "Nesoryzomys Sp_D" = "Nesoryzomys spD",
                                     "Nothrotheriops shastense" = "Nothrotheriops shastensis",
                                     "Pachyarmaterium brasiliense" = "Pachyarmatherium brasiliense",
                                     "Peroryctes SP_NOW" = "Peroryctes spA",
                                     #"Peroryctinae GEN_NOW" = "",
                                     "Tapirus copei" = "Tapirus merriami",
                                     "Valgipes deformis" = "Valgipes bucklandi"#,
                                     #"Zalophus japonicus" = ""
                                   ))


# Decide on some relevant focal years (rather than every year)

(seq(0, (126000)^(1/3), length.out = 20))^3
focal_years <- c(0, 400, 3000, 10000, 23000, 46000, 80000, 126000)
focal_years <- (seq(0, (126000)^(1/3), length.out = 40))^3 %>% round()

#plot(y = rep(1, length(focal_years)), x = (focal_years ^ (1/3)))


# Get metrics for the cells for which we were actually able to get networks.
# In some it wasn't possible because there aren't terrestrial mammalian carnivore
# interactions (e.g., oceanic islands)

cells_with_webs <- which(!lapply(web_pres_nat, is.null) %>% unlist())

n_links_fun <- function(x) dim(x)[1]
n_nodes_fun <- function(x) length(unique(unlist(x)))


grid_web_metrics <- tibble(cell = cells_with_webs,
                           n_links_pres_nat = lapply(web_pres_nat[cells_with_webs], n_links_fun) %>% unlist,
                           n_nodes_pres_nat = lapply(web_pres_nat[cells_with_webs], n_nodes_fun) %>% unlist,
                           #mean_chain_length_pres_nat = NA,
                           n_links_current = lapply(web_current[cells_with_webs], n_links_fun) %>% unlist,
                           n_nodes_current = lapply(web_current[cells_with_webs], n_nodes_fun) %>% unlist,
                           #mean_chain_length_current = NA
) %>% as.data.frame()



# Decide how many of 100 extinction years to sample
te_samp <- 1 #paste("te", 1, sep = ".")

# Get matrix of link and node metrics

past_metrics <- expand_grid(cells_with_webs, focal_years, te_samp) %>% 
  rename(cell = cells_with_webs) %>% 
  mutate(n_links = NA, n_nodes = NA,
         n_links_null = NA, n_nodes_null = NA)

for(i in cells_with_webs){
  for(k in te_samp){ # Switching the indexing a little
    
    # Create vectors to be populated with a record of the mammal species
    # that went extinct in a given cell, or simulated to go extinct
    # under a random simulation
    true_extinct <- c()
    null_extinct <- c()
    
    for(j in rev(focal_years)){
      
      ind <- which(past_metrics$cell == i & past_metrics$focal_years == j & past_metrics$te_samp == k)
      
      te_cell_samp <- paste("te", sample(1:99, 1), sep = ".")
      
      # Determine if the consumer and resource are absent at the given date
      c_bool <-  !web_pres_nat[[i]][,1] %in% te_dates$binomial[te_dates[,te_cell_samp] >= j]
      r_bool <-  !web_pres_nat[[i]][,2] %in% te_dates$binomial[te_dates[,te_cell_samp] >= j]
      
      past_web <- web_pres_nat[[i]][c_bool & r_bool,]
      
      past_metrics$n_links[ind] <- past_web %>% n_links_fun()
      past_metrics$n_nodes[ind] <- past_web %>% n_nodes_fun()
      
      # Get record of the now extinct species
      new_true_extinct <- tibble(new_true_extinct = c(web_pres_nat[[i]][,1][!c_bool], 
                                                      web_pres_nat[[i]][,2][!r_bool]) %>% unique()) %>% 
        filter(!new_true_extinct %in% true_extinct)
      
      true_extinct <- c(true_extinct, new_true_extinct$new_true_extinct)
      
      # For null simulation, choose random species
      new_null_extinct <- unique(c(web_pres_nat[[i]][,1], web_pres_nat[[i]][,2])) %>% 
        sample(dim(new_true_extinct)[1])
      
      null_extinct <- c(null_extinct, new_null_extinct)
      
      
      # Perform similar calculations for web with randomized extinctions
      c_bool_null <-  !web_pres_nat[[i]][,1] %in% null_extinct
      r_bool_null <-  !web_pres_nat[[i]][,2] %in% null_extinct
      
      past_web_null <- web_pres_nat[[i]][c_bool_null & r_bool_null,]
      
      past_metrics$n_links_null[ind] <- past_web_null %>% n_links_fun()
      past_metrics$n_nodes_null[ind] <- past_web_null %>% n_nodes_fun()
      
      
    }
  }
  print(round(i / tail(cells_with_webs, 1), digits = 3))
}

write.csv(file = "past_metrics.csv", past_metrics)

# Left join to pres_nat values

prop_past_metrics <- past_metrics %>% 
  left_join(grid_web_metrics) %>% 
  mutate(n_links_change = n_links / n_links_pres_nat,
         n_nodes_change = n_nodes / n_nodes_pres_nat,
         n_links_change_null = n_links_null / n_links_pres_nat,
         n_nodes_change_null = n_nodes_null / n_nodes_pres_nat) %>% 
  left_join(continent_cells)

prop_past_metrics <- prop_past_metrics %>%
  group_by(cell, focal_years) %>% 
  summarise(continent = first(continent),
            n_links_change = mean(n_links_change),
            n_nodes_change = mean(n_nodes_change),
            n_links_change_null = mean(n_links_change_null),
            n_nodes_change_null = mean(n_nodes_change_null),
            n_nodes_change_current = mean(n_nodes_current/n_nodes_pres_nat),
            n_links_change_current = mean(n_links_current/n_links_pres_nat))

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

pdf(file = "change over time.pdf", width = 4.5, height = 5.5)
par(mfrow = c(3,2))
counter <- 1
bottom_row_inds <- 5:6
left_col_inds <- c(1,3,5)
focal_year_labels <- c(0, 400, 3000, 10000, 23000, 46000, 80000, 126000)
focal_year_labels <- c(0, 400, 10000, 46000, 126000)

for(i in levels(prop_past_metrics$continent)){
  
  
  dat <- prop_past_metrics %>% filter(continent == i) %>% 
    mutate(y1 = n_links_change, #log(n_links_change),
           y2 = n_nodes_change, #log(n_nodes_change),
           y3 = n_links_change_null,
           y4 = n_nodes_change_null,
           x = c(focal_years ^ (1/3)))
  
  if(i %in% c("Oceanic_main", "Caribbean")) next()
  
  if(counter %in% left_col_inds){
    par(mar = c(3.8,5,1.1,0.2))
  } else{
    par(mar = c(3.8,4,1.1,1.2))
  }
  
  plot(NA, 
       #data = dat,
       xlim = (c(126000, 0) ^(1/3)),
       ylim = c(0.25, 1), #log(c(0.05, 1)),
       pch = 16,
       col = rgb(0,0,1,0.1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       frame = F,
       las = 1)
  
  text(52, 1.075, i, cex = 1.2, xpd = T, pos = 4)
  
  if(counter == 6) mtext("Thousand years before present", side = 1, line = 2.5, cex = 0.9, adj = 2.345, xpd = T)
  if(counter == 3) mtext("Percent change", side = 2, line = 3.5, cex = 0.9)
  
  if(counter %in% bottom_row_inds){
    axis(1, 
         at = rev(focal_year_labels^(1/3)), 
         labels = rev(focal_year_labels/1000))
  } else{
    axis(1, 
         at = rev(focal_year_labels^(1/3)), 
         labels = NA)
  }

  if(counter %in% left_col_inds){
    axis(2, at = c(1, 0.75, 0.5, 0.25), #log(c(1, 0.5, 0.2, 0.05)), 
         labels = c("0%", "-25%", "-50%", "-75%"), las = 1)
  } else{
    axis(2, at = c(1, 0.75, 0.5, 0.25), #log(c(1, 0.5, 0.2, 0.05)), 
         labels = NA, las = 1)
  }
    
  
  rect(xleft = (filter(arrival_kya, continent == i)$min * 1000)^(1/3),
       xright = (filter(arrival_kya, continent == i)$max * 1000)^(1/3),
       ybottom = 0.23,
       ytop = 1.02,
       border = F,
       col = "lightgrey")
  
  xx <- c(focal_years ^ (1/3))
  
  samp_mean <- function(x,i){mean(x[i])}
  
  se_fun <- function(x, dir = "hi"){
    if(dir == "hi"){
      return(mean(x) + sd(x)/sqrt(length(x)) * 1.96)
      #return(quantile(x, 0.025))
      #return(boot::boot(x,samp_mean,1000)$t %>% quantile(0.975))
    } 
    if(dir == "lo"){
      return(mean(x) - sd(x)/sqrt(length(x)) * 1.96)
      #return(quantile(x, 0.975))
      #return(boot::boot(x,samp_mean,1000)$t %>% quantile(0.025))
    } 
    
  } 
  yy1_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y1)) %>% ungroup() %>% select(mean) %>% unlist()
  yy1_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y1, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy1_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y1, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()

  yy2_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y2)) %>% ungroup() %>% select(mean) %>% unlist()
  yy2_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y2, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy2_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y2, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  yy3_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y3)) %>% ungroup() %>% select(mean) %>% unlist()
  yy3_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y3, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy3_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y3, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  yy4_mean <- dat %>% group_by(x) %>% summarize(mean = mean(y4)) %>% ungroup() %>% select(mean) %>% unlist()
  yy4_sehi <- dat %>% group_by(x) %>% summarize(sehi = se_fun(y4, dir = "hi")) %>% ungroup() %>% select(sehi) %>% unlist()
  yy4_selo <- dat %>% group_by(x) %>% summarize(selo = se_fun(y4, dir = "lo")) %>% ungroup() %>% select(selo) %>% unlist()
  
  
  
  link_col <- rgb(217,95,2, maxColorValue = 255)
  link_col_poly <- rgb(217,95,2, maxColorValue = 255, alpha = 200)
  node_col <- rgb(117,112,179, maxColorValue = 255)
  node_col_poly <- rgb(117,112,179, maxColorValue = 255, alpha = 200)
  null_col <- rgb(27,158,119, maxColorValue = 255)
  null_col_poly <- rgb(27,158,119, maxColorValue = 255, alpha = 200)
  
  polygon(x = c(xx, rev(xx)), y = c(yy1_sehi, rev(yy1_selo)),
          border = F,
          col = link_col_poly)
  lines(x = xx, y = yy1_mean, col = link_col)
  
  polygon(x = c(xx, rev(xx)), y = c(yy2_sehi, rev(yy2_selo)),
          border = F,
          col = node_col_poly)
  lines(x = xx, y = yy2_mean, col = node_col)
  
  polygon(x = c(xx, rev(xx)), y = c(yy3_sehi, rev(yy3_selo)),
          border = F,
          col = null_col_poly)
  lines(x = xx, y = yy3_mean, col = null_col)
  
  polygon(x = c(xx, rev(xx)), y = c(yy4_sehi, rev(yy4_selo)),
          border = F,
          col = "lightgrey")
  lines(x = xx, y = yy4_mean, col = "darkgrey")
  
  if(counter == 1){
    text(x = 52,#-2.5,
         y = c(0.28, 0.37),
         pos = 4,#2,
         font = 3,
         labels = c("Number of links",
                    "Number of species"),
         col = c(link_col, node_col))
    
    text(x = 15000^(1/3),
         y = 0.8,
         pos = 2,
         font = 3,
         labels = "Extinction\nonly")
    arrows(x0 = 15000^(1/3),
           x1 = 4000^(1/3),
           y0 = 0.8 + 0.03,
           y1 = 0.94,
           length = 0.07)

    
    text(x = 500^(1/3),
         y = 0.55,
         pos = 2,
         font = 3,
         labels = "Extinction &\nrange change")
    arrows(x0 = 500^(1/3),
           x1 = 4^(1/3),
           y0 = 0.6 + 0.03,
           y1 = 0.72,
           length = 0.07)
  }
  
  if(counter == 2){
    
    
    text(x = 38000^(1/3),
         y = 0.45,
         pos = 4,
         font = 3,
         labels = "Human\narrival")
    arrows(x0 = 35000^(1/3),
           x1 = 58000^(1/3),
           y0 = 0.45 + 0.03,
           y1 = 0.57,
           length = 0.07)
  }
  

  
  
  dat2 <- dat %>% filter(focal_years == 0) %>% 
    mutate(y1 = n_links_change_current, #log(n_links_change_current),
           y2 = n_nodes_change_current) #log(n_nodes_change_current))
  
  x_pos <- -1.5
  x_seg_len <- 1

  segments(x0 = x_pos, 
           lwd = 3,
           lend = "butt",
           xpd = T,
           y0 = se_fun(dat2$y1, dir = "hi"),
           y1 = se_fun(dat2$y1, dir = "lo"),
           col = rgb(252,141,98, maxColorValue = 255, alpha = 200))
  segments(x0 = x_pos - x_seg_len, x1 = x_pos + x_seg_len, y0 = mean(dat2$y1),
         lwd = 2,
         xpd = T,
         pch = "-",
         col = rgb(252,141,98, maxColorValue = 255))
  
  
  

  segments(x0 = x_pos, 
           lwd = 3,
           lend = "butt",
           xpd = T,
           y0 = se_fun(dat2$y2, dir = "hi"),
           y1 = se_fun(dat2$y2, dir = "lo"),
           col = rgb(141,160,203, maxColorValue = 255, alpha = 200))
  segments(x0 = x_pos - x_seg_len, x1 = x_pos + x_seg_len, y0 = mean(dat2$y2),
         lwd = 2,
         xpd = T,
         pch = "-",
         col = rgb(141,160,203, maxColorValue = 255))
  
  counter <- counter + 1
  
}

dev.off()



# Some data analysis ------------------------

plot.sma.fit <- function(mod){
  curve(x*1, add = T, lty = 2, xpd = F,
        from = min(mod$data, na.rm = T),
        to = max(mod$data, na.rm = T))
  
  curve(coef(mod)[1] + x * coef(mod)[2], 
        add = T, lwd = 1, col = 1,
        from = min(mod$data, na.rm = T),
        to = max(mod$data, na.rm = T),
        xpd = F)
  
  x <- seq(min(mod$data, na.rm = T),
           max(mod$data, na.rm = T), length.out = 100)
  y1 <- mod$groupsummary$Int_lowCI[1] + x * mod$groupsummary$Slope_lowCI[1]
  y2 <- mod$groupsummary$Int_highCI[1] + x * mod$groupsummary$Slope_highCI[1]
  
  xx <- c(x, rev(x))
  yy <- c(y1, rev(y2))
  
  polygon(xx, yy, col = rgb(0,0,0,0.3), border = F, xpd = F)
  
}

plot(n_links_current ~ n_links_pres_nat, 
     data = grid_web_metrics,
     xlab = "Number of links (present natural)",
     ylab = "Number of links (current)",
     asp = 1,
     frame = F,
     pch = 16,
     col = "grey",
     xlim = c(0,1500),
     ylim = c(0,1500))
n_link_mod <- smatr::sma(n_links_current ~ n_links_pres_nat, 
                         data = grid_web_metrics)
plot.sma.fit(n_link_mod)

plot(c(n_links_current/n_nodes_current) ~ c(n_links_pres_nat/n_nodes_pres_nat), 
     data = grid_web_metrics,
     xlab = "Linkage density (present natural)",
     ylab = "Linkage density (current)",
     asp = 1,
     frame = F,
     pch = 16,
     col = "lightgrey",
     xlim = c(0,15),
     ylim = c(0,15))
link_dens_mod <- smatr::sma(c(n_links_current/n_nodes_current) ~ c(n_links_pres_nat/n_nodes_pres_nat), 
                            data = grid_web_metrics)
plot.sma.fit(link_dens_mod)



# Some map versions
cell <- 1:(142*360) %>% rev()
grid_long <- left_join(tibble(cell = cell), grid_web_metrics) %>% 
  mutate(x = (cell %% 360),
         y = 360 - floor(cell/360)) %>% 
  mutate(linkage_density_current = c(n_links_current/n_nodes_current),
         linkage_density_pres_nat = c(n_links_pres_nat/n_nodes_pres_nat)) %>% 
  mutate(delta_linkage_density = linkage_density_current / linkage_density_pres_nat,
         delta_nodes = n_nodes_current / n_nodes_pres_nat,
         delta_links = n_links_current / n_links_pres_nat,) %>% 
  filter(!is.na(n_links_pres_nat))


# # Singleplot
# ggplot(grid_long, aes(x, y)) +
#   geom_point(aes(colour = delta_linkage_density)) +
#   scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 1,
#                          name = "Proportional\nchange in\nlinkage density") +
#   theme_void() +
#   coord_equal()



# Reshape to get three panel
grid_longer <- grid_long %>% 
  pivot_longer(c(delta_linkage_density, delta_nodes, delta_links), names_to = "metric", values_to = "value")

metric_names <- list(
  "delta_linkage_density"="Linkage density",
  "delta_links"="Links",
  "delta_nodes"="Nodes"
)

grid_longer$metric <- factor(grid_longer$metric, levels = c("delta_nodes",
                                                            "delta_links",
                                                            "delta_linkage_density"), 
                             labels = c("Numer of species",
                                        "Number of food web links",
                                        "Linkage density (mean links per species)"))

pdf(file = "change maps.pdf", width = 6.5, height = 6.5)
grid_longer %>% 
  mutate(value = ifelse(value > 2, 2, value)) %>%
  mutate(value = ifelse(value < 0.25, 0.25, value)) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(colour = log(value)), size = 0.1, shape = 15) +
  facet_wrap(vars(metric), nrow=3) +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0,
                         name = "Percent change\n(present / natural)\n",
                         breaks=c(log(2), log(1), log(1- 0.5), log(1 - 0.75)), 
                         labels=c(">200% gain", "0% (no change)", "50% decline", ">75% decline"),
                         limits = c(log(1 - 0.75), log(2))) +
  theme_void() +
  coord_equal()
dev.off()









# Network examples

library("tidyverse")
library("bipartite")
library("igraph")
#install.packages("intergraph")
library("intergraph")
library("ggnetwork")

set.seed(5)

i <- 9802

dt <- as.data.frame(web_current[[i]]) %>% 
  select(1,2) %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix() %>% 
  as.matrix()
n <- as.network(dt, directed = F, bipartite = F)

is_predator <- ifelse(network::get.vertex.attribute(n, attrname = "vertex.names") %in% web_current[[i]]$consumer_sp,
                      "yes", "no")
n <- network::set.vertex.attribute(n, # the name of the network object
                                   "Mammal_predator", # the name we want to reference the variable by in that object
                                   is_predator # the value we are giving that variable
)

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodelabel_repel(aes(label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  geom_nodes(aes(colour = Mammal_predator), size = 8) +
  theme_blank()



dt <- as.data.frame(web_pres_nat[[i]]) %>% 
  select(1,2) %>% 
  graph_from_data_frame() %>% 
  as_adjacency_matrix() %>% 
  as.matrix()
n <- as.network(dt, directed = F, bipartite = F)

is_predator <- ifelse(network::get.vertex.attribute(n, attrname = "vertex.names") %in% web_pres_nat[[i]]$consumer_sp,
                      "yes", "no")
n <- network::set.vertex.attribute(n, # the name of the network object
                                   "Mammal_predator", # the name we want to reference the variable by in that object
                                   is_predator # the value we are giving that variable
)

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodelabel_repel(aes(label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  geom_nodes(aes(colour = Mammal_predator), size = 8) +
  theme_blank()



