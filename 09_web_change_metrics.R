library("tidyverse")
library("tidymodels")


# Load the fit model ------------------------------

load("intx_tidymodel.RData")


# Load extinction dates from Andermann ------------------------------

te_dates <- read.table("Andermann/global_pyrate_species_list.txt", header = T) %>% tibble() %>% 
  bind_cols(read.table("Andermann/global_ts_te_dates.txt", header = T) %>% tibble())

te_dates <- te_dates %>% 
  mutate(binomial = paste(id, taxon, sep = " ")) %>% 
  select(-starts_with("ts"))

# Decide on some relevant focal years (rather than every year)

focal_years <- (seq(0, (126000)^(1/3), length.out = 8))^3


# Make world wide web predictions ------------------------------

# Load trait data
impute_trait_data <- read.csv("impute_trait_data.csv")[,-1] %>% tibble()

# Load matrix of mammal presence
load("m.mamm.pres.nat.RData")
load("m.mamm.current.RData")


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

# Lastly, remove humans
m.mamm.pres.nat <- m.mamm.pres.nat[rownames(m.mamm.pres.nat) != "Homo sapiens",]
m.mamm.current <- m.mamm.current[rownames(m.mamm.current) != "Homo sapiens",]

# Function to get all pairwise consumer / resource combination given a species list
get_combn_df <- function(x, sp_combn = T, gen_combn = T){
  
  dat <- cbind(combn(x, 2), combn(x, 2)[2:1,]) %>% t() %>% as.data.frame()
  colnames(dat) <- c("consumer_sp", "resource_sp")
  
  if(sp_combn){
    dat$sp_combn <- paste(dat$consumer_sp, dat$resource_sp)
  }
  
  if(gen_combn){
    dat$gen_combn <- paste(word(dat$consumer_sp, 1), word(dat$resource_sp, 1))
  }
  
  dat
}


xx <- seq(1, 360, by = 3)
yy <- seq(1, 142, by = 3)


# Get some cells to actually sample (not every single one at this point)
cells_to_sample <- apply(expand.grid(yy, xx),
                         1,
                         function(x) ((x[1] - 1) * 360 + x[2])) %>% sort()

# Also want to skip ones in the ocean
cells_to_sample <- cells_to_sample[cells_to_sample %in% which(colSums(m.mamm.pres.nat)>1)]

# Get webs at every single site
i <- 120 # Arctic
i <- 9802 # North America
i <- 11809 # China
i <- 23325 # SE Asia
i <- 20353 # Africa
web_pres_nat <- list()
web_current <- list()


# For loop to make webs at every study location

for(i in cells_to_sample){
  ind <- which(cells_to_sample == i)
  
  spp_pres_nat <- rownames(m.mamm.pres.nat)[m.mamm.pres.nat[,i]]
  spp_current <- rownames(m.mamm.current)[m.mamm.current[,i]]
  
  # Will skip ones where there's only one species in the pixel
  if(length(spp_pres_nat) < 2 | length(spp_current) < 2) next()
  
  # Get combination df
  dat_pres_nat <- get_combn_df(spp_pres_nat, sp_combn = F, gen_combn = F)
  dat_current <- get_combn_df(spp_current, sp_combn = F, gen_combn = F)
  
  # Join with traits
  dat_pres_nat <- dat_pres_nat %>% 
    left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
    left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>% 
    filter(det_vend_c != 0 | det_vect_c != 0 | det_scav_c != 0 | det_vunk_c != 0)
  
  dat_current <- dat_current %>% 
    left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
    left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>% 
    filter(det_vend_c != 0 | det_vect_c != 0 | det_scav_c != 0 | det_vunk_c != 0)
  
  # Will skip ones where there's only one species in the pixel
  if(dim(dat_pres_nat)[1] < 2 | dim(dat_current)[1] < 2) next()
  
  # Make predictions
  dat_pres_nat_normalized <- bake(biv_rec, new_data = dat_pres_nat, all_predictors())
  dat_current_normalized <- bake(biv_rec, new_data = dat_current, all_predictors())
  
  dat_pres_nat <- dat_pres_nat %>%
    bind_cols(predict(nnet_fit, new_data = dat_pres_nat_normalized),
              predict(nnet_fit, new_data = dat_pres_nat_normalized, type = "prob"))
  
  dat_current <- dat_current %>%
    bind_cols(predict(nnet_fit, new_data = dat_current_normalized),
              predict(nnet_fit, new_data = dat_current_normalized, type = "prob"))
  
  if(sum(as.numeric(as.character(dat_pres_nat$.pred_class))) < 2 |
     sum(as.numeric(as.character(dat_current$.pred_class))) < 2) next()
  
  # dat_current$.pred_1 %>% hist()
  # dat_pres_nat$.pred_1 %>% hist()
  
  # dat_pres_nat %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp")
  # dat_pres_nat %>% filter(.pred_1 > 0.1) %>% select("consumer_sp", "resource_sp")
  
  web_pres_nat[[i]] <- dat_pres_nat %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp") %>% unique()
  web_current[[i]] <- dat_current %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp") %>% unique()
  
  # # Turn these into a cheddar community
  # ched_pres_nat <- Community(nodes = data.frame(node = unique(unlist(web_pres_nat[[i]]))),
  #                            properties = list(title = 10137),
  #                            trophic.links = web_pres_nat[[i]] %>% rename(consumer = consumer_sp,
  #                                                          resource = resource_sp))
  # ched_current <- Community(nodes = data.frame(node = unique(unlist(web_current[[i]]))),
  #                            properties = list(title = 10137),
  #                            trophic.links = web_current[[i]] %>% rename(consumer = consumer_sp,
  #                                                                         resource = resource_sp))
  # 
  # 
  # 
  # grid_web_metrics[ind,-1] <- c(ched_pres_nat %>% NumberOfTrophicLinks(), # n_links_pres_nat
  #                               ched_pres_nat %>% NumberOfNodes(), # n_nodes_pres_nat
  #                               #TrophicChainsStats(ched_pres_nat)$chain.lengths %>% mean(), # mean_chain_length_pres_nat
  #                               ched_current %>% NumberOfTrophicLinks(), # n_links_current
  #                               ched_current %>% NumberOfNodes()#, # n_nodes_current
  #                               #TrophicChainsStats(ched_current)$chain.lengths %>% mean()
  #                               ) # mean_chain_length_current
  
  print(i / 50000)
  
}

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

range(grid_long$red, na.rm=T)
range(grid_long$blue, na.rm=T)

grid_long

# plot(y ~ x,
#      pch = ifelse(gain_loss == "gain", 
#                   17, 16),
#      col = rgb(red, 0, blue),
#      cex = log10(n_nodes_pres_nat),
#   data = grid_long,
#   asp = 1)



ggplot(grid_long, aes(x, y)) +
  geom_point(aes(colour = delta_linkage_density)) +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 1,
                         name = "Proportional\nchange in\nlinkage density") +
  theme_void() +
  coord_equal()



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
grid_longer %>% 
  mutate(value = ifelse(value > 2, 2, value)) %>%
  mutate(value = ifelse(value < 0.25, 0.25, value)) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(colour = log(value))) +
  facet_wrap(vars(metric), nrow=3) +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0,
                         name = "Percent change\n(present / natural)\n",
                         breaks=c(log(2), log(1), log(1- 0.5), log(1 - 0.75)), 
                         labels=c(">200% gain", "0% (no change)", "50% decline", ">75% decline"),
                         limits = c(log(1 - 0.75), log(2))) +
  theme_void() +
  coord_equal()


