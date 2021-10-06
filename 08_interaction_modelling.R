# install.packages("tidymodels")
# install.packages("tidypredict")
# install.packages("keras")
library("keras")
library("tidymodels")
library("tidyverse")
library("bipartite")


# Read in data -----------------------------------------------------------------
intx_short <- read.csv("data/intx_short.csv")[,-1] %>% tibble()
intx_short$consumed <- factor(intx_short$consumed)


# Split into training and test -------------------------------------------------
set.seed(4)
sp_col <- c("consumer_sp", "resource_sp", "mammal_predator_c", "mammal_predator_r")
intx_split <- intx_short %>% 
  select(-all_of(sp_col)) %>% 
  initial_split(prop = 0.75)
train_data <- training(intx_split) 
test_data <- testing(intx_split)


# Make recipe ------------------------------------------------------------------

nnet_rec <- 
  recipe(consumed ~ ., data = train_data) %>%
  step_YeoJohnson(all_numeric())%>%
  step_normalize(all_numeric()) %>%
  prep(training = train_data, retain = TRUE)

# For looking at model predicted vs observed training data
train_normalized <- bake(nnet_rec, new_data = train_data, all_predictors())
test_normalized <- bake(nnet_rec, new_data = test_data, all_predictors())


# Run model --------------------------------------------------------------------
# Here we are using tuned parameters
set.seed(4)
tensorflow::set_random_seed(4)
nnet_fit <-
  #mlp(epochs = 557, hidden_units = 10, dropout = 0.778, activation = "elu") %>%
  mlp(epochs = 223, hidden_units = 6, dropout = 0.0991, activation = "elu") %>%
  set_mode("classification") %>% 
  set_engine("keras", verbose = T) %>%
  fit(consumed ~ ., data = bake(nnet_rec, new_data = NULL))

nnet_fit


# Model performance  -----------------------------------------------------------

# Multi-trait deep learning model
tr_results <-
  train_data %>%
  bind_cols(
    predict(nnet_fit, new_data = train_normalized),
    predict(nnet_fit, new_data = train_normalized, type = "prob")
  )
tr_results

tr_results %>% roc_auc(truth = consumed, .pred_0)

tr_results %>% accuracy(truth = consumed, .pred_class)

tr_results %>% conf_mat(truth = consumed, .pred_class)

tr_conf_mat <- caret::confusionMatrix(factor(tr_results$consumed), 
                       factor(tr_results$.pred_class))
tr_conf_mat
BIOMOD::TSS.Stat(tr_conf_mat$table)


test_results <-
  test_data %>%
  bind_cols(
    predict(nnet_fit, new_data = test_normalized),
    predict(nnet_fit, new_data = test_normalized, type = "prob")
  )
test_results

test_results %>% roc_auc(truth = consumed, .pred_0)

test_results %>% accuracy(truth = consumed, .pred_class)

test_results %>% conf_mat(truth = consumed, .pred_class)

test_conf_mat <- caret::confusionMatrix(factor(test_results$consumed), 
                                      factor(test_results$.pred_class))
test_conf_mat
BIOMOD::TSS.Stat(test_conf_mat$table)



# Testing performance for a trait matching logistic regression using body mass ratio
train_data$mass_ratio <- train_data$adult_mass_g_c / train_data$adult_mass_g_r
glm_fit <- glm(consumed ~ mass_ratio, data = train_data, family = "binomial")
summary(glm_fit)

test_data$mass_ratio <- test_data$adult_mass_g_c / test_data$adult_mass_g_r
glm_val <- predict(glm_fit, newdata = test_data, type = "response")

library("pROC")
ModelMetrics::auc(actual = test_data$consumed, predicted = glm_val)
set.seed(4)
# Note that I will make estimates probabalistically since no predicted value is greater than 0.5
glm_conf_mat <- caret::confusionMatrix(test_data$consumed, factor(rbinom(length(glm_val), 1, prob = glm_val), levels = c(0,1)))
glm_conf_mat
BIOMOD::TSS.Stat(glm_conf_mat$table)

# Testing performance of genus-level phylogenetic model
# Note that I need to re-do the split so that I can include the consumer and resource names
set.seed(4)
intx_split <- intx_short %>% 
  initial_split(prop = 0.75)
train_data <- training(intx_split) 
test_data <- testing(intx_split)

genus_combos <- train_data %>% 
  filter(consumed == 1) %>% 
  tibble(consumer_gen = word(consumer_sp, 1),
         resource_gen = word(resource_sp, 1)) %>% 
  select(consumer_gen, resource_gen)
genus_combos <- paste(genus_combos$consumer_gen,
                      genus_combos$resource_gen,
                      sep = " x ")

test_combos <- paste(word(test_data$consumer_sp, 1),
                     word(test_data$resource_sp, 1),
                     sep = " x ")

test_data$consumed_predicted <- ifelse(test_combos %in% genus_combos, 1, 0)

ModelMetrics::auc(actual = test_data$consumed, predicted = test_data$consumed_predicted)
genus_conf_mat <- caret::confusionMatrix(test_data$consumed, factor(test_data$consumed_predicted))
genus_conf_mat
BIOMOD::TSS.Stat(genus_conf_mat$table)


# Finalized model ---------------------------------------------------------------

# Will fit the finalized model using all the data here
train_data <- intx_short %>% select(-all_of(sp_col))
train_normalized <- bake(nnet_rec, new_data = train_data, all_predictors())

set.seed(4)
tensorflow::set_random_seed(4)
nnet_fit <-
  #mlp(epochs = 557, hidden_units = 10, dropout = 0.778, activation = "elu") %>%
  mlp(epochs = 223, hidden_units = 6, dropout = 0.0991, activation = "elu") %>%
  set_mode("classification") %>% 
  set_engine("keras", verbose = T) %>%
  fit(consumed ~ ., data = bake(nnet_rec, new_data = NULL))

nnet_fit

tr_results <-
  train_data %>%
  bind_cols(
    predict(nnet_fit, new_data = train_normalized),
    predict(nnet_fit, new_data = train_normalized, type = "prob")
  )
tr_results

tr_results %>% roc_auc(truth = consumed, .pred_0)

tr_results %>% accuracy(truth = consumed, .pred_class)

tr_results %>% conf_mat(truth = consumed, .pred_class)

tr_conf_mat <- caret::confusionMatrix(factor(tr_results$consumed), 
                                      factor(tr_results$.pred_class))
tr_conf_mat
BIOMOD::TSS.Stat(tr_conf_mat$table)



# Make world wide web predictions ----------------------------------------------

# Load trait data (Need to run 07_prep_data_for_interaction_modelling to get this)
impute_trait_data <- read.csv("data/impute_trait_data.csv")[,-1] %>% tibble()


# Diversion to get list of primarily carnivorous species 
# First, pull potential mammals from carnidiet
potential_mammals <- RCurl::getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/Supplementary%20data/Potential%20species%20list.csv") %>% 
  read.csv(text = .) %>% tibble()

# Need to also get this for extinct species
# (NOTE that this code is replicated in '09_web_change_metrics.R')
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

# Get webs at every single site.
i <- 120 # Arctic
i <- 9802 # North America
i <- 11809 # China
i <- 23325 # SE Asia
i <- 20353 # Africa

# Will save the networks in a list object
web_pres_nat <- list()
web_current <- list()
web_no_endangered <- list()

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
      web_pres_nat[[i]] <- tmp
    } else{
      # Make predictions
      dat_pres_nat_normalized <- bake(nnet_rec, new_data = dat_pres_nat, all_predictors())
      dat_pres_nat <- dat_pres_nat %>%
        bind_cols(predict(nnet_fit, new_data = dat_pres_nat_normalized),
                  predict(nnet_fit, new_data = dat_pres_nat_normalized, type = "prob"))
      # Save this to output list
      web_pres_nat[[i]] <- dat_pres_nat %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp") %>% unique()
    }
    
  } else{ # Cases where there's only one species 
    tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
    colnames(tmp) <- c("consumer_sp", "resource_sp")
    web_pres_nat[[i]] <- tmp
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
      web_current[[i]] <- tmp
    } else{
      # Make predictions
      dat_current_normalized <- bake(nnet_rec, new_data = dat_current, all_predictors())
      dat_current <- dat_current %>%
        bind_cols(predict(nnet_fit, new_data = dat_current_normalized),
                  predict(nnet_fit, new_data = dat_current_normalized, type = "prob"))
      # Save this to output list
      web_current[[i]] <- dat_current %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp") %>% unique()
    }
    
  } else{ # Cases where there's only one species 
    tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
    colnames(tmp) <- c("consumer_sp", "resource_sp")
    web_current[[i]] <- tmp
  }
  
  
  
  
  # Get all mammal species in the no_endangered scenario
  spp_no_endangered <- rownames(m.mamm.no.endangered)[m.mamm.no.endangered[,i]]
  if(length(spp_no_endangered) >= 2){
    # Get combination df
    dat_no_endangered <- get_combn_df(spp_no_endangered, sp_combn = F, gen_combn = F)
    # Join with traits
    dat_no_endangered <- dat_no_endangered %>% 
      left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
      left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>% 
      filter(mammal_predator_c)
    if(dim(dat_no_endangered)[1] == 0){ # Cases where we dont have any predators
      tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
      colnames(tmp) <- c("consumer_sp", "resource_sp")
      web_no_endangered[[i]] <- tmp
    } else{
      # Make predictions
      dat_no_endangered_normalized <- bake(nnet_rec, new_data = dat_no_endangered, all_predictors())
      dat_no_endangered <- dat_no_endangered %>%
        bind_cols(predict(nnet_fit, new_data = dat_no_endangered_normalized),
                  predict(nnet_fit, new_data = dat_no_endangered_normalized, type = "prob"))
      # Save this to output list
      web_no_endangered[[i]] <- dat_no_endangered %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp") %>% unique()
    }
    
  } else{ # Cases where there's only one species 
    tmp <- as.data.frame(matrix(nrow=0, ncol = 2))
    colnames(tmp) <- c("consumer_sp", "resource_sp")
    web_no_endangered[[i]] <- tmp
  }
  
  print(ind / length(cells_to_sample))
  
}

# Save this 
save(web_pres_nat, web_current, web_no_endangered, file = "data/hindcast_webs.RData")
