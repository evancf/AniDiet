# A lot of this code should actually go to the 07 prep section

# install.packages("tidymodels")
# install.packages("keras")
library("keras")
library("tidymodels")
library("tidyverse")

# Load trait data
impute_trait_data <- read.csv("impute_trait_data.csv")[,-1] %>% tibble()

# Read in interaction data
intx_tidy <- read.csv("intx_tidy.csv")[,-1] %>% tibble()

# Want to have order info in here too... Should have this in the 07 code
order_dat <- impute_trait_data %>% 
  select(phylacine_binomial, order)
intx_tidy <- intx_tidy %>% 
  left_join(order_dat, by = c("consumer_sp" = "phylacine_binomial")) %>% 
  left_join(order_dat, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) 

# Couple little things to fix
# Make sure the outcome variable is a factor
intx_tidy$consumed <- factor(intx_tidy$consumed)

# Want to remove non-terrestrial species
intx_tidy <- filter(intx_tidy, foraging_stratum_r != -1)

# Couple columns where there is no variation
cols_to_skip <- which(apply(intx_tidy, 2, function(x) length(levels(factor(x)))) == 1)
intx_tidy <- intx_tidy %>% 
  select(-all_of(cols_to_skip))


# Keep only the columns used for analysis
sp_col <- c("consumer_sp", "resource_sp")
outcome_col <- "consumed"
trait_col <- colnames(intx_tidy)[which(colnames(intx_tidy) == "adult_mass_g_c"): ncol(intx_tidy)]

intx_tidy <- intx_tidy %>% 
  select(all_of(c(sp_col, outcome_col, trait_col, "order_c", "order_r")))

str(intx_tidy) # Looks good - note all numeric even for ones that are binary or ordered categories

# Get "short" version of the dataset with only a single combination of
# each consumer and resource species. In other words, reported as an interaction
# anywhere or not

intx_short <- intx_tidy[order(intx_tidy$consumed, decreasing = T),]
intx_short <- intx_short[!duplicated(select(intx_short, c(consumer_sp, resource_sp))),]
dim(intx_short)
table(intx_short$consumed)

# Lastly, we will add data showing that consumers of one species are not 
# in the diet of the other (unless otherwise known)

# First just reverse the names
intx_inverse <- intx_short %>% 
  select(consumer_sp, resource_sp, consumed) %>% 
  filter(consumed == 1) %>% 
  mutate(consumer_sp1 = resource_sp,
         resource_sp1 = consumer_sp,
         consumed1 = 0) %>% 
  select(-consumer_sp, -resource_sp, -consumed) %>% 
  mutate(consumer_sp = consumer_sp1,
         resource_sp = resource_sp1,
         consumed = consumed1) %>% 
  select(-consumer_sp1, -resource_sp1, -consumed1)

# Next remove any where the opposite is known to be true
intx_short_consumed1 <- intx_short %>% filter(consumed == 1)
observed_consumed <- select(intx_short_consumed1, consumer_sp, resource_sp) %>% apply(1, function(x) paste(x, collapse = " x "))
inverse_combos <- select(intx_inverse, consumer_sp, resource_sp) %>% apply(1, function(x) paste(x, collapse = " x "))
intx_inverse <- intx_inverse %>% 
  filter(!inverse_combos %in% observed_consumed) # Note these are mainly carnivores

intx_inverse <- intx_inverse %>% 
  left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
  left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>% 
  select(all_of(c(sp_col, outcome_col, trait_col, "order_c", "order_r"))) %>% 
  filter(order_c != "Carnivora") # Will remove the carnivores




# Now add these on to the intx_short
intx_short$consumed <- as.numeric(as.character(intx_short$consumed))
intx_short <- intx_short %>% bind_rows(intx_inverse)
intx_short$consumed <- factor(intx_short$consumed)



# Split into training and test

#intx_split <- intx_tidy %>% 
intx_split <- intx_short %>% 
  select(-all_of(sp_col)) %>% 
  initial_split()
train_data <- training(intx_split) 
test_data <- testing(intx_split)



biv_rec <- 
  recipe(consumed ~ ., data = train_data) %>%
  step_BoxCox(all_numeric())%>%
  step_normalize(all_numeric()) %>%
  prep(training = train_data, retain = TRUE)

# We will bake(new_data = NULL) to get the processed training set back

# For validation:
val_normalized <- bake(biv_rec, new_data = test_data, all_predictors())

# For testing when we arrive at a final model: 
test_normalized <- bake(biv_rec, new_data = test_data, all_predictors())

# train_normalized <- bake(biv_rec, new_data = train_data, all_predictors())


set.seed(4)
nnet_fit <-
  mlp(epochs = 50, hidden_units = 5, dropout = 0.1) %>%
  set_mode("classification") %>% 
  # Also set engine-specific `verbose` argument to prevent logging the results: 
  set_engine("keras", verbose = T) %>%
  fit(consumed ~ ., data = bake(biv_rec, new_data = NULL))

nnet_fit


val_results <- 
  test_data %>%
  bind_cols(
    predict(nnet_fit, new_data = val_normalized),
    predict(nnet_fit, new_data = val_normalized, type = "prob")
  )
#val_results

val_results %>% roc_auc(truth = consumed, .pred_0)

val_results %>% accuracy(truth = consumed, .pred_class)

val_results %>% conf_mat(truth = consumed, .pred_class)

caret::confusionMatrix(factor(val_results$consumed), 
                factor(val_results$.pred_class))


# tr_results <- 
#   train_data %>%
#   bind_cols(
#     predict(nnet_fit, new_data = train_normalized),
#     predict(nnet_fit, new_data = train_normalized, type = "prob")
#   )
# tr_results
# 
# tr_results %>% roc_auc(truth = consumed, .pred_0)
# 
# tr_results %>% accuracy(truth = consumed, .pred_class)
# 
# tr_results %>% conf_mat(truth = consumed, .pred_class)


# Make world wide web predictions ------------------------------

# Load matrix of mammal presence
load("m.mamm.pres.nat.RData")
load("m.mamm.current.RData")

rownames(m.mamm.pres.nat) <- gsub("/Users/efricke/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/Phylacine 1.2.1/Ranges/Present natural/",
                                  "", 
                                  rownames(m.mamm.pres.nat), fixed = T)


# Remove those that forage in marine environments
impute_trait_data <- impute_trait_data %>% 
  filter(foraging_stratum != -1)

# Also can skip these marine foraging species in the range matrices
m.mamm.pres.nat <- m.mamm.pres.nat[rownames(m.mamm.pres.nat) %in% impute_trait_data$phylacine_binomial,]
m.mamm.current <- m.mamm.current[rownames(m.mamm.current) %in% impute_trait_data$phylacine_binomial,]


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


xx <- seq(1, 360, by = 4)
yy <- seq(1, 142, by = 4)

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

for(i in cells_to_sample){
  
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
    left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) 
  
  dat_current <- dat_current %>% 
    left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
    left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) 
  
  
  # Make predictions
  dat_pres_nat_normalized <- bake(biv_rec, new_data = dat_pres_nat, all_predictors())
  dat_current_normalized <- bake(biv_rec, new_data = dat_current, all_predictors())
  
  dat_pres_nat <- dat_pres_nat %>%
    bind_cols(predict(nnet_fit, new_data = dat_pres_nat_normalized),
              predict(nnet_fit, new_data = dat_pres_nat_normalized, type = "prob"))
  
  dat_current <- dat_current %>%
    bind_cols(predict(nnet_fit, new_data = dat_current_normalized),
              predict(nnet_fit, new_data = dat_current_normalized, type = "prob"))
  
  dat_current$.pred_1 %>% hist()
  dat_pres_nat$.pred_1 %>% hist()
  
  web_pres_nat[[i]] <- dat_pres_nat %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp")
  web_current[[i]] <- dat_current %>% filter(.pred_class == 1) %>% select("consumer_sp", "resource_sp")
  
  
  print(i)
}






# # https://www.r-bloggers.com/2020/08/how-to-use-lightgbm-with-tidymodels/
# 
# library(doParallel)
# all_cores <- parallel::detectCores(logical = FALSE)
# registerDoParallel(cores = all_cores)
# 
# # set the random seed so we can reproduce any simulated results.
# set.seed(1234)
# # load the housing data and clean names
# ames_data <- make_ames() %>%
#   janitor::clean_names()
# 
# 
# ames_split <- rsample::initial_split(
#   ames_data,
#   prop = 0.8,
#   strata = sale_price
# )
# 
# 
# preprocessing_recipe <-
#   recipes::recipe(sale_price ~ ., data = training(ames_split)) %>%
#   # combine low frequency factor levels
#   recipes::step_other(all_nominal(), threshold = 0.01) %>%
#   # remove no variance predictors which provide no predictive information 
#   recipes::step_nzv(all_nominal()) %>%
#   # prep the recipe so it can be used on other data
#   prep()
# 
# 
# 
# ames_cv_folds <-
#   recipes::bake(
#     preprocessing_recipe,
#     new_data = training(ames_split)
#   ) %>%
#   rsample::vfold_cv(v = 5)
# 
# 
# 
# lightgbm_model<-
#   parsnip::boost_tree(
#     mode = "regression",
#     trees = 1000,
#     min_n = tune(),
#     tree_depth = tune(),
#   ) %>%
#   set_engine("lightgbm", objective = "reg:squarederror",verbose=-1)
# 
# 
# 
# lightgbm_params <-
#   dials::parameters(
#     # The parameters have sane defaults, but if you have some knowledge 
#     # of the process you can set upper and lower limits to these parameters.
#     min_n(), # 2nd important
#     tree_depth() # 3rd most important
#   )