# install.packages("tidymodels")
# install.packages("keras")
library("keras")
library("tidymodels")
library("tidyverse")

# Read in interaction data
intx_tidy <- read.csv("intx_tidy.csv")[,-1] %>% tibble()

# Couple little things to fix
# Make sure the outcome variable is a factor
intx_tidy$consumed <- factor(intx_tidy$consumed)

# Want to remove non-terrestrial species
filter(intx_tidy, foraging_stratum_r != -1)

# Couple columns where there is no variation
cols_to_skip <- which(apply(intx_tidy, 2, function(x) length(levels(factor(x)))) == 1)
intx_tidy <- intx_tidy %>% 
  select(-all_of(cols_to_skip))


# Keep only the columns used for analysis
sp_col <- c("consumer_sp", "resource_sp")
outcome_col <- "consumed"
trait_col <- colnames(intx_tidy)[which(colnames(intx_tidy) == "adult_mass_g_c"): ncol(intx_tidy)]

intx_tidy <- intx_tidy %>% 
  select(all_of(c(sp_col, outcome_col, trait_col)))

str(intx_tidy) # Looks good - note all numeric even for ones that are binary or ordered categories

# Get "short" version of the dataset with only a single combination of
# each consumer and resource species. In other words, reported as an interaction
# anywhere or not

intx_short <- intx_tidy[order(intx_tidy$consumed, decreasing = T),]
intx_short <- intx_short[!duplicated(select(intx_short, c(consumer_sp, resource_sp))),]
dim(intx_short)
table(intx_short$consumed)

# Split into training and test

#intx_split <- intx_tidy %>% 
intx_split <- intx_short %>% 
  select(-all_of(sp_col)) %>% 
  initial_split()
train_data <- training(intx_split) 
test_data <- testing(intx_split)



biv_rec <- 
  recipe(consumed ~ ., data = train_data) %>%
  step_BoxCox(all_predictors())%>%
  step_normalize(all_predictors()) %>%
  prep(training = train_data, retain = TRUE)

# We will bake(new_data = NULL) to get the processed training set back

# For validation:
val_normalized <- bake(biv_rec, new_data = test_data, all_predictors())

# For testing when we arrive at a final model: 
test_normalized <- bake(biv_rec, new_data = test_data, all_predictors())


set.seed(57974)
nnet_fit <-
  mlp(epochs = 300, hidden_units = 5, dropout = 0.1) %>%
  set_mode("classification") %>% 
  # Also set engine-specific `verbose` argument to prevent logging the results: 
  set_engine("keras", verbose = 0) %>%
  fit(consumed ~ ., data = bake(biv_rec, new_data = NULL))

nnet_fit


val_results <- 
  test_data %>%
  bind_cols(
    predict(nnet_fit, new_data = val_normalized),
    predict(nnet_fit, new_data = val_normalized, type = "prob")
  )
val_results

val_results %>% roc_auc(truth = consumed, .pred_0)

val_results %>% accuracy(truth = consumed, .pred_class)

val_results %>% conf_mat(truth = consumed, .pred_class)


# Make world wide web predictions ------------------------------

load("m.mamm.pres.nat.RData")
load("m.mamm.current.RData")

rownames(m.mamm.pres.nat) <- gsub("/Users/efricke/Dropbox/*Science/*Research/*SESYNC/1 Predicting interactions/Data/distribution/range maps/Phylacine 1.2.1/Ranges/Present natural/",
                                  "", 
                                  rownames(m.mamm.pres.nat), fixed = T)

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

# Get big tibble of cooccurrences
cr_pres_nat <- list()
for(i in 1:ncol(m.mamm.pres.nat)){
  spp <- rownames(m.mamm.pres.nat)[m.mamm.pres.nat[,i]]
  
  # Will skip ones where there's only one species in the pixel
  if(length(spp) < 2) next()
  
  cr_pres_nat[[i]] <- get_combn_df(spp, sp_combn = F, gen_combn = F)
  
  if(i %in% seq(1, ncol(m.mamm.pres.nat), by = 2000)){
    print(i)
  }
}

cr_pres_nat <- do.call(rbind, cr_pres_nat)

# Now join on trait data
impute_trait_data <- read.csv("impute_trait_data.csv")[,-1] %>% tibble()

cr_pres_nat <- cr_pres_nat %>% 
  left_join(impute_trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
  left_join(impute_trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r"))

all_intx <- tibble(consumer)



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