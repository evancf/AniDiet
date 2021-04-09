# install.packages("tidymodels")
# install.packages("keras")
library("keras")
library("tidymodels")
library("tidyverse")
library("bipartite")


# Read in data -----------------------------------------------------------------
intx_short <- read.csv("intx_short.csv")[,-1] %>% tibble()
intx_short$consumed <- factor(intx_short$consumed)


# Split into training and test -------------------------------------------------
sp_col <- c("consumer_sp", "resource_sp")
intx_split <- intx_short %>% 
  select(-all_of(sp_col)) %>% 
  initial_split(prop = 0.99)
train_data <- training(intx_split) 
test_data <- testing(intx_split)


# Make recipe ------------------------------------------------------------------

biv_rec <- 
  recipe(consumed ~ ., data = train_data) %>%
  step_BoxCox(all_numeric())%>%
  step_normalize(all_numeric()) %>%
  prep(training = train_data, retain = TRUE)


# For validation:
val_normalized <- bake(biv_rec, new_data = test_data, all_predictors())

# For testing when we arrive at a final model: 
test_normalized <- bake(biv_rec, new_data = test_data, all_predictors())

# For looking at model predicted vs observed training data
train_normalized <- bake(biv_rec, new_data = train_data, all_predictors())


# Run model --------------------------------------------------------------------
set.seed(4)
nnet_fit <-
  mlp(epochs = 50, hidden_units = 5, dropout = 0.1) %>%
  set_mode("classification") %>% 
  # Also set engine-specific `verbose` argument to prevent logging the results: 
  set_engine("keras", verbose = T) %>%
  fit(consumed ~ ., data = bake(biv_rec, new_data = NULL))

nnet_fit


# Model validation  ------------------------------------------------------------
val_results <- 
  test_data %>%
  bind_cols(
    predict(nnet_fit, new_data = val_normalized),
    predict(nnet_fit, new_data = val_normalized, type = "prob"))
#val_results

val_results %>% roc_auc(truth = consumed, .pred_0)

val_results %>% accuracy(truth = consumed, .pred_class)

val_results %>% conf_mat(truth = consumed, .pred_class)

caret::confusionMatrix(factor(val_results$consumed), 
                factor(val_results$.pred_class))


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

caret::confusionMatrix(factor(tr_results$consumed), 
                       factor(tr_results$.pred_class))


# Output tidymodel -------------------------------------------------------------

save(biv_rec, nnet_fit, file = "intx_tidymodel.RData")


# Maybe want to use a lightgbm version?



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