# This script pulls diet data from the CarniDIET database
# https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.13296
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.2v6wwpzmr

library('tidyverse')
library("RCurl")

# Pull in the CarniDIET data
carnidiet_tidy <- getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/CarniDIET%201.0.csv") %>% 
  read.csv(text = .) %>% tibble()

# Do a little renaming

carnidiet_tidy <- carnidiet_tidy %>% 
  rename(consumer_sp = scientificNameCarni,
         resource_sp = scientificNamePrey,
         data_source = sourcePrimaryReference) %>% 
  mutate(consumed = 1) %>% 
  mutate(consumer_sp = gsub("_", " ", consumer_sp, fixed = T),
         resource_sp = gsub("_", " ", resource_sp, fixed = T))


