# install.packages("rglobi")
library("rglobi")
library("tidyverse")

mam_checklist <- read.csv("mam_checklist.csv")[,-1] %>% tibble()
pred_tidy <- read.csv("pred_tidy.csv")[,-1] %>% tibble()

mp <- mam_checklist %>% 
  #filter(mammal_predator == 1) %>% 
  select(species) %>% 
  unlist() %>% sort() %>% unique()

source_is_consumer <- c("preysOn", "eats", "kills")
source_is_resource <- c("preyedUponBy", "eatenBy", "killedBy")
intx_types <- c(source_is_consumer, source_is_resource)

globi_dat <- get_interactions(taxon = mp[1], interaction.type = intx_types)

for(i in mp[-1]){ # Not necessary to put in a for loop...
  globi_dat <- globi_dat %>% bind_rows(get_interactions(taxon = i, 
                                                        interaction.type = intx_types))
  print(i)
}

globi_dat$interaction_type %>% table()

# Subset to just the relevant interaction types

globi_set <- globi_dat %>% filter(interaction_type %in% intx_types)

dim(globi_set)

globi_set$consumer_sp <- ifelse(globi_set$interaction_type %in% source_is_consumer,
                                globi_set$source_taxon_name, globi_set$target_taxon_name)

globi_set$resource_sp <- ifelse(globi_set$interaction_type %in% source_is_consumer,
                                globi_set$target_taxon_name, globi_set$source_taxon_name)

head(globi_set)

# Many of these are duplicates of each other
globi_set <- globi_set[!duplicated(globi_set[,c("consumer_sp", "resource_sp")]),]

# Want to add in a consumed = 1

globi_set <- globi_set %>% 
  mutate(consumed = 1)

# Output data
write.csv(file = "globi_tidy.csv", globi_set)

