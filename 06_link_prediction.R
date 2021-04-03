# There is a weird issue with phylopars. Probably want to get PCA first

library("tidyverse")

# require(devtools)
# install_github("ericgoolsby/Rphylopars",dependencies = TRUE)

# Load interaction data -------------------------

# Get a metaweb (including a long version)
pred_tidy <- read.csv("pred_tidy.csv")[,-1] %>% tibble()
globi_tidy <- read.csv("globi_tidy.csv")[,-1] %>% tibble()
source("02c_pull_CarniDiet_data.R")

# This is a little clunky but will join these together
pred_tidy_cols <- colnames(pred_tidy)

pred_tidy <- bind_rows(pred_tidy, globi_tidy, carnidiet_tidy) %>% 
  select(all_of(pred_tidy_cols))


# Make sure we are only using binomial (from pred_tidy)
pred_tidy$consumer_sp <- pred_tidy$consumer_sp %>% word(start = 1, end = 2)
pred_tidy$resource_sp <- pred_tidy$resource_sp %>% word(start = 1, end = 2)

# Clean a couple things up (THIS SHOULD REALLY BE ELSEWHERE!!!!!)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
pred_tidy$consumer_sp <- pred_tidy$consumer_sp %>% tolower() %>% firstup()
pred_tidy$resource_sp <- pred_tidy$resource_sp %>% tolower() %>%firstup()


# Load trait data -------------------------

# Get the COMBINE dataset for traits https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.3344
temp <- tempfile()
download.file("https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3344&file=ecy3344-sup-0001-DataS1.zip",temp)
trait_data <- read.table(unz(temp, "COMBINE_archives/trait_data_imputed.csv"), sep = ",", header = T) %>% tibble()
unlink(temp)

# Going to treat foraging stratum as a numeric
trait_data <- trait_data %>% 
  mutate(foraging_stratum = as.numeric(plyr::revalue(foraging_stratum, c("M" = -1,
                                                        "G" = 0,
                                                        "S" = 1,
                                                        "Ar" = 2,
                                                        "A" = 3))))

chr_columns <- c("order", "family", "genus", "species", 
                  "iucn2020_binomial",
                  "phylacine_binomial")
focal_traits <- c("adult_mass_g", "brain_mass_g", "adult_body_length_mm", "max_longevity_d",
                  "female_maturity_d", "age_first_reproduction_d", "gestation_length_d",
                  "litter_size_n", "litters_per_year_n", "interbirth_interval_d",
                  "weaning_age_d", "generation_length_d", "dispersal_km", "hibernation_torpor",
                  "fossoriality", "dphy_invertebrate", "dphy_vertebrate", "dphy_plant",
                  "det_inv", "det_vend", "det_vect", "det_vfish", "det_vunk", "det_scav",
                  "det_fruit", "det_nect", "det_seed", "det_plantother", "det_diet_breadth_n",
                  "trophic_level", "activity_cycle",
                  "foraging_stratum",
                  #"freshwater", "marine",
                  #"island_endemicity", "biogeographical_realm",
                  "terrestrial_non.volant", "terrestrial_volant", "habitat_breadth_n")

trait_data <- trait_data %>% 
  select(all_of(c(chr_columns,focal_traits)))


# Phylogenetic imputation -------------------------
library("Rphylopars")

temp <- tempfile()
download.file("https://github.com/MegaPast2Future/PHYLACINE_1.2/raw/master/Data/Phylogenies/Complete_phylogeny.nex",temp)
mamm_phylo <- read.nexus(temp)
unlink(temp)

# # Get a consensus tree. Probably want to do this in the future!!!!!
# consensus_mamm_phylo <- phytools::averageTree(mamm_phylo)

# For now, will just use one
mamm_phylo <- mamm_phylo[[1]]

# Get a dataset that we'll actually give to phylopars
mean.na.rm <- function(x) mean(x, na.rm=T)
impute_dat <- trait_data %>% 
  select(all_of(c("phylacine_binomial", focal_traits))) %>% 
  rename(species = phylacine_binomial) %>% 
  mutate(species = gsub(" ", "_", species, fixed = T)) %>% 
  group_by(species) %>% 
  summarize_all(mean.na.rm) %>% 
  filter(species %in% mamm_phylo$tip.label) %>% 
  as.data.frame()

impute_dat[impute_dat == "NaN"] <- NA

#impute_dat[,-1] <- apply(impute_dat[,-1], 2, as.double)

# mini.add <- function(x) x + 0.00001#runif(1, min = 0, max = 0.00000001)
# impute_dat <- impute_dat %>%
#   mutate(across(where(is.numeric), ~ sum(runif(1, min = 0, max = 0.0001), .x)))
# impute_dat


# # Working around a weird issue
# temp_impute_dat <- impute_dat
# temp_mamm_phylo <- mamm_phylo
# temp_mamm_phylo$tip.label <- 1:length(mamm_phylo$tip.label)
# temp_impute_dat$species <- sapply(X = temp_impute_dat$species, FUN = function(x){which(mamm_phylo$tip.label==x)})


bad_cols <- which(apply(impute_dat, 2, function(x) length(levels(as.factor(x)))) < 3)
for(i in bad_cols){
  impute_dat[,i] <- impute_dat[,i] + rnorm(nrow(impute_dat), sd = 0.001)
}

time0 <- Sys.time()
set.seed(4)
phylopars_dat <- phylopars(trait_data = impute_dat, 
                           tree = mamm_phylo,
                           pheno_error = F,
                           pheno_correlated = F)
time1 <- Sys.time()
time1-time0


# Join pred_tidy with trait_data by phylacine_binomial -------------------------
pred_tidy <- pred_tidy %>% 
  left_join(trait_data, by = c("consumer_sp" = "phylacine_binomial")) %>% 
  left_join(trait_data, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r"))


# A little code to figure out what traits to use
sum.is.na <- function(x) sum(is.na(x))
trait_coverage <- 1 - (trait_data %>%
  summarise(across(where(is.numeric), ~ sum.is.na(.x))) / nrow(trait_data)) %>%
  unlist()
# 
# trait_coverage %>% hist(breaks = 20)
# 
# focal_traits <- names(trait_coverage)[(trait_coverage > 0.95)]
# focal_traits <- focal_traits[!focal_traits == "consumed"]






# 


dim(pred_tidy[complete.cases(pred_tidy[,]), ]









