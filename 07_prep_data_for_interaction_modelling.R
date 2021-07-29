library("tidyverse")

# Load trait data -------------------------

# Here is the code to get a trait dataframe with no missing values. It takes
# about 20 minutes or so to run, so we will just pull in from a saved csv if possible,
# although you'll have to run through this first on your computer

if(!"impute_trait_data.csv" %in% list.files()){
  # Get the COMBINE dataset for traits https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.3344
  temp <- tempfile()
  download.file("https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3344&file=ecy3344-sup-0001-DataS1.zip",temp)
  trait_data <- read.table(unz(temp, "COMBINE_archives/trait_data_imputed.csv"), sep = ",", header = T) %>% tibble()
  unlink(temp)
  
  # Going to treat foraging stratum as a numeric (this applies )
  trait_data <- trait_data %>%
    mutate(foraging_stratum = as.numeric(plyr::revalue(foraging_stratum, c("M" = -1,
                                                                           "G" = 0,
                                                                           "S" = 1,
                                                                           "Ar" = 2,
                                                                           "A" = 3))))
  
  # # A little code to figure out what traits to use
  # sum.is.na <- function(x) sum(is.na(x))
  # trait_coverage <- 1 - (trait_data %>%
  #   summarise(across(where(is.numeric), ~ sum.is.na(.x))) / nrow(trait_data)) %>%
  #   unlist()
  #
  # trait_coverage %>% hist(breaks = 20)
  #
  # focal_traits <- names(trait_coverage)[(trait_coverage > 0.95)]
  # focal_traits <- focal_traits[!focal_traits == "consumed"]
  
  
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
  
  
  
  
  # Missing value imputation -------------------------
  
  # A small portion of values
  
  # missForest imputation (could do the phylogenetic approach below, but there's
  # a weird error maybe because of too many variables?)
  # install.packages("missForest")
  library("missForest")
  
  # For this version, we are able to use categories.
  category_columns <- c("order", "family", "genus", 
                        "hibernation_torpor", "fossoriality", "trophic_level",
                        "activity_cycle", "foraging_stratum", "terrestrial_non.volant",
                        "terrestrial_volant", "habitat_breadth_n")
  trait_data <- trait_data %>%
    mutate_at(vars(category_columns),
              list(factor)) %>%
    as.data.frame()
  
  impute_trait_data <- missForest(xmis = trait_data[,c("order",
                                                       focal_traits)]) # Takes about 20 minutes?
  impute_trait_data <- impute_trait_data$ximp %>% 
    bind_cols(trait_data[,chr_columns[-1]])
  
  # Going to save a personal copy of this
  write.csv(file = "impute_trait_data.csv", impute_trait_data)
  
  # # Phylogenetic imputation -------------------------
  # library("Rphylopars")
  # 
  # temp <- tempfile()
  # download.file("https://github.com/MegaPast2Future/PHYLACINE_1.2/raw/master/Data/Phylogenies/Complete_phylogeny.nex",temp)
  # mamm_phylo <- read.nexus(temp)
  # unlink(temp)
  # 
  # # # Get a consensus tree. Probably want to do this in the future!!!!!
  # # consensus_mamm_phylo <- phytools::averageTree(mamm_phylo)
  # 
  # # For now, will just use one
  # mamm_phylo <- mamm_phylo[[1]]
  # 
  # # Get a dataset that we'll actually give to phylopars
  # mean.na.rm <- function(x) mean(x, na.rm=T)
  # impute_dat <- trait_data %>% 
  #   select(all_of(c("phylacine_binomial", focal_traits))) %>% 
  #   rename(species = phylacine_binomial) %>% 
  #   mutate(species = gsub(" ", "_", species, fixed = T)) %>% 
  #   group_by(species) %>% 
  #   summarize_all(mean.na.rm) %>% 
  #   filter(species %in% mamm_phylo$tip.label) %>% 
  #   as.data.frame()
  # 
  # impute_dat[impute_dat == "NaN"] <- NA
  # 
  # #impute_dat[,-1] <- apply(impute_dat[,-1], 2, as.double)
  # 
  # # mini.add <- function(x) x + 0.00001#runif(1, min = 0, max = 0.00000001)
  # # impute_dat <- impute_dat %>%
  # #   mutate(across(where(is.numeric), ~ sum(runif(1, min = 0, max = 0.0001), .x)))
  # # impute_dat
  # 
  # 
  # # # Working around a weird issue
  # # temp_impute_dat <- impute_dat
  # # temp_mamm_phylo <- mamm_phylo
  # # temp_mamm_phylo$tip.label <- 1:length(mamm_phylo$tip.label)
  # # temp_impute_dat$species <- sapply(X = temp_impute_dat$species, FUN = function(x){which(mamm_phylo$tip.label==x)})
  # 
  # 
  # bad_cols <- which(apply(impute_dat, 2, function(x) length(levels(as.factor(x)))) < 3)
  # for(i in bad_cols){
  #   impute_dat[,i] <- impute_dat[,i] + rnorm(nrow(impute_dat), sd = 0.001)
  # }
  # 
  # time0 <- Sys.time()
  # set.seed(4)
  # phylopars_dat <- phylopars(trait_data = impute_dat, 
  #                            tree = mamm_phylo,
  #                            pheno_error = F,
  #                            pheno_correlated = F)
  # time1 <- Sys.time()
  # time1-time0
  
  
}

impute_trait_data <- read.csv("impute_trait_data.csv")[,-1] %>% tibble()

# Diversion to get list of primarily carnivorous species 
potential_mammals <- RCurl::getURL("https://raw.githubusercontent.com/osmiddleton/CarniDIET-Database/master/Version%201.0/Supplementary%20data/Potential%20species%20list.csv") %>% 
  read.csv(text = .) %>% tibble()

mammal_predator <- gsub("_"," ",potential_mammals$Bin.)

impute_trait_data$mammal_predator <- impute_trait_data$phylacine_binomial %in% mammal_predator


# Load interaction data -------------------------

#
pred_tidy <- read.csv("pred_tidy.csv")[,-1] %>% tibble()
globi_tidy <- read.csv("globi_tidy.csv")[,-1] %>% tibble()
source("02c_pull_CarniDiet_data.R")

# For the carniDIET data, we know exact locations. Based on the IUCN
# range maps, we will see which animal species are there

load("mamm.presence.by.carnicoords.wide.RData")

mamm.presence.by.carnicoords.wide <- mamm.presence.by.carnicoords.wide %>% 
  filter(!is.na(binomial)) %>% 
  as.data.frame()

# Rename to easier "mamm.pres" and make binomials the rownames
rownames(mamm.presence.by.carnicoords.wide) <- mamm.presence.by.carnicoords.wide$binomial
mamm.presence.by.carnicoords.wide <- mamm.presence.by.carnicoords.wide %>% select(-binomial)

# Deal with subspecies - want the analyses to be at the species level
mamm.pres.bin <- unique(word(rownames(mamm.presence.by.carnicoords.wide), 1, 2))
mamm.pres <- matrix(NA, nrow = length(mamm.pres.bin),
                    ncol = ncol(mamm.presence.by.carnicoords.wide))
colnames(mamm.pres) <- colnames(mamm.presence.by.carnicoords.wide)
rownames(mamm.pres) <- mamm.pres.bin

# There's got to be a better way to do this. Probably with dplyr::summarize
# after grouping by species...
for(i in mamm.pres.bin){
  inds <- which(word(rownames(mamm.presence.by.carnicoords.wide), 1, 2) == i)
  mamm.pres[i,] <- apply(mamm.presence.by.carnicoords.wide[inds,], 2, any)
  print(paste(which(rownames(mamm.pres) == i), i))
}


carnidiet_tidy <- carnidiet_tidy %>% 
  mutate(carni.id = paste(decimalLatitude,
                          decimalLongitude))

carnidiet_list <- list()

for(i in unique(carnidiet_tidy$carni.id)){
  
  dat <- carnidiet_tidy %>% 
    filter(carni.id == i)
  c_spp <- dat$consumer_sp %>% unique()
  
  zero_dat <- tibble(consumer_sp = c(),
                     consumed = c(),
                     resource_sp = c())
  
  # First, check if there's presence / absense data for this site
  if(i %in% colnames(mamm.pres)){
    
    for(j in c_spp){
      
      r_spp <-  dat %>% filter(consumer_sp == j) %>% 
        select(genusPrey, resource_sp)
      
      # Get all the species reported at this site (according to IUCN range maps)
      zero_spp <- rownames(mamm.pres)[mamm.pres[,i]]
      
      # Remove species observed as resources from zero_spp
      zero_spp <- zero_spp[!zero_spp %in% r_spp$resource_sp]
      
      # Remove species potentially in diet but where only genus-level records exist
      r_gen_only <- r_spp$genusPrey[is.na(r_spp$resource_sp)]
      zero_spp <- zero_spp[!word(zero_spp, 1) %in% r_gen_only]
      
      zero_dat <- zero_dat %>% 
        bind_rows(tibble(consumer_sp = j,
                         consumed = 0,
                         resource_sp = zero_spp))
      
    }
    
    carnidiet_list[[i]] <- bind_rows(dat, zero_dat)
    
  } else{ # If there is not, will only remember the consumed = 1 values
    carnidiet_list[[i]] <- dat
  }
  
}

carnidiet_tidy_with_zeros <- do.call(rbind, carnidiet_list)




# This is a little clunky but will join these together
pred_tidy_cols <- colnames(pred_tidy)

intx_data <- bind_rows(pred_tidy, globi_tidy, carnidiet_tidy_with_zeros) %>%  
  select(all_of(pred_tidy_cols))


# Make sure we are only using binomial (from intx_data)
intx_data$consumer_sp <- intx_data$consumer_sp %>% word(start = 1, end = 2)
intx_data$resource_sp <- intx_data$resource_sp %>% word(start = 1, end = 2)

# Clean a couple things up (THIS SHOULD REALLY BE ELSEWHERE!!!!!)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
intx_data$consumer_sp <- intx_data$consumer_sp %>% tolower() %>% firstup()
intx_data$resource_sp <- intx_data$resource_sp %>% tolower() %>% firstup()




# Join intx_data with trait_data by phylacine_binomial -------------------------
impute_trait_data_phyla <- impute_trait_data[!duplicated(impute_trait_data$phylacine_binomial),] %>% 
  select(-iucn2020_binomial)
impute_trait_data_iucn <- impute_trait_data[!duplicated(impute_trait_data$iucn2020_binomial),] %>% 
  select(-phylacine_binomial)

intx_data_phyla <- intx_data %>% 
  left_join(impute_trait_data_phyla, by = c("consumer_sp" = "phylacine_binomial")) %>% 
  left_join(impute_trait_data_phyla, by = c("resource_sp" = "phylacine_binomial"), suffix = c("_c", "_r")) %>% 
  mutate(zz = 1:nrow(.))

intx_data_iucn <- intx_data %>% 
  left_join(impute_trait_data_iucn, by = c("consumer_sp" = "iucn2020_binomial")) %>% 
  left_join(impute_trait_data_iucn, by = c("resource_sp" = "iucn2020_binomial"), suffix = c("_c", "_r")) %>% 
  mutate(zz = 1:nrow(.))

intx_data_both <- bind_rows(intx_data_phyla, intx_data_iucn) %>% 
  unique() %>% 
  select(-zz, 
         -order_c, -order_r,
         -family_c, -family_r,
         -genus_c, -genus_r,
         -species_c, -species_r)


# 

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

focal_traits2 <- c(paste0(focal_traits, "_r"),
                   paste0(focal_traits, "_c"))

intx_data_both %>% select(all_of(focal_traits2)) %>% complete.cases() %>% table()

# THIS REMOVES A LOT! Especially birds and other diet items!!!!!
intx_tidy <- filter(intx_data_both, complete.cases(select(intx_data_both, all_of(focal_traits2))))



# Couple little things to fix
# Make sure the outcome variable is a factor
intx_tidy$consumed <- factor(intx_tidy$consumed)

# Want to remove non-terrestrial species
intx_tidy <- filter(intx_tidy, foraging_stratum_r != -1)

# Couple columns where there is no variation
cols_to_skip <- which(apply(intx_tidy, 2, function(x) length(levels(factor(x)))) == 1)
intx_tidy <- intx_tidy %>% 
  select(-all_of(cols_to_skip[!names(cols_to_skip) == "mammal_predator_c"]))


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
  select(all_of(c(sp_col, outcome_col, trait_col))) %>% 
  filter(det_vend_c != 0 | det_vect_c != 0 | det_scav_c != 0 | det_vunk_c != 0) # Will remove "consumer" specie with no evidence of being vertebrate carnivores
  #filter(mammal_predator_c)


# Now add these on to the intx_short
intx_short$consumed <- as.numeric(as.character(intx_short$consumed))
intx_short <- intx_short %>% bind_rows(intx_inverse)
intx_short$consumed <- factor(intx_short$consumed)



# Write this to csv
write.csv(file = "intx_short.csv", intx_short)




