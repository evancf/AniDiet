library("tidyverse")

#install.packages("wosr")
library("wosr")


# Get community data ------------------------------------
# This is from Rowan et al. 2019 
# https://www.pnas.org/content/117/3/1559/tab-figures-data

# Will download to this temp file, which gets removed with unlink below.
tmp <- tempfile(fileext = ".xlsx")
sheet <- "Trait Data"
download.file(url = "https://www.pnas.org/highwire/filestream/902133/field_highwire_adjunct_files/2/pnas.1910489116.sd01.xlsx", 
              destfile = tmp, mode="wb")
comm_afr <- readxl::read_excel(path = tmp, sheet = sheet)

download.file(url = "https://www.pnas.org/highwire/filestream/902133/field_highwire_adjunct_files/2/pnas.1910489116.sd02.xlsx", 
              destfile = tmp, mode="wb")
comm_ind <- readxl::read_excel(path = tmp, sheet = sheet)

download.file(url = "https://www.pnas.org/highwire/filestream/902133/field_highwire_adjunct_files/3/pnas.1910489116.sd03.xlsx", 
              destfile = tmp, mode="wb")
comm_mad <- readxl::read_excel(path = tmp, sheet = sheet)

download.file(url = "https://www.pnas.org/highwire/filestream/902133/field_highwire_adjunct_files/4/pnas.1910489116.sd04.xlsx", 
              destfile = tmp, mode="wb")
comm_neo <- readxl::read_excel(path = tmp, sheet = sheet)

unlink(tmp) # Remove the temp file



# Get some species name vectors

# Just the species considered carnivores of mammals
mampred_sp_afr <- filter(comm_afr, Mammal > 0)$Species %>% gsub("_", " ", .)
mampred_sp_ind <- filter(comm_ind, Mammal > 0)$Species %>% gsub("_", " ", .)
mampred_sp_mad <- filter(comm_mad, Mammal > 0)$Species %>% gsub("_", " ", .)
mampred_sp_neo <- filter(comm_neo, Mammal > 0)$Species %>% gsub("_", " ", .)

mampred_sp <- c(mampred_sp_afr,
                mampred_sp_ind,
                mampred_sp_mad,
                mampred_sp_neo)

# All species
all_sp_afr <- comm_afr$Species %>% gsub("_", " ", .)
all_sp_ind <- comm_ind$Species %>% gsub("_", " ", .)
all_sp_mad <- comm_mad$Species %>% gsub("_", " ", .)
all_sp_neo <- comm_neo$Species %>% gsub("_", " ", .)

all_sp <- c(all_sp_afr,
            all_sp_ind,
            all_sp_mad,
            all_sp_neo)


# Make a dataframe of species names to output
mam_checklist <- tibble(species = all_sp,
                        region = rep(c("African",
                                       "Indomalayan",
                                       "Madagascan",
                                       "Neotropical"), 
                                     times = lapply(list(all_sp_afr,
                                                         all_sp_ind,
                                                         all_sp_mad,
                                                         all_sp_neo),length)),
                        mampred = ifelse(all_sp %in% mampred_sp, "yes", "no"))
write.csv(mam_checklist, "mam_checklist.csv")


# Get Web of Science queries -------------------------------------------
sid <- auth(username = NULL, password = NULL)

wos_terms <- paste("TS = (\"", all_sp, "\" AND diet)", sep = "")

wos_res <- pull_wos(wos_terms[1], sid = sid)

wos_x <- wos_res[[1]][0,]
wos_x$species <- c()
wos_x$mammal_predator <- c()
wos_x$counter <- c()


for(i in 1:length(wos_terms)){
  dat <- pull_wos(wos_terms[i], sid = sid)[[1]]
  if(nrow(dat) == 0) next
  dat$species <- all_sp[i]
  dat$mammal_predator <- ifelse(all_sp[i] %in% mampred_sp, 1, 0)
  dat <- dat[sample(1:nrow(dat)),]
  dat$counter <- 1:nrow(dat)
  
  wos_x <- rbind(wos_x, dat)
  print(paste(i, "of", length(wos_terms)))
}



wos_afr <- wos_x %>% filter(species %in% all_sp_afr)
wos_ind <- wos_x %>% filter(species %in% all_sp_ind)
wos_mad <- wos_x %>% filter(species %in% all_sp_mad)
wos_neo <- wos_x %>% filter(species %in% all_sp_neo)

write.csv(wos_afr, file = "wos_afr.csv")
write.csv(wos_ind, file = "wos_ind.csv")
write.csv(wos_mad, file = "wos_mad.csv")
write.csv(wos_neo, file = "wos_neo.csv")


