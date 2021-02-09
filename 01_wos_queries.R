library("tidyverse")

install.packages("wosr")
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

# All species
all_sp_afr <- comm_afr$Species %>% gsub("_", " ", .)
all_sp_ind <- comm_ind$Species %>% gsub("_", " ", .)
all_sp_mad <- comm_mad$Species %>% gsub("_", " ", .)
all_sp_neo <- comm_neo$Species %>% gsub("_", " ", .)



# Get Web of Science queries -------------------------------------------
sid <- auth(username = NULL, password = NULL)

paste("TS = (\"", all_sp_ind, "\" AND diet)", sep = "")
query_wos("TS = (\"Puma concolor\" AND diet)", sid = sid)
