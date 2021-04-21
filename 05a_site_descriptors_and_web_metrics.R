library("tidyverse")


# Get the Rowan community data
tmp <- tempfile(fileext = ".xlsx")
sheet <- "Community Data"
download.file(url = "https://www.pnas.org/highwire/filestream/902133/field_highwire_adjunct_files/1/pnas.1910489116.sd01.xlsx", 
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



comm_data <- bind_rows(comm_afr, comm_ind, comm_mad, comm_neo)





# Read local web metrics
local_web_metrics <- read.csv(file = "local_web_metrics.csv")[,-1] %>% tibble()

local_web_metrics

# Read in site descriptors
site_descriptors <- read.csv(file = "./Site descriptor data/Rowan_2020_Community_Env_partialRoad.csv")[,-1] %>% tibble()

site_descriptors_afr <- site_descriptors %>% 
  filter(Realm == "Afrotropic") %>% left_join(select(comm_data, 1:2), 
                               by = c("Full.Name" = "Full Name")) %>% 
  mutate(site = Community, .before = Full.Name) %>% 
  select(-Full.Name, -Community)
site_descriptors_ind <- site_descriptors %>% 
  filter(Realm == "IndoMalay") %>% left_join(select(comm_data, 1:2), 
                                              by = c("Full.Name" = "Full Name")) %>% 
  mutate(site = Community, .before = Full.Name) %>% 
  select(-Full.Name, -Community)
site_descriptors_mad <- site_descriptors %>% 
  filter(Realm == "Madagascar") %>% left_join(select(comm_data, 1:2), 
                                              by = c("Full.Name" = "Community")) %>% 
  mutate(site = Full.Name, .before = Country) %>% 
  select(-Full.Name, -"Full Name")

site_descriptors_neo <- site_descriptors %>% 
  filter(Realm == "Neotropic") %>% left_join(select(comm_data, 1:2), 
                                              by = c("Full.Name" = "Full Name")) %>% 
  mutate(site = Community, .before = Full.Name) %>% 
  select(-Full.Name, -Community)


site_descriptors <- bind_rows(site_descriptors_afr,
                              site_descriptors_ind,
                              site_descriptors_mad,
                              site_descriptors_neo) %>% 
  select(-Realm)


# Join network metric and site descriptors together ----------------------

local_web_data <- local_web_metrics %>% 
  left_join(site_descriptors, by = "site")

# Write this out --------------------------

write.csv(file = "local_web_data.csv", local_web_data)


# Testing

library("tidyverse")

# Read in data
local_web_data <- read.csv("local_web_data.csv")[,-1] %>% tibble()

# Explore the data a little
local_web_data

# Look at the structure
str(local_web_data)

# There are many column names - let's look at them all
colnames(local_web_data)

# Filter to just your dataset
my_data <- filter(local_web_data,
                  reg == "neo")

# Produce a scatterplot 
ggplot(my_data, aes(x = Lat,
                    y = connectance)) +
  geom_point() +
  geom_smooth(method = "lm")

# Fit a linear regression model
mod <- lm(connectance ~ Lat, data = my_data)

# Examine the model output
summary(mod)





