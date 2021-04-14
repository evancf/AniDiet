#install.packages("googlesheets4")
library("googlesheets4")
library("tidyverse")


# Pull in data from google sheets -----------
neo_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1jJuUNBwvfQ7RzHK035MUOOkwAtMimeXvOnZhBY1FAjk/edit#gid=581553521",
                      sheet = "Neotropical") %>% 
  mutate(region = "Neotropical")
afr_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1jJuUNBwvfQ7RzHK035MUOOkwAtMimeXvOnZhBY1FAjk/edit#gid=581553521",
                      sheet = "African") %>% 
  mutate(region = "African")
mad_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1jJuUNBwvfQ7RzHK035MUOOkwAtMimeXvOnZhBY1FAjk/edit#gid=581553521",
                      sheet = "Madagascan") %>% 
  mutate(region = "Madagascan")
ind_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1jJuUNBwvfQ7RzHK035MUOOkwAtMimeXvOnZhBY1FAjk/edit#gid=581553521",
                      sheet = "Indomalayan") %>% 
  mutate(region = "Indomalayan")
ind_dat$mammal_predator <- unlist(ind_dat$mammal_predator) %>% as.numeric()


# Bind rows
paper_data <- bind_rows(neo_dat, afr_dat, mad_dat, ind_dat)

# Removing duplicates --------------

# Want to keep just ones where there's useful dat

paper_data <- paper_data[order(paper_data$`contains useful data? (y/n)`, decreasing = T),]
paper_data <- paper_data[!duplicated(select(paper_data, ut, species)),]

paper_data$`contains useful data? (y/n)` <- gsub("yes", "y", paper_data$`contains useful data? (y/n)`)

# Outputing as csv file -----------------------

write.csv(file = "diet_paper_classification.csv", paper_data)



