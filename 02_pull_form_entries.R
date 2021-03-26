#install.packages("googlesheets4")
library("googlesheets4")
library("tidyverse")
library("taxize")


# Pull in data from google sheets -----------
pred_dat <- read_sheet("https://docs.google.com/spreadsheets/d/1CEhz9fN6OzqyDeRIYdbyt5LsOsZUY1MsmAuuNlqgVKk/edit#gid=1739630829")

# Fix up colnames 
pred_dat <- rename(pred_dat, 
                   researcher_id = "What's your researcher ID?",
                   data_source = "What is the data source?",
                   consumer_sp = "What is the consumer species?",
                   diet_items = "What items are reported in its diet?",
                   items_not_in_diet = "Are there items reported to NOT be in its diet?",
                   includes_scavanged = "Diet items includes scavenged items?",
                   replaces_previous = "Does this replace a form entry with an error?")

colnames(pred_dat) <- tolower(colnames(pred_dat))


# Deal with ones that replace a previous -----------

replaced_rows <- which(!is.na(pred_dat$replaces_previous))
cols_to_match <- c("researcher_id", "data_source", "consumer_sp")
rows_to_remove <- c()

# Very clunky way to match across columns - there must be a better way!
for(i in replaced_rows){
  rows_to_remove <- c(rows_to_remove,
                      which(apply(pred_dat[, cols_to_match], 
                                  1, function(x) paste(x, collapse = " ")) %in%
                              apply(pred_dat[i,cols_to_match], 
                                    1, function(x) paste(x, collapse = " ")))[1])
  
}
# Now actually remove these rows
pred_dat <- pred_dat[-rows_to_remove, ]
# No longer need the column of replaces previous
pred_dat <- pred_dat %>% select(-replaces_previous)


# Split up the comma separated prey items -----------

# First, there are some issues where commas aren't used
pred_dat$diet_items <- gsub(" and ", ", ", pred_dat$diet_items)
pred_dat$items_not_in_diet <- gsub(" and ", ", ", pred_dat$items_not_in_diet)

pred_dat$diet_items <- gsub("\n", ", ", pred_dat$diet_items, fixed = T)
pred_dat$items_not_in_diet <- gsub("\n", ", ", pred_dat$items_not_in_diet, fixed = T)


pred_long <- tibble(timestamp = c(),
                    researcher_id = c(),
                    data_source = c(),
                    consumer_sp = c(),
                    resource_sp = c(),
                    consumed = c(),
                    includes_scavanged = c(),
                    notes = c())

for(i in 1:nrow(pred_dat)){
  diet <- strsplit(pred_dat$diet_items[i], ",", fixed = T) %>% unlist() %>% trimws()
  not_diet <- strsplit(pred_dat$items_not_in_diet[i], ",", fixed = T) %>% unlist() %>% trimws() %>% na.omit()
  
  pred_long <- pred_long %>% 
    rbind(tibble(timestamp = pred_dat$timestamp[i],
                 researcher_id = pred_dat$researcher_id[i],
                 data_source = pred_dat$data_source[i],
                 consumer_sp = pred_dat$consumer_sp[i],
                 resource_sp = c(diet, not_diet),
                 consumed = rep(c(1,0), 
                                times = c(length(diet),
                                          length(not_diet))),
                 includes_scavanged = pred_dat$includes_scavanged[i],
                 notes = pred_dat$notes[i]))
  
  if(i %in% round(seq(1, nrow(pred_dat), length.out = 10))){
    print(paste(i, "of", nrow(pred_dat)))
  }
}



# Fix up some issues with diet item names
unique(sort(pred_long$resource_sp))
pred_long$resource_sp <- pred_long$resource_sp %>% tolower()

# This is a little tricky because the order matters
descriptors_to_remove <- c("adult ",
                           " fawns",
                           " fawn",
                           "domestic ", "wild ",
                           "large ammounts of ",
                           "large ", 
                           " cf. ", "occasional ", "other ",
                           "particularly ",
                           "rarely attack",
                           "calves ",
                           " eggs",
                           "small ", "specialize in ",
                           "young ")

for(i in descriptors_to_remove){
  pred_long$resource_sp <- gsub(i, "", pred_long$resource_sp, fixed = T)
}

words_to_remove <- c("a", "sp.", "spp.",
                     "sp", "spp")

inds <- which(sapply(strsplit(pred_long$resource_sp, " "), length) != 1)

for(i in words_to_remove){
  pred_long$resource_sp[inds] <- ifelse(word(pred_long$resource_sp[inds], 1) == i,
                                        word(pred_long$resource_sp[inds], 2),
                                        pred_long$resource_sp[inds])
  
  pred_long$resource_sp[inds] <- ifelse(word(pred_long$resource_sp[inds], 2) == i,
                                        word(pred_long$resource_sp[inds], 1),
                                        pred_long$resource_sp[inds])
}

pred_long$resource_sp <- pred_long$resource_sp %>% trimws()

pred_long

# Translate from common name to scientific name -----------
c2s_names <- comm2sci(unique(sort(pred_long$resource_sp)))
c2s_names1 <- unlist(c2s_names)
c2s_names1[which(nchar(c2s_names1) > 0)] # Check if these seem right

c2s_names1 <- c2s_names1 %>% bind_rows() %>% t() %>% # Why on earth is this not easier...
  as.data.frame() %>% 
  tibble(resource_sp = rownames(.)) %>% 
  rename(new_resource_sp = V1)

pred_long <- pred_long %>% left_join(c2s_names1) %>% 
  mutate(resource_sp = ifelse(is.na(new_resource_sp),
                              resource_sp,
                              new_resource_sp)) %>% 
  select(-new_resource_sp)


# Clean things up a little bit

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

pred_long$resource_sp <- pred_long$resource_sp %>% firstup()


# Output to a csv file 

write.csv(pred_long, "pred_long.csv")


# Some summary stuff

# How many predators are prey of other species?
pred_long$consumer_sp[which(pred_long$consumer_sp %in% pred_long$resource_sp)] %>% unique()
