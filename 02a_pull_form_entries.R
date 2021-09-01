# This script pulls diet data from our new assembly of literature data.
#install.packages("googlesheets4")
library("googlesheets4")
library("tidyverse")


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

# Check out if there are any other issues that should be addressed at this stage
unique(sort(pred_long$resource_sp))

# Output to a csv file 
write.csv(pred_long, "pred_long.csv")
