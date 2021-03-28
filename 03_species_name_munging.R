library("tidyverse")
library("taxize")

# Pull in the pred_long data
pred_long <- read.csv("pred_long.csv")[,-1] %>% tibble()

# Fix up some issues with diet item names
unique(sort(pred_long$resource_sp))
pred_long$resource_sp <- pred_long$resource_sp %>% tolower()

# This is a little tricky because the order matters
descriptors_to_remove <- c("adult ",
                           " fawns",
                           " fawn",
                           "domestic ", "wild ",
                           "large ammounts of ",
                           "ammounts of ", 
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


# Translate from common name to scientific name -----------

# c2s_names <- comm2sci(unique(sort(pred_long$resource_sp)))
# comm2sci_fixes <- unlist(c2s_names)
# comm2sci_fixes[which(nchar(comm2sci_fixes) > 0)] # Check if these seem right
# 
# comm2sci_fixes <- comm2sci_fixes %>% bind_rows() %>% t() %>% # Why on earth is this not easier...
#   as.data.frame() %>% 
#   tibble(resource_sp = rownames(.)) %>% 
#   rename(new_resource_sp = V1)
# 
# write.csv(comm2sci_fixes, "comm2sci_fixes.csv")

comm2sci_fixes <- read.csv("comm2sci_fixes.csv")[,-1] %>% tibble()


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


# Manually change problems -----------

# These next couple steps are for creating a csv file to manually add 
# corrections, filling in new ones to correct
new_manual_resource_sp_changes <- tibble(resource_sp = sort(unique(pred_long$resource_sp)),
                                         corrected_resource_sp = NA)
prev_manual_resource_sp_changes <- read.csv("manual_resource_sp_changes.csv", row.names = 1) %>% tibble()

manual_resource_sp_changes <- rbind(prev_manual_resource_sp_changes,
                                    filter(new_manual_resource_sp_changes, 
                                           !resource_sp %in% prev_manual_resource_sp_changes$resource_sp))

write.csv(manual_resource_sp_changes, "manual_resource_sp_changes.csv")


# Now go to that file to correct things, save, and proceed

manual_resource_sp_changes <- read.csv("manual_resource_sp_changes.csv", row.names = 1) %>% tibble()

# Do the switcheroo
pred_long <- pred_long %>% left_join(manual_resource_sp_changes) %>% 
  mutate(resource_sp = ifelse(is.na(corrected_resource_sp),
                              resource_sp,
                              corrected_resource_sp)) %>% 
  select(-corrected_resource_sp)



# Some summary stuff -----------

# How many predators are prey of other species?
pred_long$consumer_sp[which(pred_long$consumer_sp %in% pred_long$resource_sp)] %>% unique()
