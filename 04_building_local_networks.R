library("tidyverse")

# Get community data ------------------------------------
# This is from Rowan et al. 2019 
# https://www.pnas.org/content/117/3/1559/tab-figures-data

# Will download to this temp file, which gets removed with unlink below.
tmp <- tempfile(fileext = ".xlsx")
sheet <- "PA Data"
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

comm_list <- list(comm_afr, comm_ind, comm_mad, comm_neo)

# Assemble networks at each site ------------

# Get a metaweb (including a long version)
pred_long <- read.csv("pred_long.csv")[,-1] %>% tibble()

pred_long$consumer_sp <- pred_long$consumer_sp %>% word(start = 1, end = 2)
pred_long$resource_sp <- pred_long$resource_sp %>% word(start = 1, end = 2)


pl <- pred_long[,c("consumer_sp", "resource_sp")]
pl <- pl[which(sapply(strsplit(pl$resource_sp, " "), length)==2),]
pl$consumer_gen <- word(pl$consumer_sp, 1)
pl$resource_gen <- word(pl$resource_sp, 1)

metaweb <- table(pl$consumer_sp[row(pl[-1])], unlist(pl[-1])) %>% t()



# 

web_list <- list()

reg <- c("afr", "ind", "mad", "neo")

for(i in 1:length(reg)){
  web_list[[reg[i]]] <- list()
  
  comm <- comm_list[[i]][,-(1:3)]
  
  for(j in colnames(comm[,-1])){
    
    local_sp <- comm$Species[which(comm[,j] == 1)] %>% gsub("_", " ", . , fixed = T)
    local_gen <- comm$Species[which(comm[,j] == 1)] %>% gsub("_", " ", . , fixed = T) %>% word(1)
    
    sp_level_inds <- which(pl$consumer_sp %in% local_sp & pl$resource_sp %in% local_sp)
    gen_level_inds <- which(word(pl$consumer_sp, 1) %in% local_gen & word(pl$resource_gen) %in% local_gen)
    
    inds <- unique(c(sp_level_inds, gen_level_inds))
    
    local_pl <- pl[inds, 
                   c("consumer_sp", "resource_sp")]
    
    web_list[[reg[i]]][[j]] <- table(local_pl$consumer_sp[row(local_pl[-1])], 
                                     unlist(local_pl[-1])) %>% t()
  }
  print(paste(i, "out of", length(reg)))
}


save(file = "web_list.RData", web_list)


