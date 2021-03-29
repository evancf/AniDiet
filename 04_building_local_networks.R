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
pred_tidy <- read.csv("pred_tidy.csv")[,-1] %>% tibble()

pred_tidy$consumer_sp <- pred_tidy$consumer_sp %>% word(start = 1, end = 2)
pred_tidy$resource_sp <- pred_tidy$resource_sp %>% word(start = 1, end = 2)


pl <- pred_tidy[,c("consumer_sp", "resource_sp")]
pl <- pl[which(sapply(strsplit(pl$resource_sp, " "), length)==2),]
pl$consumer_gen <- word(pl$consumer_sp, 1)
pl$resource_gen <- word(pl$resource_sp, 1)

pl$intx_sp <- paste(pl$consumer_sp, pl$resource_sp)
pl$intx_gen <- paste(pl$consumer_gen, pl$resource_gen)

# Get a metaweb (not sure if useful at this point)
metaweb <- table(pl$consumer_sp[row(pl[-1])], unlist(pl[-1])) %>% t()


# Now develop a web at every local site, given species presence
web_list <- list()

reg <- c("afr", "ind", "mad", "neo")

get_combn_df <- function(x){
  dat <- cbind(combn(x, 2), combn(x, 2)[2:1,]) %>% t() %>% as.data.frame()
  colnames(dat) <- c("consumer_sp", "resource_sp")
  dat$sp_combn <- paste(dat$consumer_sp, dat$resource_sp)
  dat$gen_combn <- paste(word(dat$consumer_sp, 1), word(dat$resource_sp, 1))
  dat
}

for(i in 1:length(reg)){
  web_list[[reg[i]]] <- list()
  
  comm <- comm_list[[i]][,-(1:3)]
  
  for(j in colnames(comm[,-1])){
    
    local_sp <- comm$Species[which(comm[,j] == 1)] %>% gsub("_", " ", . , fixed = T)
    
    combn_df <- get_combn_df(local_sp)
    
    sp_level_inds <- which(combn_df$sp_combn %in% pl$intx_sp)
    gen_level_inds <- which(combn_df$gen_combn %in% pl$intx_gen)
    
    inds <- unique(c(sp_level_inds, gen_level_inds))
    
    local_pl <- combn_df[inds, 
                   c("consumer_sp", "resource_sp")]
    
    web_list[[reg[i]]][[j]] <- table(local_pl$consumer_sp[row(local_pl[-1])], 
                                     unlist(local_pl[-1])) %>% t()
    #print(j)
  }
  print(paste(i, "out of", length(reg)))
}


save(file = "web_list.RData", web_list)


