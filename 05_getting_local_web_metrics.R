library("tidyverse")
library("bipartite")
library("igraph")
#install.packages("intergraph")
library("intergraph")
library("ggnetwork")


# Load the list of networks
load(file = "web_list.RData")

# Decide what species level indices
sp_indices <- c("degree", "normalised degree")

# Get vectors of names for the network level and species level metrics
web_names <- names(networklevel(Safariland, weighted=F))
consumer_names <- paste("consumer", names(specieslevel(Safariland, index = sp_indices)[[1]]), sep="_")
resource_names <- paste("resource", names(specieslevel(Safariland, index = sp_indices)[[2]]), sep="_")

# Prepare a tibble to store all these data
regs <- rep(names(web_list), 
            times = unlist(lapply(web_list, 
                          function(x) length(names(x)))))

sites <- lapply(web_list, names) %>% unlist() %>% as.vector()

local_web_metrics <- tibble(reg = regs,
                            site = sites)

metric_colnames <- c(web_names,
                     consumer_names,
                     resource_names) 

local_web_metrics <- bind_cols(local_web_metrics,
                               matrix(NA, nrow = nrow(local_web_metrics),
                                      ncol = length(metric_colnames)) %>% 
                                 as.data.frame() %>% 
                                 rename_at(vars(colnames(.)), ~ metric_colnames) %>% 
                                 tibble())
local_web_metrics <- local_web_metrics %>% mutate_at(metric_colnames, as.numeric)


# Now use a loop to extract metrics (a tricky lapply could do this...)

my_median <- function(x){
  x[x == "NaN"] <- NA
  
  median(x, na.rm = T)
}

reg <- c("afr", "ind", "mad", "neo")
for(i in reg){
  for(j in names(web_list[[i]])){
    
    # Get the focal web
    web <- web_list[[i]][[j]]
    
    # Skip cases where there aren't at least two species
    if(all(dim(web)<2)){
      next()
    }
    
    # Get the network- and species-level metrics from the bipartite package
    web_metrics <- networklevel(web, weighted=F)
    sp_metrics <- specieslevel(web, index = sp_indices)
    
    # Get median values across species
    consumer_metrics <- sp_metrics[[1]] %>% apply(., 2, my_median)
    resource_metrics <- sp_metrics[[2]] %>% apply(., 2, my_median)
    
    # Set up to add these to the big tibble
    ind <- which(local_web_metrics$reg == i & local_web_metrics$site == j)
    
    all_metrics <- as.numeric(c(web_metrics, consumer_metrics, resource_metrics))
    all_metrics[all_metrics == "NaN"] <- NA
    
    local_web_metrics[ind, metric_colnames] <- as.list(all_metrics)
    
    print(j)
  }
  print(i)
}

# Write this out to csv
write.csv(file = "local_web_metrics.csv", local_web_metrics)



#
dim(local_web_metrics)

sum(complete.cases(local_web_metrics))

plot(connectance ~ factor(reg), data = local_web_metrics)

plot(local_web_metrics$'number of compartments' ~ local_web_metrics$connectance)
