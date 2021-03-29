library("tidyverse")
library("bipartite")
library("igraph")
#install.packages("intergraph")
library("intergraph")
library("ggnetwork")

load(file = "web_list.RData")
pred_tidy <- read.csv("pred_tidy.csv", header = T)[,-1] %>% tibble()

# Look at an example network
i <- "neo"
j <- 5

# n <- graph_from_incidence_matrix(web_list[[i]][[j]], directed = T, mode = "out")
# 
# n <- 
# 
# n <- network(x = web_list[[i]][[j]],
#              directed = TRUE,
#              matrix.type = "incidence")
# 
# plot(n)
# 
# asNetwork()

set.seed(5)
n <- as.network(web_list[[i]][[j]], directed = T, bipartite = T)

is_predator <- ifelse(network::get.vertex.attribute(n, attrname = "vertex.names") %in% pred_tidy$consumer_sp,
                      "predator", "prey")
n <- network::set.vertex.attribute(n, # the name of the network object
                     "Level", # the name we want to reference the variable by in that object
                     is_predator # the value we are giving that variable
) 

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodelabel_repel(aes(label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  geom_nodes(aes(colour = Level), size = 8) +
  theme_blank()



# Testing
for(i in reg[-1]){
  for(j in names(web_list[[i]])){
    web <- web_list[[i]][[j]]
    
    if(dim(web)[1] == 0){
      next()
    } 
    
    if(colnames(web) %in% rownames(web)){
      stop()
    }
  }
}
