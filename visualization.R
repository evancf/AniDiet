library("tidyverse")
library("bipartite")
library("igraph")
#install.packages("intergraph")
library("intergraph")
library("ggnetwork")

load(file = "web_list.RData")

# Look at an example network
i <- "neo"
j <- 3
n <- as.network(web_list[[i]][[j]], directed = T, bipartite = T)

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodelabel_repel(aes(label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  geom_nodes(color = "black", size = 6) +
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
