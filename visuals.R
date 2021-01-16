
if(!require("igraph"))install.packages("igraph")
library(igraph)

build_network <- function(distance){
  # build the graph object
  network <- graph_from_adjacency_matrix(distance, weighted=T, mode="undirected", diag=F)  
  return(network)
}
############################3
if(!require("RColorBrewer"))install.packages("RColorBrewer")
library(RColorBrewer)

get_trait_from_id <- function(data, names){
  ill <- c()
  for (id in names){
    ill <- c(ill, data$Trait[data$Study.id == id][1])
  }
  return(ill)
}

plot_cfg <- function(network, dist_names){
  ill <- get_trait_from_id(data, dist_names)
  colr <- brewer.pal(nlevels(as.factor(ill)), "Set2")
  #colr <- colorRampPalette(colr)(25)
  my_color <- colr[as.numeric(as.factor(ill))]
  
  #ceb <- cluster_edge_betweenness(network)
  #plot(ceb, network)
  cfg <- cluster_fast_greedy(as.undirected(network))
  plot(cfg, as.undirected(network),
       vertex.size=12,
       vertex.color=my_color, 
       vertex.label.cex=0.7,
       vertex.label.color="black")
       #vertex.frame.color="transparent")
  
  # title and legend
  text(0,1.5,"Study.id correlation",col="black", cex=1.5)
  legend(x=-1.5, y=1.5, 
         legend=paste( levels(as.factor(ill)), sep=""), 
         col = colr, 
         bty = "n", pch=20 , pt.cex = 2, cex = 1,
         text.col="black" , horiz = F)
}

plot_network <- function(network, dist_names){
  ill <- get_trait_from_id(data, dist_names)
  colr <- colorRampPalette(brewer.pal(8, "Set2"))(nlevels(as.factor(ill)))
  colr <- colorRampPalette(colr)(25)
  my_color <- colr[as.numeric(as.factor(ill))]
  
  plot(network,
       vertex.size=12,
       vertex.color=my_color, 
       vertex.label.cex=0.7,
       vertex.label.color="black")
  #vertex.frame.color="transparent")
  
  # title and legend
  text(0,1.5,"Study.id correlation",col="black", cex=1.5)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", title = "Colors that represent diseases",
         legend=paste( levels(as.factor(ill)), sep=""), 
         col = colr, 
         bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
         text.col="black" , horiz = F)
}
#network

if(!require(ggraph))install.packages("ggraph")
library(ggraph)
###########
n <- ncol(dist1)
ggraph(network, layout="fr") + 
  #geom_edge_density(edge_fill="#69b3a2") +
  geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.3) +
  geom_node_point(aes(size=n, alpha=n)) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(rep(1,4), "cm")
  ) 

################


#dist1 <- dist[1:10, 1:10]
plot_network_with_deg<- function(distance){
  network <- graph_from_adjacency_matrix(distance, weighted = T, mode="undirected", diag=F)
  deg <- degree(network, mode = "all")
  plot(network, vertex.size= deg*6, vertex.color=rgb(0.1,0.7,0.8,0.5))
}
#dist1

#K-core decomposition
###########################
plot_K_core <- function(network){
  kc<- coreness(network, mode="all")
  plot(network,  vertex.size=kc*6, vertex.label=kc, vertex.color=colr[kc])

}

get_split_ids <- function(network){
  #components(network)
  split_names<- split(names(V(network)), components(network)$membership)
  #split_names
  
  ids <- list()
  i <- 1
  for(id in split_names){
    if(length(id) > 1){
      ids[[i]]<- id
      i <- i+1
    }
  }
  return(ids)
  
}

get_trait_list <- function(ids, data){
  trait_list <- list()
  j <- 1
  for (i in 1:length(ids)){
    traits <- get_trait_from_id(data, ids[[i]])
    trait_list[[j]] <- traits
    j <- j+1
  }
  return(trait_list)
}


get_max_id <- function(ids){
  max <- c()
  for(i in 1:length(ids))
    if (length(ids[[i]]) > length(max) ){
      max <- ids[[i]]
    }
  return(max)
}
####################################
####################################
####################################
#dist <- readRDS("../Desktop/stud/7sem/kursinis2/study_distance.rds")
#data <- read.csv("../Desktop/stud/7sem/kursinis2/asc_stud_coh.csv", as.is = TRUE, sep = ",")
#sum(dist[5,] != 0)
#length(dist[1,])
#dist[1,] != 0

#table(dist)
#kurie atitinka > 30%
#kurie atitinka > 50%

#proc <- ncol(dist) *30 /100
#proc

#arr_id <- c()
#for( id in rownames(dist)){
#  num <- sum(dist[id,] != 0)
#  if (num > proc){
#    arr_id <- c(arr_id, id)    
#  }
#}
#dist1 <- dist[arr_id,arr_id]
#dist1 <- dist[1:50, 1:50]


############
############

fisher <- readRDS("../Desktop/stud/7sem/kursinis2/adjusted_pval_more_1_odds_5.rds")
odds <- readRDS("../Desktop/stud/7sem/kursinis2/adjusted_pval_less_1_odds_5.rds")

#dist1
#dist1[1:5, 1:5]

dist <- fisher[1:200,1:200]

#####1 make network
net <- build_network(dist)
#####2 plot K core
plot_K_core(net)
#####3 get split ID's -> ids that has similarities
ids <- get_split_ids(net)
ids
#get cluster with max study.id
max<- get_max_id(ids)
dist <- dist[colnames(dist) %in% max, colnames(dist) %in% max]
net <- build_network(dist)
#plot_cfg(net,colnames(dist))

plot_network(net, colnames(dist))
#plot_network_with_deg(dist)
#plot(net, layout=layout.circle, main="circle")
#plot(net, layout=layout.fruchterman.reingold, main="fruchterman.reingold")
#####4 check traits from them
trait_list <- get_trait_list(ids, data)
#trait_list
#components(net)
max_traits <- get_trait_from_id(data, max)
unique(max_traits)
####

unique(max_traits)
length(max_traits)
length(unique(max_traits))
#par(bg="black")
# plot it

############
#########
######
dist <- odds[1:200,1:200]
#####1 make network
net <- build_network(dist)
#####2 plot K core
plot_K_core(net)
#####3 get split ID's -> ids that has similarities
ids <- get_split_ids(net)
ids
max<- get_max_id(ids)
dist <- dist[colnames(dist) %in% max, colnames(dist) %in% max]
net <- build_network(dist)
#plot_cfg(net,colnames(dist))

plot_network(net, colnames(dist))
#plot_network_with_deg(dist)
#plot(net, layout=layout.circle, main="circle")
#plot(net, layout=layout.fruchterman.reingold, main="fruchterman.reingold")
#####4 check traits from them
trait_list <- get_trait_list(ids, data)
#trait_list
#components(net)
max_traits <- get_trait_from_id(data, max)
unique(max_traits)
####

unique(max_traits)
length(max_traits)
length(unique(max_traits))
#par(bg="black")
# plot it


#id_array <- get_trait_from_id(data, ids)

#length(id_array)
#unique(id_array)
#dim(dist1)

#dist2 <- dist1[colnames(dist1) %in% ids ,colnames(dist1) %in% ids ]
#dim(dist2)
#network2 <- graph_from_adjacency_matrix(dist[1:100,1:100], weighted=T, mode="undirected", diag=F)

#kc<- coreness(network2, mode="all")
#plot(network2,  vertex.size=kc*6, vertex.label=kc, vertex.color=colr[kc])


#split_names<- split(names(V(network2)), components(network2)$membership)


##########################
############################
#############################

#library(tidyverse)
#library(viridis)
#library(patchwork)
#library(hrbrthemes)
#library(ggraph)
#library(igraph)
#if(!require("networkD3"))install.packages("networkD3")
#library(networkD3)


# Load researcher data
#dataUU <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyUndirectedUnweighted.csv", header=TRUE)

# Transform the adjacency matrix in a long format
#connect <- dataUU %>% 
#  gather(key="to", value="value", -1) %>%
#  na.omit()

# Number of connection per person
#c( as.character(connect$from), as.character(connect$to)) %>%
#  as.tibble() %>%
#  group_by(value) %>%
#  summarize(n=n()) -> coauth
#colnames(coauth) <- c("name", "n")

# NetworkD3 format
#graph=simpleNetwork(connect)

# Plot
#simpleNetwork(connect,     
#              Source = 1,                 # column number of source
#              Target = 2,                 # column number of target
#              height = 880,               # height of frame area in pixels
#              width = 1980,
#              linkDistance = 10,         # distance between node. Increase this value to have more space between nodes
#              charge = -4,              # numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)
#              fontSize = 5,              # size of the node names
#              fontFamily = "serif",       # font og node names
#              linkColour = "#666",        # colour of edges, MUST be a common colour for the whole graph
#              nodeColour = "#69b3a2",     # colour of nodes, MUST be a common colour for the whole graph
#              opacity = 0.9,              # opacity of nodes. 0=transparent. 1=no transparency
#              zoom = T                    # Can you zoom on the figure?
#)
