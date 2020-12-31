## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(gretel)
BuchDarrah19

## ------------------------------------------------------------------------
opt_gpv(BuchDarrah19, source = 1, target = 5, p = 0)

## ------------------------------------------------------------------------
opt_gpv(BuchDarrah19, source = 1, target = 5, p = 1)

## ------------------------------------------------------------------------
opt_gpv(BuchDarrah19, source = 1, target = 5, p = Inf)

## ------------------------------------------------------------------------
optimal_path <- opt_ppv(BuchDarrah19, source = 1, target = 5, odds_scale = 4)
print(optimal_path)

## ------------------------------------------------------------------------
prob_path_value <- ppv(BuchDarrah19, path = optimal_path, odds_scale = 4)
transmission_odds <- prob_path_value/4
print(prob_path_value)
print(transmission_odds)

## ------------------------------------------------------------------------
generate_proximities(BuchDarrah19, mode = "ogpv", p = 1)

generate_proximities(BuchDarrah19, mode = "oppv", odds_scale = 4)

generate_proximities(BuchDarrah19, mode = "sconductivity")

## ------------------------------------------------------------------------
# Suppose we wish to calculate the betweenness centrality of node 3 in 
# the example sociomatrix 'BuchDarrah19'

# There are n-1 choose 2 pairs of nodes for which neither the source nor the target is 3
# Since n = 5 in this case, there are 6 pairs of nodes to consider.
# Therefore, node 3 mediates at most 6 shortest paths.
# (If this were a directed graph, this number would be 12)

all_paths <- all_opt_gpv(BuchDarrah19, p = 1)
paths_mediated <- 0
# We can get away with just looking at half the sociomatrix (below) because this
# one is undirected.
for(i in 1:5){
  for(j in (i+1:5)){
    shortest_ij <- unpack(all_paths[[i]], i, j)
    if(3 %in% shortest_ij) paths_mediated <- paths_mediated + 1
  }
}

# Print Betweenness Centrality
print(paths_mediated/6)


## What about if p = 0?
all_paths <- all_opt_gpv(BuchDarrah19, p = 0)
paths_mediated <- 0
# We can get away with just looking at half the sociomatrix (below) because this
# one is undirected.
for(i in 1:5){
  for(j in (i+1:5)){
    shortest_ij <- unpack(all_paths[[i]], i, j)
    if(3 %in% shortest_ij) paths_mediated <- paths_mediated + 1
  }
}

# Print Betweenness Centrality
print(paths_mediated/6)


