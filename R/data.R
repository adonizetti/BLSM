#' Example Adjacency Matrix
#' Adjacency matrix of a 10 nodes random network for testing purposes
#' 
#' @name example_adjacency_matrix
#' @docType data
#' @format A binary adjacency matrix representing links between nodes.
"example_adjacency_matrix"

#' Example Weights Matrix
#' 
#' "BLSM weights" matrix of a 10 nodes random network for testing purposes
#' 
#' @name example_weights_matrix
#' @docType data
#' @format A matrix containing positive weights for all pairs of nodes. 
#' 
#' Given a couple of nodes, a weight expresses the importance of the distance between the 
#' coordinates associated to the two nodes in the latent space in terms of the overall likelihood of the graph. 
#' For this reason, even missing links must have a coefficient, otherwise the relative positioning of disconnected nodes
#' would have no effect at all on the graph likelihood.
#'  
#' The exact probability equation is described in \link[BLSM]{BLSM}, as well as the notation used.
#' 
#' A few examples:
#' \itemize{
#' \item for unweighted networks, the "BLSM weights" matrix has all the values set to 1. 
#' \item if two nodes share a strong connection, then
#' the weight coefficient should be greater than 1 so that their positions in the latent space will be closer than they would be in an unweighted framework. 
#' \item if two nodes share a weak connection, a coefficient smaller than 1 will allow the latent coordinates to be pretty far from each other even though the nodes are connected. 
#' }
"example_weights_matrix"