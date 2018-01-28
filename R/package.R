#' @name BLSM
#' @title Bayesian Latent Space Model
#' @docType package
#' 
#' @details R package allowing the computation of a Bayesian Latent Space Model for Complex Networks.
#' 
#' Given a network characterized by its adjacency \eqn{Y} and weights \eqn{W} matrices, the model assigns a binary random 
#' variable to each tie: e.g., \eqn{Y_ij} is related to the tie between nodes \eqn{i} and \eqn{j} and its value is 1
#' if the tie exists, 0 otherwise. 
#'  
#' The model assumes the independence of \eqn{Y_ij | x_i,x_j\alpha}, where \eqn{x_i} and \eqn{x_j} are the coordinates
#' of the nodes in the multidimensional latent space and \eqn{\alpha} is an additional parameter such that 
#' \eqn{P(Y_ij = 1) = \alpha - ||x_i -x_j||}.
#' 
#' The latent space coordinates are estimated by maximizing the overall likelihood induced by the above equation. 
#' 
#' The output of the model allows the user to inspect the MCMC simulation, create insightful graphical representation or 
#' apply further clustering techniques to better describe the latent space. 
#' See \link[BLSM]{estimate_latent_positions} or \link[BLSM]{plot_latent_positions} for further information.
#' 
#' @import Rcpp RcppEigen 
#' @importFrom Rcpp evalCpp
#' @useDynLib BLSM
#' @references A. Donizetti, A Latent Space Model Approach for Clustering Complex Network Data, 
#' Master's Thesis, Politecnico di Milano, (2017).
#' 
#' P. D. Hoff, A. E. Raftery, M. S. Handcock, Latent Space Approaches to Social Network Analysis, 
#' Journal of the American Statistical Association, Vol. 97, No. 460, (2002), pp. 1090-1098.

NULL
