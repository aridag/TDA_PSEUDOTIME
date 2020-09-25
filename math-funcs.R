#' Calculate the pairwise similarity between rows of input matrix
#'
#' @param input_mat Input matrix.
#'
#' @return
#' @export Square matrix of the similarity. Diagonal values should be 1.
#'
#' @examples
#' m <- matrix(1:8, ncol = 2) # a 4x2 matrix
#' cosine_sim_func(m) # a 4x4 similarity matrix
cosine_sim_func <- function(input_mat){
  cosine_sim <- input_mat / sqrt(rowSums(input_mat * input_mat))
  cosine_sim <- cosine_sim %*% t(cosine_sim)
  cosine_sim[cosine_sim > 1] <- 1.0 
  cosine_sim
}

# ==============================================================
# SIMILARITY FUNCTIONS

#' Compute Jaccard similarity to assign each subject to 
#' most similar trajectory
#'
#' @param x Numeric vector of node numbers
#' indicating an individual's trajectory
#' @param y Numeric vector of node numbers
#' indicating a general trajectory
#'
#' @return A scalar as measure of similarity
#' @export
#'
#' @examples 
#' sim_jaccard(c(1,2), c(1,3,4))
sim_jaccard  <- function(x, y) {
  (length(intersect(x, y))) / length(union(x, y))
}

#' Compute intersection similarity to assign each subject to 
#' most similar trajectory
#'
#' @param x Numeric vector of node numbers
#' indicating an individual's trajectory
#' @param y Numeric vector of node numbers
#' indicating a general trajectory
#'
#' @return A scalar as measure of similarity
#' @export
#'
#' @examples 
#' sim_intersection(c(1,2), c(1,3,4))
sim_intersection  <- function(x, y) {
  length(intersect(x, y)) / length(x)
}

#' Compute similarity based on length to assign each subject to 
#' most similar trajectory
#'
#' @param x Numeric vector of node numbers
#' indicating an individual's trajectory
#' @param y Numeric vector of node numbers
#' indicating a general trajectory
#'
#' @return A scalar as measure of similarity
#' @export
#'
#' @examples 
#' sim_length(c(1,2), c(1,3,4))
sim_length  <- function(x, y) {
  abs(length(x) - length(y)) # should we do exp(-...)?
}


#' Find the shortest paths.
#'
#' @param x Pair of node ids from start node to end node.
#'
#' @return
#' @export
#'
#' @examples
shortest_paths_func <- function(x, mst_weights){
  all_shortest_paths(
    mst_weights,
    from = x[1],
    to = x[2],
    mode = c("out"),
    weights = NULL
  )
}
