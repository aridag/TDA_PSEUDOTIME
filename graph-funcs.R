#' Make igraph object by calculating adjacency, resizing nodes,
#' weighting edges and color nodes.
#'
#' @param f_sim_map TDAmapper object
#' @param f_time Dataframe of id and enrichment variable. 
#' Default is time (performing time enrichment).
#' However, f_time can contain any other variable.
#'
#' @return igraph object of the graph output.
#' @export
#'
#' @examples
make_graph <- function(f_sim_map, f_time, color_method = 'clust_color') {
  f_graph <- igraph::graph.adjacency(f_sim_map$adjacency, mode = "undirected")
  f_graph <- resize_nodes(f_sim_map, f_graph)
  f_graph <- weight_edges(f_sim_map, f_graph, f_time)
  f_graph <- color_graph(f_sim_map, f_graph, f_time, color_method)
  f_graph
}


#' Color nodes of graph given coloring method.
#'
#' @param f_sim_map TDAmapper object
#' @param f_graph igraph object, out put from graph.adjacency
#' @param f_time Data frame of the original data but 
#' with only two columns: ID and val (time)
#' @param method Character string specifying the coloring method.
#' Can be 'basic', 'clust_shade' or 'clust_color'.
#'
#' @return Modified graph with colors at nodes.
#' @export
#'
#' @examples
color_graph <- function(
  f_sim_map, f_graph, f_time, 
  method = c('basic', 'clust_shade', 'clust_color')) {
  
  method <- match.arg(method)
  
  if (method == 'basic'){
    color_map <- color_vertex(f_sim_map, f_time)
    V(f_graph)$color <- color_map$colors
    plot(f_graph)
  } else if (method == 'clust_shade'){
    my_clusters <- commu_clus(f_graph) # community detection
    plot(my_clusters, f_graph) # highlight community clusters
    # enrich graph with cluster colors
    V(f_graph)$color <- paste(my_clusters$membership)
    f_graph$clusters <- my_clusters
  } else {
    my_clusters <- commu_clus(f_graph) # community detection
    node_color <- color_clust(f_sim_map, my_clusters)
    V(f_graph)$color <- node_color$color
    pal <- node_color %>% select(- node) %>% distinct()
    plot(f_graph) # plot with assigned palette
    legend(
      "topleft",
      legend = pal$cluster,
      col = pal$color,
      fill = pal$color,
      horiz = TRUE,
      box.lty = 0,
      cex = 0.8
    )
    f_graph$node_color <- node_color
    f_graph$clusters <- my_clusters
    f_graph$pal <- pal
  }
  f_graph
}



#' Resize the nodes of graph given the number of points it contains
#'
#' @param f_sim_map TDAmapper object
#' @param f_graph igraph object, out put from graph.adjacency
#' @param mult Scalar as a multiplicative size rescale parameter.
#' @param add Scalar as an additive size rescale parameter.
#'
#' @return An igraph object of the modified f_graph.
#' @export
#'
#' @examples
resize_nodes <- function(f_sim_map, f_graph, mult = 6, add = 8){
  n_vertices <- f_sim_map$num_vertices
  vertex.size <- vector('numeric', n_vertices)
  for (i in seq.int(n_vertices)) {
    points.in.vertex <- f_sim_map$points_in_vertex[[i]]
    vertex.size[i] <- length(f_sim_map$points_in_vertex[[i]])
  }
  min_size <- min(vertex.size)
  max_size <- max(vertex.size)
  V(f_graph)$size <- (vertex.size - min_size) / (max_size - min_size) * mult + add
  f_graph
}

#' Weight edges of a graph based on mean time of each edge.
#'
#' @param f_sim_map TDAmapper object
#' @param f_graph igraph object, out put from graph.adjacency
#' @param f_time Data frame of the original data but 
#' with only two columns: ID and val (time)
#'
#' @return An igraph object of the modified f_graph.
#' @export
#'
#' @examples
weight_edges <- function(f_sim_map, f_graph, f_time){
  n_edges <- length(E(f_graph))
  
  for (j in seq.int(n_edges)) {
    tail <- igraph::tail_of(f_graph, E(f_graph)[j])
    head <- igraph::head_of(f_graph, E(f_graph)[j])
    pointInTail <- f_sim_map$points_in_vertex[[tail]]
    pointInHead <- f_sim_map$points_in_vertex[[head]]
    commonIDS <- intersect(pointInTail, pointInHead)
    
    # ========  TIME WEIGHT
    E(f_graph)$weight[j] <- f_time %>% 
      filter(ID %in% f_time$ID[commonIDS]) %>% 
      pull(val) %>% 
      mean()
  }
  f_graph
}

#' Adjust colors of nodes using enrichment function.
#'
#' @param f_sim_map TDAmapper object
#' @param f_time Data frame of the original data but 
#' with only two columns: ID and val (time)
#'
#' @return Data frame of vertex ids and colors.
#' @export
#'
#' @examples
color_vertex <- function(f_sim_map, f_time){
  y.mean.vertex <- list()
  for (i in 1:f_sim_map$num_vertices) {
    points.in.vertex <- f_sim_map$points_in_vertex[[i]]
    y.mean.vertex[[i]] <- 
      data.frame(id = paste(i),
                 value = f_time %>% 
                   filter(ID %in% f_time$ID[points.in.vertex]) %>% 
                   pull(val) %>% 
                   as.numeric() %>% 
                   mean())
  }
  
  color_map <- bind_rows(y.mean.vertex) %>% 
    arrange(value) %>% 
    mutate(colors = unique(.$value) %>% length() %>% colfunc()) %>% 
    # arrange(as.numeric(id)) # should it be this instead???
    arrange(id)
}

#' Detect community structure based on edge betweeness.
#'
#' @param f_graph igraph object, out put from graph.adjacency
#' @param directed Logical. whether to calculate directed edge betweenness 
#' for directed graphs. Ignored for undirected graphs.
#' @param bridges Logical. whether to return a list the edge removals 
#' which actually a component of the graph.
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
commu_clus <- function(f_graph, directed = FALSE, bridges = TRUE, ...){
  igraph::edge.betweenness.community(
    f_graph,
    weights = E(f_graph)$value,
    directed = directed,
    bridges = bridges,
    ...
  )
}


#' Get cluster colors for nodes.
#'
#' @param f_sim_map TDAmapper object
#' @param my_clusters 
#' 
#' @return Data frame of nodes and corresponding colors
#' based on the cluster each node belongs.
#' @export
#'
#' @examples
color_clust <- function(f_sim_map, my_clusters) {
  cluster_vec <- as.factor(unique(my_clusters$membership))
  # Make a palette of cluster colors
  my_palette <- data.frame(
    color = RColorBrewer::brewer.pal(length(cluster_vec), "Set1"),
    cluster = cluster_vec
  )
  
  # Create data frame of nodes and cluster
  node_color <- data.frame(
    node = f_sim_map$level_of_vertex,
    cluster = as.factor(my_clusters$membership)
  ) %>% 
    left_join(my_palette, by = 'cluster') %>% 
    arrange(node)
}


