

# setwd("/Users/user/TDA_PSEUDOTIME/")

library(TDAmapper)
library(igraph)
# library(lsa)
# library(irlba)
library(dplyr)

library(mice)
library(ggplot2)
library(RColorBrewer)

library(gridExtra)
library(stringdist)

# To check - not essential
# library(GrpString)
# library(heemod)
# library(corpcor) #fast svd
# library(CINNA)

#==========
#==========PREPROCESS
#==========
#Import Data
# ========== (blue > green > yellow > orange > red)
set.seed(1234)

my_colors = c("#00A3DD", "#60C659", "#FFBC21", "#FF7F1E", "#EF2B2D")
my_file = "MyDataSim.csv"

# Create palette for enrichment 
colfunc <- colorRampPalette(my_colors)


# Run TDA > grid search parameters
intrvls_seq <- seq(5, 10, 1)
prgntg_seq <- seq(50, 60, 10)
clust_bins <- c(6, 8, 10)
ii <- 2
p <- 2
b <- 2
num_intervals = c(intrvls_seq[ii], intrvls_seq[ii])
percent_overlap = prgntg_seq[p]
num_bins_when_clustering = clust_bins[b]

FupData <- read.csv(my_file, header = TRUE, colClasses = "character")

# ====== Check the format of the file
# first col > row id "id" (can be any value, i.e. unique row number)
# second col > pts id "covid_id"
# third col > date "day"

#Extract Lab Values to be used for TDA distance Matrix
non_lab_value_names <- c('id', 'covid_id', 'day')
lab_value_names <- setdiff(names(FupData), non_lab_value_names)

#======Create a time line for each subject (being time zero the first observation)
#Check and Change Data Format

processed_data <- FupData %>% 
  mutate(day = as.Date(day, format = '%Y-%m-%d')) %>% 
  group_by(covid_id) %>% 
  mutate(first_date = min(day)) %>% 
  ungroup() %>% 
  mutate(time = day - first_date) %>% 
  arrange(covid_id, time) %>% 
  mutate_at(vars(lab_value_names), as.numeric)

lab_values_df <- processed_data[, lab_value_names] %>% 
  mice(m = 5, maxit = 50, meth = 'pmm', seed = 500) %>% # imputation
  complete(1)
lab_values_mat <- lab_values_df  %>% 
  scale() # prepare for cosine similarity calculation


# check:
# > sum(processed_data$first_date != FupData$FirstDate)
# [1] 0
# > sum(processed_data$time != FupData$time)
# [1] 0
# > sum(diff(as.numeric(processed_data$id))!=1)
# [1] 0
# > sum(diff(as.numeric(FupData$id))!=1)
# [1] 0
# str(Labvalues)
# str(lab_values_df)

#==========
#==========TDA SET-UP
#==========

# Matrix <- as.matrix(matrix)
cosine_sim_func <- function(input_mat){
  cosine_sim <- input_mat / sqrt(rowSums(input_mat * input_mat))
  cosine_sim <- cosine_sim %*% t(cosine_sim)
  cosine_sim[cosine_sim > 1] <- 1.0 
  cosine_sim
}
cosine_sim <- cosine_sim_func(lab_values_mat)
#Rename cols and rows
colnames(cosine_sim) <- 
  rownames(cosine_sim) <- 
  cosine_id <- 
  processed_data$id #rows ids

# ========== Create Lens FUNCTIONS (works better with EUCLIDEAN)
hist(apply(cosine_sim, 1, mean)) # mean
hist(apply(cosine_sim, 1, max)) # L-infinity
hist(apply(cosine_sim, 1, min))

# ========== Create Lens FUNCTIONS (works better with COSINE)

# SVD - principal and secondary metric

trunc_svds <- RSpectra::svds(cosine_sim, k = 2, nu = 2, nv = 2)
svd1 <- - trunc_svds$u[, 1]
svd2 <- - trunc_svds$u[, 2]

# takes about 70s to run
# system.time(svd_list <- svd(cosine_sim))
# SVD1 <- svd_list$u[, 1]
# SVD2 <- svd_list$u[, 2]
# hist(SVD1)
# hist(SVD2)
# tolerance = 1^(-10)
# sum(svd2 - SVD2 > tolerance)
# sum(svd1 - SVD1 > tolerance)

# ==========
# ========== Create enrichment functions - time, ids and lab values..
# ==========
f.time <- processed_data %>% 
  select(ID = id, val = time)
hist(as.numeric(f.time$val))

# f.idss <- processed_data %>% 
#   select(id, patient_id = covid_id)

#Create as many function as needed to enrich the topology
f.crp <- processed_data %>% 
  select(ID = id, val = CRP)


#====== Loop though different set of parameters (remember to un comment brackets after computing the MST)
# for(ii in 1:length(intrvls_seq)){
#   for (p in 1:length(prgntg_seq) ){
#     for (b in 1:length(clust_bins)){

#====== Chose a fixed set
# ==== num_intervals = c(6,6),
# ==== percent_overlap = 60,
# ==== num_bins_when_clustering = 8

#==============
#============== Run Mapper > set parameters and lens function
#==============

#====== Run the Mapper Function


F2 <- mapper2D(
  distance_matrix = cosine_sim,
  #============== choose the functions
  filter_values = list(svd1 , svd2),
  #choose the functions
  num_intervals = c(intrvls_seq[ii], intrvls_seq[ii]),
  percent_overlap = prgntg_seq[p],
  num_bins_when_clustering = clust_bins[b]
)

F2.graph <- graph.adjacency(F2$adjacency, mode = "undirected")

#==============
#============== Adjust Colors (Enrichment Function), and Size of Nodes
#============== Note time is set for enrichment as default - need to change the code to use another function
#==============


#===========COLORS > Enrichment (here the function f.time perform a TIME ENRICHMENT)
y.mean.vertex <- list()
for (i in 1:F2$num_vertices) {
  points.in.vertex <- F2$points_in_vertex[[i]]
  y.mean.vertex[[i]] <- 
    data.frame(
      id = paste(i),
      value = f.time %>% 
        filter(ID %in% cosine_id[points.in.vertex]) %>% 
        pull(val) %>% 
    as.numeric() %>% 
    mean())
}

clrMap <- bind_rows(y.mean.vertex) %>% 
  arrange(value) %>% 
  mutate(clrs = unique(.$value) %>% length() %>% colfunc()) %>% 
  # arrange(as.numeric(id)) # should it be this instead???
  arrange(id)

V(F2.graph)$color <- clrMap$clrs
plot(F2.graph)

#===========SIZE of the Vertex

resize_nodes <- function(F2, F2.graph, mult = 6, add = 8){
  n_vertices <- F2$num_vertices
  vertex.size <- vector('numeric', n_vertices)
  for (i in seq.int(n_vertices)) {
    points.in.vertex <- F2$points_in_vertex[[i]]
    vertex.size[i] <- length((F2$points_in_vertex[[i]]))
  }
  min_size <- min(vertex.size)
  max_size <- max(vertex.size)
  V(F2.graph)$size <- (vertex.size - min_size) / (max_size - min_size) * mult + add
  F2.graph
}
F2.graph <- resize_nodes(F2, F2.graph)
# plot(F2.graph)

#============
#============WEIGHT EDGES - mean time in the Edges
#============
n_edges <- length(E(F2.graph))

for (j in seq.int(n_edges)) {
  tail <- tail_of(F2.graph, E(F2.graph)[j])
  head <- head_of(F2.graph, E(F2.graph)[j])
  pointInTail <- F2$points_in_vertex[[tail]]
  pointInHead <- F2$points_in_vertex[[head]]
  
  commonIDS <- intersect(pointInTail, pointInHead)
  
  # ========  TIME WEIGHT
  E(F2.graph)$weight[j] <- f.time %>% 
    filter(ID %in% cosine_id[commonIDS]) %>% 
    pull(val) %>% 
    mean()
}
plot(F2.graph)


#======================
#======================cluster - community detection
#======================

# highlight community clusters
commu_clus <- function(F2.graph, directed = FALSE, bridges = TRUE, ...){
  edge.betweenness.community(
    F2.graph,
    weights = E(F2.graph)$value,
    directed = directed,
    bridges = bridges,
    ...
  )
}
my_clusters <- commu_clus(F2.graph)
plot(my_clusters, F2.graph)
# dev.off()

# Enrich the graph with  cluster colors
V(F2.graph)$color <- paste(my_clusters$membership)

cluster_vec <- as.factor(unique(my_clusters$membership))
# Make a palette of cluster colors
Palette <-
  data.frame(
    color = brewer.pal(length(cluster_vec), "Set1"),
    cluster = cluster_vec
  )

# Create data frame nodes - cluster

node_color <-
  data.frame(node = F2$level_of_vertex,
             cluster = as.factor(my_clusters$membership)) %>% 
  left_join(Palette, by = 'cluster') %>% 
  arrange(node)


# Plot with assigned palette
V(F2.graph)$color <- node_color$color
plot(F2.graph)
legend(
  "topleft",
  legend = Palette$cluster,
  col = Palette$color,
  fill = Palette$color,
  horiz = TRUE,
  box.lty = 0,
  cex = 0.8
)



#======================
#======================Create and plot the MST
#======================
minspantreeweights <- mst(F2.graph, weights =   my_clusters$edge.betweenness)
plot(minspantreeweights)
legend(
  "topright",
  legend = Palette$cluster,
  col = Palette$color,
  fill = Palette$color,
  horiz = TRUE,
  cex = 0.4
)

#}}}


#========
#========Assign Observations to a node in the network
#========
# processed_data$node <- NA
n_vertices <- F2$num_vertices
for (i in seq.int(n_vertices)) {
  points.in.vertex <- F2$points_in_vertex[[i]]
  processed_data[processed_data$id %in% cosine_id[points.in.vertex], "node"] <-
    i
}
# processed_data$node <- as.numeric(processed_data$node)

# processed_data_ori = processed_data

processed_data <- processed_data%>% 
  left_join(node_color, by = 'node')

#========
#========Plot Values in clusters
#========

#Plot Time
processed_data %>% 
  ggplot(aes(x = as.factor(cluster), y = as.numeric(time), fill = as.factor(cluster))) +
  geom_boxplot(alpha = 0.8) + 
  scale_fill_manual(values = Palette$color) + 
  scale_color_manual(values = Palette$color) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5))

processed_data %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  select(cluster, lab_value_names) %>% 
  tidyr::pivot_longer(- cluster, names_to = 'Lab', values_to = 'lab_value') %>% 
  ggplot(aes(x = cluster, y = lab_value, fill = cluster)) +
  geom_boxplot(alpha = 0.8) + 
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = Palette$color) + 
  scale_color_manual(values = Palette$color) +
  theme(legend.position = "none") +
  facet_wrap(~ Lab, scales = 'free_y')


#========
#======== Find trajectories in the MST >>> All the trajectory or choose the nodes
#======== THIS STEP NEED MANUAL REVIEW > check the MST plot and TIME boxplots
#========
mst_between <- betweenness(minspantreeweights)
starting_nodes <- V(minspantreeweights)[mst_between == min(mst_between)]
ending_nodes <- V(minspantreeweights)[mst_between == max(mst_between)]
starts_and_ends <- expand.grid(starting_nodes, ending_nodes)

#=========== OR chose specific nodes in the MST
#starting_nodes<-c(3,4,8,20,23)
#ending_nodes<-c(25,18)

shortest_paths_func <- function(x){
  # x: pair from start node to end node
  all_shortest_paths(
  minspantreeweights,
  from = x[1],
  to = x[2],
  mode = c("out"),
  weights = NULL
)
}



#========
#========Jaccard Similarity to assign each subject to most similar trajectory
#========
simlJaccard  <- function(x, y) {
    (length(intersect(x, y))) / length(union(x, y))
  }
simlIntersection  <- function(x, y) {
    length(intersect(x, y)) / length(x)
  }
simlLength  <- function(x, y) {
  abs(length(x) - length(y))
}

trajectories <- apply(starts_and_ends, 1, shortest_paths_func)

# trajectories[[3]]$res

unique_patients <- unique(processed_data$covid_id)

similarity_df <- data.frame()

for (i in unique_patients) {
  patient_i <- processed_data %>% filter(covid_id == i)
  traj <- unique(patient_i$node)
  
  # Define trajectoris among clusters (coarse granualrity)
  trajcluster <- unique(patient_i$cluster)
  
  temp <- list()
  
  for (t in 1:length(trajectories)) {
    if (length(trajectories[[t]]$res) > 0) {
      traj_res_1 <- trajectories[[t]]$res[[1]]
      temp[[t]] <- data.frame(
        covid_id = i,
        trajPaz = paste(traj, collapse = " "),
        trajPazclusters = paste(na.omit(trajcluster), collapse = " "),
        trajNumb = t,
        trajElmnts = paste(traj_res_1, collapse = " "),
        trajLenght = length(traj_res_1),
        SJ = simlJaccard(traj, traj_res_1),
        SI = simlIntersection(traj, traj_res_1),
        SL = simlLength(traj, traj_res_1),
        # Jaro-Winkler distance > take into account order
        # Values are comparable to Jaccard but formally is more correct (TBD)
        JW =  stringsim(
          paste(traj_res_1, collapse = "_"),
          paste(traj, collapse = "_"),
          method = 'jw'
        )
      )
    }
    tempRow <- bind_rows(temp)
    
  }
  similarity_df <- rbind(similarity_df, tempRow)
}

#Chose the most similar trajectory (JW similiarity is selected)
MostSimilarTraj <- similarity_df %>%
  group_by(covid_id) %>%
  #  slice(which.max(SJ))
  slice(which.max(JW))

# TO REPLACE
#Manual Regrouping of trajectories
# MostSimilarTraj$trajNumbManual <-
#   ifelse(MostSimilarTraj$trajNumb %in% c(1, 2, 6), "R>G", "")
# MostSimilarTraj$trajNumbManual <-
#   ifelse(
#     MostSimilarTraj$trajNumb %in% c(3, 4, 5, 7, 8),
#     "B>G",
#     paste(MostSimilarTraj$trajNumbManual)
#   )
# MostSimilarTraj$trajNumbManual <-
#   ifelse(MostSimilarTraj$trajNumb %in% c(9),
#          "P>G",
#          paste(MostSimilarTraj$trajNumbManual))

#Add info about clusters
traj_cluster <- list()

for(t in 1:length(trajectories)) {
  
  if (length(trajectories[[t]]$res) > 0) {
    traj_res_1 <- trajectories[[t]]$res[[1]]
    temp <- data.frame(node = as.vector(traj_res_1)) %>% 
      left_join(node_color, by = 'node')
    
    traj_cluster[[t]] <- data.frame(
      trajElmnts = paste(traj_res_1, collapse = " "),
      clusterTraj = paste(unique(temp$cluster), collapse = ">")
    )
  }
}
traj_clusters <- bind_rows(traj_cluster)

#Add cluster information into trajectories
similarity_df <- left_join(similarity_df, traj_clusters, by=c("trajElmnts"))


#========
#========OUTPUTs
#========

MyDataOutput <- similarity_df[, c("covid_id", "clusterTraj")]
names(MyDataOutput) <- c("covid_id", "assignedTrajectory")

table(MyDataOutput$assignedTrajectory)

#========Merge with data and Plot trajectories
processed_data_traj <-
  left_join(processed_data, similarity_df, by = c("covid_id"))

processed_data_traj %>% 
  mutate(clusterTraj = as.factor(clusterTraj),
         time = as.numeric(time)) %>%
  select(time, clusterTraj, lab_value_names) %>% 
  tidyr::pivot_longer(- c(time, clusterTraj), 
                      names_to = 'Lab', values_to = 'lab_value') %>% 
  ggplot(aes(time, lab_value, colour = clusterTraj, 
             group = clusterTraj, fill = clusterTraj)) +
  geom_smooth(method = "loess") +
  theme(legend.position = "none") +
  facet_wrap(~ Lab, scales = 'free_y')


#=====
