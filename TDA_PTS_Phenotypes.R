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

# brewer.pal(length(unique(CommunityCluster$membership)), "Set1")
my_colors = c("#00A3DD", "#60C659", "#FFBC21", "#FF7F1E", "#Ef_sim_mapB2D")
my_file = "MyDataSim.csv"

# Create palette for enrichment 
colfunc <- colorRampPalette(my_colors)


# Run TDA > grid search parameters
# hyperparameters
intervals_seq <- seq(5, 10, 1)
overlaps_seq <- seq(50, 60, 10)
clust_bins <- c(6, 8, 10)
ii <- 2
p <- 2
b <- 2
num_intervals = c(intervals_seq[ii], intervals_seq[ii])
percent_overlap = overlaps_seq[p]
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
  mutate(time = as.numeric(day - first_date, units = 'days')) %>% 
  arrange(covid_id, time) %>% 
  mutate_at(vars(lab_value_names), as.numeric)

lab_values_mat <- processed_data[, lab_value_names] %>% 
  mice(m = 5, maxit = 50, meth = 'pmm', seed = 500) %>% # imputation
  complete(1) %>% 
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

#==========
#==========TDA SET-UP
#==========

# Matrix <- as.matrix(matrix)

cosine_sim <- cosine_sim_func(lab_values_mat)

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
# enrich topology by any variable
# use 'time' for now.
f_time <- processed_data %>% 
  select(ID = id, val = time) # or val = CRP, etc.
hist(as.numeric(f_time$val))


#==============
#============== Run Mapper > set parameters and lens function
#==============

#====== Run the Mapper Function


f_sim_map <- mapper2D(
  distance_matrix = cosine_sim,
  filter_values = list(svd1 , svd2),
  num_intervals = c(intervals_seq[ii], intervals_seq[ii]),
  percent_overlap = overlaps_seq[p],
  num_bins_when_clustering = clust_bins[b]
)

f_graph <- make_graph(f_sim_map, f_time, 'clust_color')


#======================
#======================Create and plot the MST
#======================
minspantreeweights <- igraph::mst(f_graph, weights = f_graph$clusters$edge.betweenness)
plot(minspantreeweights)
legend(
  "topright",
  legend = f_graph$pal$cluster,
  col = f_graph$pal$color,
  fill = f_graph$pal$color,
  horiz = TRUE,
  cex = 0.4
)

#========
#========Assign Observations to a node in the network
#========
# processed_data$node <- NA
### ASK ARIANNA HERE!!!!!!!!!!!!!!!!!!!!!!
# what if a patient belongs to multiple nodes??????
n_vertices <- f_sim_map$num_vertices
for (i in seq.int(n_vertices)) {
  points.in.vertex <- f_sim_map$points_in_vertex[[i]]
  processed_data$node[processed_data$id %in% points.in.vertex] <- i
}
processed_data <- processed_data %>% 
  left_join(f_graph$node_color, by = 'node')

#========
#========Plot Values in clusters
#========

#Plot Time
processed_data %>% 
  ggplot(aes(x = cluster, y = time, fill = cluster)) +
  geom_boxplot(alpha = 0.8) + 
  scale_fill_manual(values = Palette$color) + 
  scale_color_manual(values = Palette$color) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5))

processed_data %>% 
  mutate(cluster = cluster) %>% 
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



trajectories <- apply(starts_and_ends, 1, shortest_paths_func, mst_weights = minspantreeweights)

# trajectories[[3]]$res

unique_patients <- unique(processed_data$covid_id)

similarity_ls <- list()

for (i in unique_patients) {
  patient_i <- processed_data %>% filter(covid_id == i)
  traj <- unique(patient_i$node)
  
  # Define trajectoris among clusters (coarse granularity)
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
        SJ = sim_jaccard(traj, traj_res_1),
        SI = sim_intersection(traj, traj_res_1),
        SL = sim_length(traj, traj_res_1),
        # Jaro-Winkler distance > take into account order
        # Values are comparable to Jaccard but formally is more correct (TBD)
        JW =  stringsim(
          paste(traj_res_1, collapse = "_"),
          paste(traj, collapse = "_"),
          method = 'jw'
        )
      )
    }
  }
  similarity_ls[[i]] <- bind_rows(temp)
}

similarity_df <- bind_rows(similarity_ls)

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
         time = time) %>%
  select(time, clusterTraj, lab_value_names) %>% 
  tidyr::pivot_longer(- c(time, clusterTraj), 
                      names_to = 'Lab', values_to = 'lab_value') %>% 
  ggplot(aes(time, lab_value, colour = clusterTraj, 
             group = clusterTraj, fill = clusterTraj)) +
  geom_smooth(method = "loess") +
  theme(legend.position = "none") +
  facet_wrap(~ Lab, scales = 'free_y')


#=====
