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
library(tidyr)

source('graph-funcs.R')
source('math-funcs.R')
source('processing.R')
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
my_colors = c("#00A3DD", "#60C659", "#FFBC21", "#FF7F1E", "#EF2B2D")
colfunc <- colorRampPalette(my_colors)


my_file = "MyDataSim.csv"

# Create palette for enrichment 


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
  mutate(
    id = as.integer(id),
    time = as.numeric(day - first_date, units = 'days')) %>% 
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

# Compute cosine similarity and SVD
# Preparing for running TDA Mapper
cosine_sim <- cosine_sim_func(lab_values_mat)



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

#====== Run the TDA Mapper Function

f_sim_map <- mapper2D(
  distance_matrix = cosine_sim,
  filter_values = list(svd1 , svd2),
  num_intervals = c(intervals_seq[ii], intervals_seq[ii]),
  percent_overlap = overlaps_seq[p],
  num_bins_when_clustering = clust_bins[b]
)

f_graph <- make_graph(f_sim_map, f_time, 'clust_color')
# color_map <- color_vertex(f_sim_map, f_time)
# V(f_graph)

#======================
#======================Create and plot the MST
#======================
minspantreeweights <- igraph::mst(f_graph, weights = f_graph$clusters$edge.betweenness)

#========
#========Assign Observations to a node in the network
#========

#' THIS STEP NEED MANUAL REVIEW > check the MST plot and TIME boxplots
out_trajectories <- find_trajectories(
  minspantreeweights, processed_data, f_sim_map, f_graph)
similarity_df <- out_trajectories[[1]]
id_node_cluster <- out_trajectories[[2]]

most_similar_traj <- similarity_df %>%
  group_by(covid_id) %>%
  # slice(which.max(JW)) %>% # use Jaro-Winkler similarity
  slice(which.max(SJ)) %>% # use Jaccard similarity for now
  mutate(trajNumbManual = case_when(
    trajNumb %in% c(1, 2, 6) ~ "R>G",
    trajNumb %in% c(3, 4, 5, 7, 8) ~ "B>G",
    trajNumb == 9 ~ "P>G",
    TRUE ~ ""
  ))
# most_similar_traj <- similarity_df %>%
#   group_by(covid_id) %>%
#   slice(which.max(JW)) %>% # use Jaccard similarity for now
#   mutate(trajNumbManual = case_when(
#     trajNumb %in% c(1, 4) ~ "R>B",
#     trajNumb %in% c(2, 3) ~ "B",
#     trajNumb %in% c(5,6 ) ~ "G>P>B",
#     trajNumb %in% c(7,8 ) ~ "P>G>P>B",
#     TRUE ~ ""
#   ))


#========
#========OUTPUTs
#========

data_out <- most_similar_traj %>% 
  select(covid_id, trajNumbManual)

table(data_out$trajNumbManual)

## Visualization

# ========== Create Lens FUNCTIONS (works better with EUCLIDEAN)
hist(apply(cosine_sim, 1, mean)) # mean
hist(apply(cosine_sim, 1, max)) # L-infinity
hist(apply(cosine_sim, 1, min))

plot_dat <- processed_data %>% 
  left_join(id_node_cluster %>% distinct(covid_id, id, cluster), 
            by = c('id', 'covid_id'))
# rowid_node_df <- convert_id_to_node(f_sim_map, processed_data)
# processed_data <- processed_data %>%
#   left_join(rowid_node_df, by = c('id', 'covid_id')) %>% 
#   left_join(f_graph$node_color, by = 'node') %>% 
#   select(covid_id, node, cluster, time) %>% 
#   distinct()


plot_dat %>% 
  ggplot(aes(x = cluster, y = time, fill = cluster)) +
  geom_boxplot(alpha = 0.8) + 
  scale_fill_manual(values = f_graph$pal$color) + 
  scale_color_manual(values = f_graph$pal$color) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5))

plot_dat %>% 
  select(cluster, lab_value_names) %>% 
  tidyr::pivot_longer(- cluster, names_to = 'Lab', values_to = 'lab_value') %>% 
  ggplot(aes(x = cluster, y = lab_value, fill = cluster)) +
  geom_boxplot(alpha = 0.8) + 
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = f_graph$pal$color) + 
  scale_color_manual(values = f_graph$pal$color) +
  theme(legend.position = "none") +
  facet_wrap(~ Lab, scales = 'free_y')


plot(minspantreeweights)
legend(
  "topright",
  legend = f_graph$pal$cluster,
  col = f_graph$pal$color,
  fill = f_graph$pal$color,
  horiz = TRUE,
  cex = 0.4
)

#========Merge with data and Plot trajectories
processed_data_traj <- processed_data %>% 
  left_join(most_similar_traj, by = c("covid_id")) %>% 
  mutate(trajNumbManual = as.factor(trajNumbManual), time) %>%
  select(time, trajNumbManual, lab_value_names) %>% 
  distinct()

processed_data_traj  %>% 
  tidyr::pivot_longer(- c(time, trajNumbManual), 
                      names_to = 'Lab', values_to = 'lab_value') %>% 
  ggplot(aes(time, lab_value, colour = trajNumbManual, 
             group = trajNumbManual, fill = trajNumbManual)) +
  geom_smooth(method = "loess") +
  theme(legend.position = "none") +
  facet_wrap(~ Lab, scales = 'free_y')


#=====
