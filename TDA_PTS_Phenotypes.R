rm(list=ls())

setwd("/Users/ariannadagliati/Desktop/TDA/")
library(TDAmapper)
library(igraph)

library(lsa)
library(irlba)
library(dplyr)
library(philentropy)

#==========
#==========PREPROCESS
#==========
#Import Data 
FupData<-read.csv("MyDataSim.csv",header = TRUE, row.names = NULL,  colClasses="character",na.strings="NA" )

# ====== Check the format of the file
# first col > row id "id" (can be any value, i.e. unique row number)
# second col > pts id "covid_id"
# third col > date "day"

#Extract Lab Values to be used for TDA distance Matrix
LabvaluesNames<-names(FupData[,4:ncol(FupData)])


#======Create a time line for each subject (being time zero the first observation)
#Check and Change Data Format
FupData$day<-as.Date(FupData$day, format = "%Y-%m-%d")

#Find First Date
FirstDate<-FupData[,c("covid_id","day")] %>% 
  group_by(covid_id) %>% 
  slice(which.min(day))
names(FirstDate)<-c("covid_id","FirstDate")
#Merge and compute difference
FupData<-merge(FupData,FirstDate, by=c("covid_id"), all.x = TRUE)
FupData$time<-as.numeric(FupData$day - FupData$FirstDate)

#Order Data
FupData <- FupData[order(FupData$covid_id, FupData$time),]


#==========
#==========TDA SET-UP 
#==========
Labvalues <- as.data.frame(apply(FupData[,names(FupData) %in% LabvaluesNames], 2, as.numeric)) 

#Scale continuous values
Labvalues<-as.data.frame(lapply(Labvalues, as.numeric))
matrix<-scale(Labvalues)

# ========== 
# ========== Create Distance Matrix 
# ========== 

# ========== Euclidean distance
# Dist.Mtrx<-as.matrix(dist(matrix, method = "euclidean"))

# ========== Cosine distance
Dist.Mtrx<-philentropy::distance(matrix, method = "cosine")


#Rename cols and rows
colnames(Dist.Mtrx)=FupData$id #rows ids
rownames(Dist.Mtrx)=FupData$id



# ========== Create Lens FUNCTIONS (works better with EUCLIDEAN)
# Mean
f.mean<-apply(Dist.Mtrx, 1,mean)
hist(f.mean)

#L-infinity
f.max<- apply(Dist.Mtrx, 1,max)
hist(f.max)

# ========== Create Lens FUNCTIONS (works better with COSINE)

#L-infinity centrality
maxcosine<-Dist.Mtrx
maxcosine[maxcosine==1] <- 0
maxcosine<-apply(maxcosine,1,max)
LinfC<-as.data.frame(maxcosine)
LinfC<-LinfC$maxcosine
hist(LinfC)


#SVD - principal and secondary metric 
SVDval<-svd(Dist.Mtrx)
SVD1<-SVDval$u[,1]
SVD2<-SVDval$u[,2]

hist(SVD1)
hist(SVD2)


# ========== 
# ========== Create enrichment functions - time, ids and lab values..
# ========== 
f.time<-data.frame(ID=FupData$id,val=FupData$time)
hist(as.numeric(f.time$val))

f.idss<- data.frame(ID=FupData$id,PATID=FupData$covid_id)

#Create as many function as needed to enrich the topology
f.crp<-data.frame(ID=FupData$id,val=FupData$C.reactive.protein..CRP...Normal.Sensitivity.)

# ==========Create palette for enrichment (blue > green > yellow > orange > red)
colfunc<-colorRampPalette(c("#00A3DD","#60C659","#FFBC21","#FF7F1E","#EF2B2D"))



# Run TDA > grid search parameters 
INTRVLS.SEQ<-seq(5,10,1)
PRGNTG.SEQ<-seq(50,60,10)
CLUST.BINS<-c(6,8,10)


#==============
#============== Run Mapper > set parameters and lens function
#==============

ii<-2
p<-2
b<-2

F2<-mapper2D( distance_matrix= Dist.Mtrx,
              #============== choose the functions
              filter_values = list(LinfC,SVD2 ), #choose the functions
              num_intervals = c(INTRVLS.SEQ[ii],INTRVLS.SEQ[ii]),
              percent_overlap = PRGNTG.SEQ[p],
              num_bins_when_clustering = CLUST.BINS[b])

F2.graph <- graph.adjacency(F2$adjacency, mode="undirected")

#==============
#============== Adjust Colors (Enrichment Function), and Size of Nodes  
#============== Note time is set for enrichment as default - need to change the code to use another function
#==============


#===========COLORS > Enrichment (TIME ENRICHMENT)
y.mean.vertex <- data.frame()
for (i in 1:F2$num_vertices){
  points.in.vertex <- F2$points_in_vertex[[i]]
  y.mean.vertex<-rbind(y.mean.vertex,
                       data.frame( id= paste(i),
                                   #value = mean(as.numeric(f.crp[f.crp$ID %in%  rownames((Dist.Mtrx))[points.in.vertex],c("val")]))
                                   value = mean(as.numeric(f.time[f.time$ID %in%  rownames((Dist.Mtrx))[points.in.vertex],c("val")]))
                       ))  }
names(y.mean.vertex)<-c("id","value")

clrlookup<-data.frame()
clrlookup<-data.frame(clrs= colfunc(length(unique(y.mean.vertex$value))),
                      value=sort(unique(y.mean.vertex$value)))

clrMap<-merge(y.mean.vertex,clrlookup, by=c("value"), all.x=TRUE)
clrMap <- clrMap[order(clrMap$id),]
V(F2.graph)$color <- paste(clrMap$clrs)

#===========SIZE of the Vertex
vertex.size <- rep(0,F2$num_vertices)
for (i in 1:F2$num_vertices){
  points.in.vertex <- F2$points_in_vertex[[i]]
  vertex.size[i] <- length((F2$points_in_vertex[[i]]))
}

V(F2.graph)$size <- ( ((vertex.size-min(vertex.size))/(max(vertex.size)-min(vertex.size)) ) *6)+8
plot(F2.graph)




#============
#============WEIGHT EDGES - mean time in the Edges
#============
for(j in 1:length(E(F2.graph))){
  
  tail<-tail_of(F2.graph, E(F2.graph)[j])
  head<-head_of(F2.graph, E(F2.graph)[j])
  
  pointInTail<-F2$points_in_vertex[[tail]]
  pointInHead<-F2$points_in_vertex[[head]]
  
  commonIDS<-intersect(pointInTail,pointInHead)
  
  
  # ========  TIME WEIGHT
  E(F2.graph)$weight[j]<-mean( f.time[f.time$ID %in%  rownames((Dist.Mtrx))[commonIDS],c("val")])
  
  
 }
plot(F2.graph)


#======================
#======================Cluster - community detection
#======================

CommunityCluster <- edge.betweenness.community(F2.graph, weights = E(F2.graph)$value, directed = FALSE, bridges=TRUE)
plot(CommunityCluster, F2.graph)


#======================
#======================Create and plot the MST
#======================
minspantreeweights = mst(F2.graph, weights =   CommunityCluster$edge.betweenness)
plot(minspantreeweights)




#========
#========Assign Observations to a node in the network
#========
FupData$node<-NA
for (i in 1:F2$num_vertices){
  points.in.vertex <- F2$points_in_vertex[[i]]
  FupData[FupData$id %in%  rownames((Dist.Mtrx))[points.in.vertex],c("node")]<-paste(i)
}
FupData$node<-as.numeric(FupData$node)




#========
#========Find Trajectories in the MST >>> All the trajectory or choose the nodes 
#========
StartingNodes<-V(minspantreeweights)[degree(minspantreeweights)==1]
#EndingNodes<-V(minspantreeweights)[degree(minspantreeweights)>median(degree(minspantreeweights))]

#=========== OR chose specific nodes in the MST
#StartingNodes<-c(3,4,8,20,23)
EndingNodes<-c(25,18)

Trajectories<-list()
count<-1
for (i in 1:length(StartingNodes)){
  for (j in 1:length(EndingNodes)){

    Trajectories[[count]]<-all_shortest_paths(minspantreeweights,from = StartingNodes[i], to=EndingNodes[j],     mode = c("out"),                                         weights = NULL)

     count<-count+1
  }
}



#========
#========Jaccard Similarity to assign each subject to most similar trajectory
#========
simlJaccard  <- function(x,y) { (length(intersect(x,y)))/length(union(x,y))}
simlIntersection  <- function(x,y) {  length(intersect(x,y))/length(x)}
simlLength  <- function(x,y) {  abs(length(x) -  length(y))      }

uniquePts<-unique(FupData$covid_id)
Similraty<-data.frame()

for (i in 1:length(uniquePts)){
  temp<-FupData[FupData$covid_id==uniquePts[i],]
  trajOrig<-temp$node
  traj<- unique(temp$node)
  
  tempRow<-data.frame()
  
  for(t in 1:length(Trajectories)){
    if(length(Trajectories[[t]]$res)>0){
      
      temp=data.frame(covid_id=uniquePts[i],
                     trajPaz=paste(traj, collapse =" "),
                      trajNumb=t,
                      trajElmnts=paste(Trajectories[[t]]$res[[1]], collapse=" "),
                      trajLenght=length(Trajectories[[t]]$res[[1]]),
                      SJ =simlJaccard(traj,Trajectories[[t]]$res[[1]]),
                      SI= simlIntersection(traj,Trajectories[[t]]$res[[1]]),
                      SL= simlLength(traj,Trajectories[[t]]$res[[1]]))
      tempRow<-rbind(tempRow,temp)
    }
  }
  Similraty<-rbind(Similraty,tempRow)
}


MostSimilarTraj<-Similraty %>% 
  group_by(covid_id) %>% 
  slice(which.max(SJ))

#========
#========OUTPUT 
#========

MyDataOutput<-MostSimilarTraj[,c("covid_id","trajNumb")]
names(MyDataOutput)<-c("covid_id","assignedTrajectory")

table(MyDataOutput$assignedTrajectory)


#=====

