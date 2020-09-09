
setwd("/Users/user/TDA_PSEUDOTIME/")

library(TDAmapper)
library(igraph)
library(lsa)
library(irlba)
library(dplyr)
library(philentropy)

library(mice)
library(ggplot2)
library(RColorBrewer)

library(gridExtra)
library(stringdist)

# To check - not essential
library(GrpString)
library(heemod)
library(CINNA)

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

#==========
#==========Multiple Imputation (OPTIONAL - TBD)
#==========
tempData <- mice(Labvalues,m=5,maxit=50,meth='pmm',seed=500)
Labvalues <- complete(tempData,1)

#==========
#==========Scale continuous values
#==========

Labvalues<-as.data.frame(lapply(Labvalues, as.numeric))
matrix<-scale(Labvalues)

# ========== 
# ========== Create Distance Matrix (Cosine Selected) 
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
f.crp<-data.frame(ID=FupData$id,val=FupData$CRP)

# ==========Create palette for enrichment (blue > green > yellow > orange > red)
colfunc<-colorRampPalette(c("#00A3DD","#60C659","#FFBC21","#FF7F1E","#EF2B2D"))

#==============
#============== Run Mapper > set parameters and lens function
#==============

# Run TDA > grid search parameters 
INTRVLS.SEQ<-seq(5,10,1)
PRGNTG.SEQ<-seq(50,60,10)
CLUST.BINS<-c(6,8,10)


#====== Loop though different set of parameters (remember to un comment brackets after computing the MST)
# for(ii in 1:length(INTRVLS.SEQ)){
#   for (p in 1:length(PRGNTG.SEQ) ){
#     for (b in 1:length(CLUST.BINS)){

#====== Chose a fixed set
# ==== num_intervals = c(6,6),
# ==== percent_overlap = 60,
# ==== num_bins_when_clustering = 8

ii<-2
p<-2
b<-2

#====== Run the Mapper Function


F2<-mapper2D( distance_matrix= Dist.Mtrx,
              #============== choose the functions
              filter_values = list(SVD1 , SVD2 ), #choose the functions
             
               num_intervals = c(INTRVLS.SEQ[ii],INTRVLS.SEQ[ii]),
              percent_overlap = PRGNTG.SEQ[p],
              num_bins_when_clustering = CLUST.BINS[b])

F2.graph <- graph.adjacency(F2$adjacency, mode="undirected")

#==============
#============== Adjust Colors (Enrichment Function), and Size of Nodes  
#============== Note time is set for enrichment as default - need to change the code to use another function
#==============


#===========COLORS > Enrichment (here the function f.time perform a TIME ENRICHMENT)
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

CommunityCluster <- edge.betweenness.community(F2.graph, 
                                               weights = E(F2.graph)$value, 
                                               directed = FALSE, bridges=TRUE)
# Plot 
plot(CommunityCluster, F2.graph)


# Enrich the graph with  cluster colors 
V(F2.graph)$color <- paste(CommunityCluster)

# Make a palette of cluster colors
Palette<-data.frame(Color= brewer.pal(length(unique(CommunityCluster$membership)), "Set1") ,
                    Cluster=unique(CommunityCluster$membership))
NodeColor<-data.frame(Node=F2$level_of_vertex, Cluster=CommunityCluster$membership)
NodeColor<-merge(NodeColor,Palette, by=c("Cluster"))
NodeColor<-NodeColor[order(NodeColor$Node),]

# Plot with assigned palette
V(F2.graph)$color <- NodeColor$Color
plot(F2.graph)
legend("topleft",
       legend=Palette$Cluster,
       col=Palette$Color, fill=Palette$Color, 
       horiz=TRUE, box.lty=0,cex=0.8)

# Create data frame nodes - cluster
NodeClustes<-data.frame(Node=F2$level_of_vertex, Cluster=CommunityCluster$membership)




#======================
#======================Create and plot the MST
#======================
minspantreeweights = mst(F2.graph, weights =   CommunityCluster$edge.betweenness)
plot(minspantreeweights)
legend("topright",legend=Palette$Cluster,
       col=Palette$Color, fill=Palette$Color, horiz=TRUE, cex=0.8)

#}}}


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
#========Plot Values in Clusters
#========

#Plot Time
ggplot(FupData[is.na(FupData$Cluster)==FALSE,],  aes(x=as.factor(Cluster),     y=as.numeric(time),      fill=as.factor(Cluster))) + 
  geom_boxplot(alpha=0.8)+scale_fill_manual(values=Palette$Color)+ scale_color_manual(values=Palette$Color)+
theme(legend.position="none",plot.title=element_text(size=8,hjust = 0.5))


#Plot Lab Values
pltList <- lapply(LabvaluesNames, function(i){
  ggplot(FupData[is.na(FupData$Cluster)==FALSE,],aes(x=as.factor(Cluster), y=as.numeric(FupData[is.na(FupData$Cluster)==FALSE,i]),fill=as.factor(Cluster))) + 
    geom_boxplot( alpha=0.8)+ ggtitle(i)+ylab("")+xlab("")+  scale_fill_manual(values=Palette$Color)+ scale_color_manual(values=Palette$Color)+
    theme(legend.position="none",plot.title=element_text(size=8,hjust = 0.5))+
     scale_y_continuous(limits = quantile(as.numeric(na.omit(FupData[is.na(FupData$Cluster)==FALSE,i])), c(0.1, 0.9)))
  
})

grid.arrange(grobs=pltList, ncol=4)



#========
#======== Find Trajectories in the MST >>> All the trajectory or choose the nodes 
#======== THIS STEP NEED MANUAL REVIEW > check the MST plot and TIME boxplots 
#========
StartingNodes<-V(minspantreeweights)[betweenness(minspantreeweights)==min(betweenness(minspantreeweights))]
EndingNodes<-V(minspantreeweights)[betweenness(minspantreeweights)==max(betweenness(minspantreeweights))]

#=========== OR chose specific nodes in the MST
#StartingNodes<-c(3,4,8,20,23)
#EndingNodes<-c(25,18)

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
  
  # Define trajectoris among clusters (coarse granualrity)
  trajCluster<-unique(temp$Cluster)

  tempRow<-data.frame()
  
  for(t in 1:length(Trajectories)){
    if(length(Trajectories[[t]]$res)>0){
      
      temp=data.frame(covid_id=uniquePts[i],
                     trajPaz=paste(traj, collapse =" "),
                     trajPazClusters=paste(na.omit(trajCluster), collapse =" "),

                      trajNumb=t,
                      trajElmnts=paste(Trajectories[[t]]$res[[1]], collapse=" "),
                      trajLenght=length(Trajectories[[t]]$res[[1]]),
                      
                      SJ =simlJaccard(traj,Trajectories[[t]]$res[[1]]),
                      SI= simlIntersection(traj,Trajectories[[t]]$res[[1]]),
                      SL= simlLength(traj,Trajectories[[t]]$res[[1]])
                     
                      # Jaro-Winkler distance > take into account order
                      # Values are comparable to Jaccard but formally is more correct (TBD)
                      JW=  stringsim(paste(Trajectories[[t]]$res[[1]], collapse="_"),
                                     paste(traj, collapse="_"),
                                     method='jw')
                     )
      tempRow<-rbind(tempRow,temp)
    }
  }
  Similraty<-rbind(Similraty,tempRow)
}

#Chose the most similar trajectory (JW similiarity is selected)
MostSimilarTraj<-Similraty %>% 
  group_by(covid_id) %>% 
#  slice(which.max(SJ))
slice(which.max(JW))

#Manual Regrouping of Trajectories
MostSimilarTraj$trajNumbManual<-ifelse(MostSimilarTraj$trajNumb %in% c(1,2,6),"R>G","")
MostSimilarTraj$trajNumbManual<-ifelse(MostSimilarTraj$trajNumb %in% c(3,4,5,7,8),"B>G",paste(MostSimilarTraj$trajNumbManual))
MostSimilarTraj$trajNumbManual<-ifelse(MostSimilarTraj$trajNumb %in% c(9),"P>G",paste(MostSimilarTraj$trajNumbManual))



#========
#========OUTPUTs 
#========

MyDataOutput<-MostSimilarTraj[,c("covid_id","trajNumb")]
names(MyDataOutput)<-c("covid_id","assignedTrajectory")

table(MyDataOutput$assignedTrajectory)

#========Merge with data and Plot Trajectories 
FupDataTraj<-merge(FupData,MostSimilarTraj, by=c("covid_id"),all.x = TRUE)

pltList <- lapply(LabvaluesNames, function(i){
 
  ggplot(FupDataTraj, aes(time, as.numeric(FupDataTraj[,i]), colour=as.factor(trajNumbManual))) + 
    geom_smooth(aes(group=as.factor(trajNumbManual), fill=as.factor(trajNumbManual),
                    method="loess", se=FALSE) )+
    ggtitle(i)+ylab("")+xlab("")+  theme(legend.position="none",plot.title=element_text(size=8,hjust = 0.5))
  
})
grid.arrange(grobs=pltList, ncol=4)


#=====



