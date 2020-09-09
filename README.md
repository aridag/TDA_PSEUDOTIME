# TDA_PSEUDOTIME
## R code for temporal phenotyping via Topological Data Analysis and Pseudo Time 

### Instruction for running TDA + PTS on the example file [MyData](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MyDataSim.csv)
*Note that data are randomly created on the basis of fake ids, fake dates and within real Lab value ranges, thus they don't have any clinical meaning*

1. Upload and pre-processing
- The example file is in the correct format:
  - The first 3 columns identify ids and dates 
  - Dates are in the format YYYY-mm-dd
- Create a time line for each subject starting from the first observation (minimum date)

- Imputation - this is optional (e.g. the example file has no missingness)

2. Compute Distance Matrix and Function
- Choose a distance metric (Euclidean or Cosine), default is Cosine  
- Create Lens Functions (Mean, Max, L-infinity centrality, Single Value Decomposition 1st and 2nd component)
- Create Enrichment Function (Time, Single subject IDs, CRP..)
- Create a Colour palette (list of colours)

3. Run the TDA - Parameter Setting
Note that the code is also meant to be used for a grid search (testing different combination of parameters)
- Create sequences of parameters:
  - INTRVLS.SEQ number of bins (index is ii)
  - PRGNTG.SEQ is how much the bins overlap (index is p)
  - CLUST.BINS indicate how to cut the hierarchical clustering dendrogram of (index is b) 
  
 - Run the 2d Mapper (two dimensions)
  - Use 2 selected functions (default are Single Value Decomposition 1st and 2nd components) 


4. Enrich the topology
- Enrich nodes with functions (default is f.time, the mean time in each node). 

The results of these first steps is a topological map as this one

![Topological Map](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/TopologicalMap.png)


5. Weight the edges of the network with mean time (mean time of the observation in the edges)


6. Find Clusters
- Apply community detection clustering to the topology using the weighted edges
- Enrich the topology with clustering infromation (colours)


7. Create the Minimum Spanning Tree (MST) and retrieve Trajectories
- Create clusters based on the weighted edges
- Create and plot the Minimum Spanning Tree

The results of is a graph like this one

![Minimum Spanning Tree](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MSTExample.png)



8. Find Trajectories in the MST
There are two possibilities:
-	Compute all the trajectories starting form starting nodes (i.e. with degree equal to one) to ending nodes (in a very heuristic way I have defined them as those with a degree > median degree)
-	Choose the nodes, in this example it was clear that two good Ending Nodes candidates are node 18 and 25

9. Compute Jaccard similarity and assign subjects to the most similar trajectory
- Compute the Jaccard similarity between the “real trajectory” of the subject and all the mined trajectories, assign the subject to the most similar one

10. The output is a data frame with two columns: subject id (covid_id) and assigned trajectory.



