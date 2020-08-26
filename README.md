# TDA_PSEUDOTIME
## R code for temporal phenotyping via Topological Data Analysis and Pseudo Time 

### Instruction for running TDA + PTS on the example file [MyData](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MyDataSim.csv)

1. Upload and preporcessing
- The example file is in the correct format:
  - The first 3 columns idetify ids and dates 
  - Dates are in the format YYYY-mm-dd
- Create a time line for each subject starting from the first observation (minimum date)

2. Compute Distance Matrix and Function
- Choose a distance metric (Euclidean or Cosine), defualt is Cosine  
- Create Lens Functions (Mean, Max, L-infinity centrality, Single Value Decomposition 1st and 2nd compontent)
- Create Enrichment Function (Time, Single subject IDs, CRP..)
- Create a Color palette (list of colors)

3. Run the TDA - Parameter Setting
Note that the code is also menat to be used for a grid search (testing differnt combination of parameters)
- Create sequences of parameters:
  - INTRVLS.SEQ number of bins (index is ii)
  - PRGNTG.SEQ is how much the bins overlap (index is p)
  - CLUST.BINS indicate how to cut the hierarchical clustering dendogram of (index is b) 
  
 - Run the 2d Mapper (two dimensions)
  - Use 2 selected functions (defualt are L-infinity centrality + Single Value Decomposition 1st component) 
 
4. Enrich the topology


![Topological Map](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/TopologicalMap.png)

5. Create the Minimum Spanning Tree and retrive Trajectories

![Minimum Spanning Tree](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MSTExample.png)

6. Assign each subject to the most similar trajectory
