# TDA_PSEUDOTIME
## R code for temporal phenotyping via Topological Data Analysis and Pseudo Time 


### Basic Instruction for running TDA + PTS on the example file [MyData](https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MyDataSim.csv)
1. Upload and preporcessing
- The example file is in the correct format:
  - The first 3 columns idetify ids and dates 
  - Dates are in the format YYYY-mm-dd
- Create a time line for each subject starting from the first observation (minimum date)

2. Compute Distance Matrix and Function
- Choose a distance metric (Euclidean or Cosine), defualt is Cosine  
- Create Lens Functions (Mean, Max, L-infinity centrality, Single Value Decomposition 1st and 2nd compontent)
- Create Enrichment Function (Time, Single subject IDs, CRP..)

3. Run the TDA - Parameter Setting

