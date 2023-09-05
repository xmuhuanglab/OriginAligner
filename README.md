# SIMarker
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/new_pipeline.png)
## Introduction:
Here we designed a modified statistical pipeline SIMarker based on the hypothesis that similar characteristics in different samples might 
have some clinical applications.

## Pull the docker (the env about R version)
```
docker pull sluo112211/sim-docker
```

## Workflow:
### Part1: Discovering potential Cirrhosis-like subpopulations
![Image text](https://github.com/xmuhuanglab/SIMarker/blob/main/images/OriginAligner.PNG)
##### Step 1 ： Calculate the hvg based on dispersion and weight（function--Screen_hvg）
##### Step 2 ： Select the hvg by yourself （function--Plot_hvg） 
##### Step 3 ： Calculate KNN and visualize  (function--OriginAligner), you had better choice the best result to visualize
##### Step 4 ： Calculate the overlap of markers and the Jaccard similarity（function--Venn_Jaccard）

### Part2: Screening and intepreting potential cirrhosis-like signatues
![Image text](https://github.com/xmuhuanglab/SIMarker/blob/main/images/part2.PNG)
##### Step 1 ： Unsupervised matrix factorization for extracting cirrhosis-like signatures （function--cNMF）
##### Step 2 ： Explainability of cirrhosis-like signatures

### Part3: Clinical applications of Cirrhosis-like signatures
![Image text](https://github.com/xmuhuanglab/SIMarker/blob/main/images/part3.PNG)
##### Step 1 ： Reversal matrix of gene pairs （function--k-TSP）
##### Step 2 ： Select of best gene pair number
##### Step 3 ： XGBoost: Training the best predictive model （function--XGBoost）
##### Step 4 ： Evaluate the capability of predictive model

## Example for code
```
source('SIMarker/Part1: Discovering potential Cirrhosis-like subpopulations/OriginAligner_main')
Screen_hvg(seurat_obj,subtype,gene_number)
Plot_hvg(disp_threshold,weight_threshold)
```
```
OriginAligner(seurat_obj,subtype,sample_1,sample_2,sample_3,type_origin,bg.col,type_origin_1,type_origin_2,best_k)
```


