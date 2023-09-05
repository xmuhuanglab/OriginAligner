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
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/OriginAligner.PNG)
##### load R package: Required R package can be found in identify the cirrhosis-like subpopulations/R package
##### Step 1 ： calculate the hvg based on dispersion and weight（function--Screen_hvg）
##### Step 2 ： select the hvg by yourself （function--Plot_hvg） 
##### Step 3 ： calculate KNN and visualize  (function--OriginAligner), you had better choice the best result to visualize
##### Step 4 ： calculate the overlap of markers and the Jaccard similarity（function--Venn_Jaccard）

### Part2: Discovering potential Cirrhosis-like subpopulations
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/OriginAligner.PNG)



## Run the code
```
source('identify the cirrhosis-like subpopulations/OriginAligner_main')
Screen_hvg(seurat_obj,subtype,gene_number)
Plot_hvg(disp_threshold,weight_threshold)
```
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/example_hvg_screen.png)
```
OriginAligner(seurat_obj,subtype,sample_1,sample_2,sample_3,type_origin,bg.col,type_origin_1,type_origin_2,best_k)
```
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/SIM_example.png)


