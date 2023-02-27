# SIM-scRNA
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/figure1.png)
## OriginAligner: a method to phenotypic comparison of different sample sources
### Introduction:
Here we designed a modified statistical method OriginAligner based on the hypothesis that cells resemble transcriptome characteristics in different samples might 
have similar microenvironments or niches.
### Workflow:
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/the%20pipeline%20of%20OriginAligner.png)

### How to cite OrignAligner
Integrated analysis of single-cell transcriptome of liver cancer and cirrhosis reveals cell lineage similarity and prognostic-associated subpopulations. bioRxiv, 2022.2011.2003.515124. doi:10.1101/2022.11.03.515124 (Preprint)

### Step
#### Step 1 ： calculate the hvg based on dispersion and weight（function--Screen_hvg）
#### Step 2 ： select the hvg by yourself （function--Plot_hvg） 
#### Step 3 ： calculate KNN and visualize  (function--OriginAligner), you had better choice the best result to visualize
#### Step 4 ： calculate the overlap of markers and the Jaccard similarity（function--Venn_Jaccard）
#### load R package:Required R package can be found in identify the cirrhosis-like subpopulations/R_package_library

### Pull the docker (the env about R version)
```
docker pull sluo112211/sim-docker
```

```
#####  the requir of input
##  R version 4.1.0
#  seurat_obj :  total seurat obj 
#  subtype :  subset of the seurat object (the KNN best result or you want to visualize）   the colnames must be called celltype    eg: 'Hepatocyte'
#  sample_1 & sample_2 & sample_3:   sample, the colnames must be called data_type   eg： 'Cirrhosis'
#  k : the number of marker you want to overlap
#  type_origin_1 & type_origin_2: the subset you want to alignment, eg: 'HCC_PGA5+ Hepa' 
#  gene_number : you should set the number of hvg you want based on diff dataset conditions
#  disp_threshold、weight_threshold :  you should set the threshold based on diff dataset conditions, and could judge based on the results of Plot_hvg
#  row.col : the order of circle plot 
#  grid.col : each circle linecolor of each subtype

load('Example data/GSE149614_Hepa_cirr_sample.RData')
seurat_obj=GSE149614_Hepa_cirr_sample
subtype="Hepatocyte"
sample_1='HCC'
sample_2='Cirrhosis'
sample_3='Healthy'   ####   circle plot need
k=100
type_origin_1='HCC_PGA5+ Hepa'
type_origin_2='Cirrhosis_Hepa KNG1'
type_origin='type_originII'    ####  the colnames
bg.col = c("#ffa940","#ffa940","#ffa940","#fa541c","#a0d911","#a0d911","#a0d911")
gene_number=500
best_k=5 ## by knn cross validation
disp_threshold=1.3     ## you could judge based on the results of Plot_hvg
weight_threshold=0.002  ##  you could judge based on the results of Plot_hvg
```

#### Run the code
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


