# SIM-scRNA
![Image text](https://github.com/xmuhuanglab/SIM-scRNA/blob/main/images/figure_pipeline.png)
## OriginAligner: a method to phenotypic comparison of different sample sources
### Introduction:
Here we designed a modified statistical method OriginAligner based on the hypothesis that cells resemble transcriptome characteristics in different samples might 
have similar microenvironments or niches.
### Workflow:
![Image text](https://github.com/xmuhuanglab/OriginAligner/blob/main/images/OriginAligner.png)

### How to cite OrignAligner
Integrated analysis of single-cell transcriptome of liver cancer and cirrhosis reveals cell lineage similarity and prognostic-associated subpopulations. bioRxiv, 2022.2011.2003.515124. doi:10.1101/2022.11.03.515124 (Preprint)

### Step
#### Step 1 ： calculate the hvg based on dispersion and weight（function--Screen_hvg）
#### Step 2 ： select the hvg by yourself （function--Plot_hvg） 
#### Step 3 ： calculate KNN and visualize  （function--OriginAligner）， you had better choice the best result to visualize
#### Step 4 ： calculate the overlap of marker and calculate the Jaccard similarity to validate above found  （function--Venn_Jaccard）
#### load R package

### Pull the docker (the env about R version)
```
docker pull sluo112211/sim-docker
```
```
library(ggpointdensity)
library(cowplot)
library(viridis)
library(Seurat)
library(dplyr)
library(tibble)
library(reshape2)
library(patchwork)
library(stringr)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(data.table)
library(Hmisc)
library(tidyr)
library(class)    
library(RColorBrewer)
library(ggvenn)
library(rgl)
library(ggforce)
library(reshape2)
library(reshape)
library(gplots)
library(circlize)
```
```
#####  the requir of input
##  R version 4.1.0
#  seurat_obj :  total seurat obj 
#  subtype :  subset of the seurat object (the KNN best result or you want to visualize）   the colnames must be called celltype    eg: 'Hepatocyte'
#  sample_1 & sample_2 & sample_3:   sample, the colnames must be called data_type   eg： 'Cirrhosis'
#  k : the number of marker you want to overlap
#  type_originII_1 & type_originII_2: the subset you want to alignment, eg: 'HCC_PGA5+ Hepa' 
#  gene_number : you should set the number of hvg you want based on diff dataset conditions
#  disp_threshold、weight_threshold :  you should set the threshold based on diff dataset conditions, and could judge based on the results of Plot_hvg
#  row.col : the order of circle plot 
#  grid.col : each circle linecolor of each subtype

load('Example data/example_dataset.RData')
seurat_obj=sub
subtype="Hepatocyte"
sample_1='HCC'
sample_2='Cirrhosis'
sample_3='Healthy'   ####   circle plot need
k=100
type_originII_1='HCC_PGA5+ Hepa'
type_originII_2='Cirrhosis_Hepa KNG1'
type_origin='type_originII'    ####  the colnames
row.col=c("white","white","#40a9ff")
grid.col = c('HCC_COX7A1+ Hepa' = "#fa541c", 'HCC_CYP2E1+ Hepa' = "#fa541c", 'HCC_HSPA1B+ Hepa' = "#fa541c",'HCC_IGFBP3+ Hepa' = "#fa541c",'HCC_PGA5+ Hepa' = "#fa541c",
             'Cirrhosis_Hepa CRP'="#ffa940",'Cirrhosis_Hepa IGLC3'="#ffa940",'Cirrhosis_Hepa KNG1'="#ffa940",'Healthy_Hepa CRP' = "#a0d911", 'Healthy_Hepa IGLC3' = "#a0d911", 'Healthy_Hepa KNG1' = "#a0d911")
gene_number=300
disp_threshold=1.5     ## you could judge based on the results of Plot_hvg
weight_threshold=0.003  ##  you could judge based on the results of Plot_hvg
```

#### Run the code
```
source('R/R_function')
Screen_hvg(seurat_obj,subtype,gene_number)
Plot_hvg(disp_threshold,weight_threshold)
OriginAligner(seurat_obj,subtype,sample_1,sample_2,sample_3,type_origin,grid.col,row.col)
Venn_Jaccard(seurat_obj,subtype,sample_1,sample_2,type_originII,k,type_originII_1,type_originII_2)
```
