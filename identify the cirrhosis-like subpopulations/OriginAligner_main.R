# Introduce the step
##Step 1 ： calculate the hvg based on dispersion and weight（function--Screen_hvg）
##Step 2 ： select the hvg by yourself （function--Plot_hvg） 
##Step 3 ： calculate KNN and visualize  （function--OriginAligner）， you had better choice the best result to visualize
##Step 4 ： calculate the overlap of marker and calculate the Jaccard similarity to validate above found  （function--Venn_Jaccard）


##      Step 1 and 2: Disp and weight calculation

Screen_hvg = function(seurat_obj,subtype,gene_number){
  cells=rownames(seurat_obj@meta.data)[which(seurat_obj$celltype==subtype)]
  sub=subset(seurat_obj,cells=cells)
  exp=as.matrix(sub@assays$RNA@counts)
  
  #normalize
  n=gene_number
  matrix=exp
  head(matrix,50)
  matrix= t(t(matrix)/apply(matrix, 2, sum))*1000000
  head(matrix,50)
  
  # filter low quality cell
  pqcells = is.na(apply(matrix>0, 2, sum)) | apply(matrix>0, 2, sum) <= 10  
  num_pqcells = length(which(pqcells == TRUE))
  matrix = matrix[,!pqcells]
  
  # log2
  matrix = log(matrix+1,2)
  matrix = data.matrix(matrix)
  
  # calculate pre-batch corrected gene counts
  counts = apply(matrix>0, 2, sum)
  matrix[which(matrix<0)] = 0
  census_normalize = function(matrix, counts) {
    xnl = 2^data.matrix(matrix) - 1
    rs = apply(xnl, 2, sum)
    rnorm = t(t(xnl) * counts/rs)
    A = log(rnorm+1,2)
    return(A)
  }
  matrix2 = census_normalize(matrix, counts)
  A=matrix2 
  n_expr = rowSums(A > 0);
  A_filt = A[n_expr >= 0.05 * ncol(A),];
  vars = apply(A_filt, 1, var);
  means = apply(A_filt, 1, mean);
  disp = vars / means    ###  vars devide means
  last_disp = tail(sort(disp), n)[1]
  A_filt = A_filt[disp >= last_disp,];
  matrix2.mvg = A_filt
  dim(matrix2.mvg)
  
  # weight calculation
  metadata=sub@meta.data
  anno=metadata[,c('barcode','type_originII')]
  aa=t(exp)
  combine=cbind(anno,aa)  
  mat=matrix(data=0,nrow=length(table(combine$type_originII)),ncol=nrow(exp))
  for(i in 1:length(table(combine$type_originII))){
    mat[i,]=colMeans(combine[which(combine$type_originII==names(table(combine$type_originII))[i]),3:(nrow(exp)+2)])
  }
  colnames(mat)=names(colMeans(combine[which(combine$type_originII==names(table(combine$type_originII))[2]),3:(nrow(exp)+2)]))
  rownames(mat)=names(table(combine$type_originII))
  mat=as.matrix(t(mat))
  head(mat,50)
  min_max_norm = function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  p_value = function(x){
    x / sum(x)
  }
  entropy = function(x){
    n = length(x)
    (-1 / log2(n)) * (sum( x * ifelse(log2(x)==-Inf, 0, log2(x)) ))
  }
  weight = function(x){
    (1-x) / (length(x)-sum(x))
  }
  head(mat)
  tb.dt = as_tibble(t(mat[intersect(rownames(mat),rownames(matrix2.mvg)),]))
  tb.dt = tb.dt %>% mutate(across(c(1:length(intersect(rownames(mat),rownames(matrix2.mvg)))), min_max_norm)) %>% 
    mutate(across(c(1:ncol(tb.dt)), p_value)) %>%
    summarise(across(c(1:ncol(tb.dt)), entropy))
  w_dat = tb.dt %>% weight
  gene_screen=cbind(t(w_dat),as.data.frame(disp)[rownames(t(w_dat)),])
  colnames(gene_screen)=c("weight","disp")
  gene_screen=as.data.frame(gene_screen)
  # return(gene_screen)
  write.csv(gene_screen,file=paste('gene_screen',subtype,'.csv',sep = '_'))
}
Plot_hvg = function(disp_threshold,weight_threshold){
  gene_screen=read.csv(paste('gene_screen',subtype,'.csv',sep = '_'))
  gene_screen$color='1'
  gene_screen$color[which(gene_screen$disp>disp_threshold & gene_screen$weight>weight_threshold)]=gsub('1','hvg',gene_screen$color[which(gene_screen$disp>disp_threshold & gene_screen$weight>weight_threshold)])
  gene_screen$color[which(gene_screen$color=='1')]="Lowly weight"
  fig=ggplot(gene_screen, aes(x=weight, y=disp,color=color)) + 
    geom_point(size=2,shape = 2)+  scale_color_manual(
      labels = c(
        paste("Hvg (",as.character(table(gene_screen$color)["hvg"]),")",sep =''), 
        paste("Lowly weight (",as.character(table(gene_screen$color)["Lowly weight"]),")",sep ='')),
      values = c("#f5222d","#8c8c8c"))+
    theme_bw()+ theme(panel.grid=element_blank())+  labs(x = "Each Gene Weight", y = "Each Gene disp",title = "Screen the gene set")+
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(),
      axis.text.x = element_text(hjust = 0.5, color = "black",size = 12),
      axis.text.y = element_text(angle = 90,hjust=0.5, color = "black",size = 12),
      axis.title.x = element_text(vjust = 0, size = 16),
      axis.title.y = element_text(vjust = 2, size = 16),
      plot.title = element_text(hjust = 0.4,size = 16),
      legend.key = element_blank(),
      legend.title = element_text(size = 14, hjust = 0.2),
      legend.text = element_text(size = 13),
      legend.position = c(0.2,0.88)
    )
  return(fig)
}

##      Step 3: OriginAligner 

OriginAligner = function(seurat_obj,subtype,sample_1,sample_2,sample_3,type_origin,bg.col,type_origin_1,type_origin_2,best_k){
  ####  gene
  gene_screen=read.csv(paste('gene_screen',subtype,'.csv',sep = '_'))
  gene=gene_screen$X
  ####  subet
  cells=rownames(seurat_obj@meta.data)[which(seurat_obj$celltype==subtype)]
  sub=subset(seurat_obj,cells=cells)
  kk=sub@meta.data
  matrix=as.matrix(sub@assays$RNA@data[gene,])  #  ncol(mat)
  metadata=sub@meta.data
  #### test
  tt=rownames(kk)[which(kk$data_type==sample_1)]
  test=t(matrix[gene,intersect(tt,colnames(matrix))])
  test_lab=as.factor(metadata[intersect(tt,colnames(matrix)),type_origin])
  #### train
  cc=rownames(kk)[which(kk$data_type==sample_2|kk$data_type==sample_3)]
  train=t(matrix[gene,intersect(cc,colnames(matrix))])
  train_lab =as.factor(metadata[intersect(cc,colnames(matrix)),type_origin])
  
  # KNN calculation
  pre=knn(train, test,
          cl=train_lab, k = best_k,
          l = 0, prob =FALSE, use.all = TRUE)
  
  data=table(pre,test_lab)
  input=as.matrix(data)
  
  df = data.frame(factors = rep(c(rownames(input),type_origin_1),100001),x =rep(seq(0,1,0.00001),nrow(input)+1), y = runif(100001*(nrow(input)+1)))

  
  par(mar = c(1, 1, 1, 1) )
  circos.initialize(factors = df$factors, x = df$x)
  circos.track(factors = df$factors,x = df$x, y = df$y,bg.col = bg.col,panel.fun = function(x,y){
  circos.text(x = get.cell.meta.data("xcenter"),y = get.cell.meta.data("cell.ylim")[2] + uy(3,"mm"),labels = get.cell.meta.data("sector.index"))})
  a=0
  for (i in 1:nrow(input)){
      ref_name=rownames(input)[i]
      if (ref_name==type_origin_2){
          head_num=a
          end_num=round((input[i,type_origin_1]/sum(input[,type_origin_1])),6)+head_num
          ref_num=round((input[i,type_origin_1]/sum(input[i,])),6)
          circos.link(type_origin_1, c(head_num,end_num), ref_name,c(0,ref_num),col='#40a9ff',h = 0.8,border="black",lty = 1,lwd = 4)
          a=end_num
       }   
      else {
     head_num=a
     end_num=round((input[i,type_origin_1]/sum(input[,type_origin_1])),6)+head_num
     ref_num=round((input[i,type_origin_1]/sum(input[i,])),6)
     # circos.link(type_origin_1, c(head_num,end_num), ref_name,c(0,ref_num),col='#FFF0F5',h = 0.8)
     a=end_num
     }
  }
}


##      Step 4:  Marker overlap (Venn plot) and Jaccard similarity
Venn_Jaccard=function(seurat_obj,subtype,sample_1,sample_2,type_originII,k,type_origin_1,type_origin_2){
  kk=seurat_obj@meta.data[which(seurat_obj$celltype==subtype),]
  cells=rownames(kk)[which(kk$data_type==sample_1)]
  sub=subset(seurat_obj,cells=cells)
  Idents(sub)=sub$type_originII
  marker_1 = FindAllMarkers(sub, logfc.threshold = 0.25, min.pct = 0, 
                            only.pos = TRUE, test.use = "wilcox")
  table(marker_1$cluster)
  marker_1 = marker_1 %>% group_by(cluster) %>% top_n(n = k, wt = avg_log2FC)
  kk=seurat_obj@meta.data[which(seurat_obj$celltype==subtype),]
  cells=rownames(kk)[which(kk$data_type==sample_2)]
  sub=subset(seurat_obj,cells=cells)
  Idents(sub)=sub$type_originII
  marker_2 = FindAllMarkers(sub, logfc.threshold = 0.25, min.pct = 0, 
                            only.pos = TRUE, test.use = "wilcox")
  table(marker_1$cluster)
  marker_2 = marker_2 %>% group_by(cluster) %>% top_n(n = k, wt = avg_log2FC)
  x = list(
    subset_1 = marker_1[which(marker_1$cluster==type_origin_1),]$gene,
    subset_2 = marker_2[which(marker_2$cluster==type_origin_2),]$gene
  )
  
  opar = par(family = "Roboto Condensed")
  mypal=c("#ff7875","#ffc069")
  pig=ggvenn(x,fill_color=mypal,fill_alpha = .6,stroke_linetype = "longdash",set_name_size = 7,text_size=9)
  pdf(paste("Venn",type_origin_1,type_origin_2,'.pdf'))
  pig
  dev.off()
  
  mat=matrix(data=0,nrow = length(names(table(marker_1$cluster))),ncol=length(names(table(marker_2$cluster))))
  for (i in 1:length(names(table(marker_1$cluster)))){
    for (j in 1:length(names(table(marker_2$cluster)))) {
      signature_1=marker_1$gene[which(marker_1$cluster==names(table(marker_1$cluster)[i]))] 
      signature_2=marker_2$gene[which(marker_2$cluster==names(table(marker_2$cluster)[j]))]
      hhh=function(a, b) {
        intersection = length(intersect(a, b))
        union = length(a) + length(b) - intersection
        return (intersection/union)
      }
      mat[i,j]=hhh(signature_1,signature_2)               
    }
  }
  rownames(mat)=names(table(marker_1$cluster)[1:length(names(table(marker_1$cluster)))])
  colnames(mat)=names(table(marker_2$cluster)[1:length(names(table(marker_2$cluster)))])
  class(mat)
  mat=as.matrix(mat)
  write.csv(mat,file=paste('Jaccard',type_originII_1,type_originII_2,".csv"))
  library(RColorBrewer)
  my_group = as.numeric(as.factor(substr(rownames(mat), 1 , 1)))
  colSide = brewer.pal(9, "Set1")[my_group]
  colMain = colorRampPalette(brewer.pal(8, "Blues"))(25)
  nn=melt(mat)
  colnames(nn)=c('type_origin_1','type_origin_2',"Jaccard_similarity")
  pic_name=ggplot(nn, aes(type_origin_2, type_origin_1, fill= Jaccard_similarity)) + 
    geom_tile()+
    theme(axis.text.x = element_text(colour="grey20",size=20,angle=45,hjust=0.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=22,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=22,angle=90,hjust=.5,vjust=.5,face="plain"))
  ggsave(pic_name, file=paste('jaccard_smilarity_',type_originII_1,type_originII_2,'.pdf'), width=20, height=10) 
  return(pig)
}
