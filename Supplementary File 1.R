
setwd('..')

library(Seurat)
library(NMF)
library(Seurat)
library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)
library(ggplot2)

scRNA=readRDS('./scRNA.rds')

# Quality Control
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
col.num <- length(levels(scRNA@active.ident))
minGene=200; maxGene=4000; pctMT=15
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))

# Standardization
# Dimensionality and Clustering
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
pc.num=1:11
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
scRNA <- FindClusters(scRNA, resolution = 1.2)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

# tSNE UMAP
scRNA = RunTSNE(scRNA, dims = pc.num)
scRNA <- RunUMAP(scRNA, dims = pc.num)
DimPlot(scRNA, reduction = "tsne",label = T) 
DimPlot(scRNA, reduction = "umap",label=T) 
scRNA=readRDS('scRNA_anno_self.rds')

# pheatmap
gene=read.table('Disulfidptosis.csv')
DoHeatmap(subset(scRNA,downsample=50,),features = gene,group.by = 'celltype',assay='RNA',slot = 'data',angle = 0,hjust = 0.6,group.bar.height = 0.04,
          group.colors =c('#313c63','#b42e20','#ebc03e','#377b4c','#7bc7cd','#5d84a4'),lines.width = 10,size = 3)+
  scale_fill_gradientn(colors=c('white','firebrick3'),na.value = 'white')
dev.off()


# Subgroup categorization, exemplified by CAFs
scRNA=readRDS('./scRNA_anno_self.rds')

#
scRNA_stromal=subset(scRNA_stromal,tissue_type %in% 'Tumor')
scRNA_fibro=subset(scRNA_stromal, stromal_cell_subtype %in% 'Fibroblasts')
gene=read.table('Disulfidptosis.csv')
scRNA_fibro$id=colnames(scRNA_fibro)
scRNA_fibro_aa=scRNA_fibro[rownames(scRNA_fibro) %in% gene,]
df <- scRNA_fibro_aa@assays$RNA@data
df=as.data.frame(df)
df <- df[rowMeans(df) !=0,  ]
df <- df[,colMeans(df) !=0 ]
scRNA_fibro=subset(scRNA_fibro,id %in% colnames(df))

# NMF
res <- nmf(df, 10, method = "snmf/r", seed = 'nndsvd')
scRNA_fibro@reductions$nmf <- scRNA_fibro@reductions$pca
scRNA_fibro@reductions$nmf@cell.embeddings <- t(coef(res))    
scRNA_fibro@reductions$nmf@feature.loadings <- basis(res)  
scRNA_fibro.nmf <- RunUMAP(scRNA_fibro, reduction = 'nmf', dims = 1:10) 
scRNA_fibro.nmf <- FindNeighbors(scRNA_fibro.nmf,reduction = 'nmf', dims = 1:10)
scRNA_fibro.nmf <- FindClusters(scRNA_fibro.nmf)
saveRDS(scRNA_fibro.nmf,file ='scRNA_fibro_NMF.RDS')


# Summary of key methodologies
# CellChat
library(CellChat)
scRNA=readRDS('./scRNA.rds')
scRNA_chat <- subset(scRNA, tissue_type == 'Tumor')
meta = scRNA_chat@meta.data
data_input <- as.matrix(scRNA_chat@assays$RNA@data)
identical(colnames(data_input),rownames(meta))
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")
CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net<- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()


## monocle
scRNA=readRDS('./scRNA.rds')
data=as.matrix(scRNA@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
plot_pseudotime_heatmap(mycds[ss,], cluster_rows = T,show_rownames = TRUE,return_heatmap = TRUE)
plot_cell_trajectory(mycds, color_by = "State")
plot_cell_trajectory(mycds, color_by = "NMF_cluster")


## SCENIC
library(SCENIC)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(cowplot)

scRNAsub=readRDS('./scRNA.rds')
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
cellInfo <- scRNAsub@meta.data
setwd("./SCENIC") 
Idents(scRNAsub)=scRNAsub$NMF_celltype
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
Idents(scRNAsub) <- "NMF_cluster"

#
org='hgnc'
dbDir="./SCENIC"
myDatasetTitle="SCENIC Challenge"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v8
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs = dbs)
scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "colVars.Rds"
scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr" = "hg19-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v9"
saveRDS(scenicOptions, file="scenicOptions.Rds") 
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.1*ncol(exprMat),minSamples=ncol(exprMat)*.1)
exprMat_filtered <- exprMat[genesKept,]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50"))

#
library(foreach)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="scenicOptions.Rds")
scenicOptions@settings$seed <- 123
runSCENIC_4_aucell_binarize(scenicOptions)
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC
aucell_regulonAUC.t <- t(aucell_regulonAUC@assays@data$AUC)
fibro.scenic <- scRNAsub
fibro.scenic@meta.data <- cbind(fibro.scenic@meta.data, aucell_regulonAUC.t[rownames(fibro.scenic@meta.data),])
DimPlot(fibro.scenic, reduction = "umap")
FeaturePlot(fibro.scenic, reduction = "umap", features = colnames(fibro.scenic@meta.data)[20:27], cols = c("yellow", "red"))

#
Idents(fibro.scenic)=fibro.scenic$NMF_celltype
cells.ord.cluster <- fibro.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
regulon.scores <- t(aucell_regulonAUC.t[names(cells.ord.cluster),])
regulon.scores.log <- log(regulon.scores +1)
regulon.scores.scaled <- scale(regulon.scores)
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(regulon.scores, 1, cal_z_score))
cluster.col <- data.frame(fibro.scenic@active.ident, row.names = names(fibro.scenic@active.ident))
colnames(cluster.col) <- "group"

# 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$NMF_celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100))


## KEGG
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)
library(R.utils)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
gene <- read.csv("deg_fibro2_kegg_0.csv")
kegg <- list()

for (i in 0:10) {
  a <- gene %>% 
    filter(cluster == i)
  genelist <- bitr(a$gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
  gene_f <- a %>% 
    right_join(genelist, by = c("gene"="SYMBOL" ))
  kegg[[i+1]] <- enrichKEGG(gene =  gene_f$ENTREZID,
                            organism     = 'hsa',
                            pvalueCutoff = 0.1,
                            qvalueCutoff = 0.1)
}

names(kegg) <- paste0('C',1:11)
celltype_pathway <- data.frame()

# 
for (i in 1:10) {
  kegg_single <- kegg[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column("name")
  kegg_single<-  kegg_single %>% 
    mutate(log.p.adjust = -log( kegg_single$p.adjust)) %>% 
    mutate(celltype_pathway = names(kegg)[i])
  celltype_pathway <- rbind(celltype_pathway, kegg_single)               
}

# 
celltype_pathway=celltype_pathway[celltype_pathway$p.adjust<0.05,]
rownames(celltype_pathway)=NULL
celltype_pathway <-celltype_pathway %>% dplyr::select(3,11,12) 
celltype_pathway$Description <- factor(celltype_pathway$Description)
ggplot(celltype_pathway,aes(x=celltype_pathway,y=Description))+
  geom_tile(aes(fill=log.p.adjust),color = "grey")+
  scale_fill_gradient(low = "#F7ED7A",high = "#E03220")+
  theme_classic()+theme(axis.text.x = element_text(angle=0))+
  labs(title = "KEGG Pathways",
       x = element_blank(),
       y=element_blank()
	   )


## scMetabolism
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(pheatmap)

# 
scRNA_mac.nmf_meta<-sc.metabolism.Seurat(obj = scRNA_mac.nmf, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
df = data.frame(t(scRNA_mac.nmf_meta@assays[["METABOLISM"]][["score"]]))
names(scRNA_mac.nmf_meta$NMF_celltype)
df = df[names(scRNA_mac.nmf_meta$NMF_celltype),]
df$NMF_celltype <- scRNA_mac.nmf_meta$NMF_celltype
avg_df =aggregate(df[,1:ncol(df)-1],list(df$NMF_celltype),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
avg_df <- avg_df[,c(1,3,4,2)]
order_index <- order(avg_df[,2])
avg_df <- avg_df[order_index, ]
avg_df_sec <- sample_n(avg_df, 30)
pheatmap(avg_df_sec, show_colnames = T,scale='row', cluster_rows = F,
         color=colorRampPalette(c('blue','white',"red"))(100),
         cluster_cols = F)
pheatmap::pheatmap(avg_df,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100))
input.pathway <- rownames(scRNA_mac.nmf_meta@assays[["METABOLISM"]][["score"]])[1:30]
DotPlot.metabolism(obj = scRNA_mac.nmf_meta, pathway = input.pathway, phenotype = "NMF_celltype", norm = "y")


## Bulk
library(genefilter)
library(Biobase)
library(edgeR)
library(ggpubr)
library(tidyverse)
library(GSVA)
library(Biobase)
library(survminer)
library(survival)
library(dplyr)
library(caret)
library(cowplot)

# Take TCGA as an example
load('./KIRC.Rdata') 

# gsva
gsva_matrix<- gsva(as.matrix(exprSet_tcga_mRNA), gene_set, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix=as.data.frame(gsva_matrix)

#
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))
rt=as.data.frame(t(gsva_matrix))
rt$group=group_list
colnames(rt)

ggboxplot(rt,x = "group",
          y = "Disulfidptosis",
          color = "black",
          fill = "group",
          xlab = "group",
          ylab = "Disulfidptosis", palette=c('#b42e20','#ebc03e')
)+stat_compare_means()

gene_set = read.table(file='elltype.txt', header = T, sep = '\t',stringsAsFactors = F)

list <- list()
for(i in 1:length(unique(gene_set$Aggre_celltype))){
  list[[i]] <- gene_set$marker[gene_set$Aggre_celltype== (unique(gene_set$Aggre_celltype)[i])]
}
names(list)<- unique(gene_set$Aggre_celltype)

gsva_matrix<- gsva(as.matrix(exprSet_tcga_mRNA), list, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix=as.data.frame(gsva_matrix)
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))
rt=as.data.frame(t(gsva_matrix))
rt$group=group_list
rt=tidyr::pivot_longer(rt,cols = -c('group'),names_to = "Aggre_celltype",values_to = 'Abundance')
ggboxplot(rt,x = "Aggre_celltype",y = "Abundance",color = "black",fill = "group", xlab = "group", ylab = "Abundance", palette=c('#b42e20','#ebc03e')) +
  stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test", hide.ns = T, size = 4.5) +
  theme(axis.text.x = element_text(angle =45, hjust = 1, vjust = 1))

# Cox
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
data=as.data.frame(t(data))
suv=read.table('./survival.tsv',row.names = 1,header = T,check.names = F)
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
pFilter=1
outTab=data.frame()
sigGenes=c("futime","fustat")

#
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
rt$futime=rt$futime/30
immune_p=c()
immune_figure=list()

dir.create('k_m')
for (i in rownames(gsva_matrix)) {
  res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(data[,i]<= cutoff, "Low", "High")
  data=rt
  data$group=Type
  data$group=factor(data$group, levels=c("Low", "High"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = data)
  length=length(levels(factor(data[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  immune_p=c(immune_p,pValue)
  if(pValue<0.05){
    fit <- survfit(Surv(futime, fustat) ~ group, data = data)
    bioCol=c('#ebc03e','#b42e20')
    bioCol=bioCol[1:length]
    p=ggsurvplot(fit, 
                 data=data,
                 conf.int=F,
                 pval=pValue,
                 pval.size=6,
                 legend.title=i,
                 legend.labs=levels(factor(data[,"group"])),
                 legend = c(0.7, 0.8),
                 font.legend=12,
                 xlab="Time(Months)",
                 palette = bioCol,
                 #surv.median.line = "hv",
                 risk.table= T,
                 cumevents=F,
                 risk.table.height=.15
                 )
    ggsave2(filename = paste0('./k_m/',i,'.pdf'),width = 4,height = 4)
  }
}

# Bubble heatmap
rt1=read.table('./uniCox_tcga.txt',header = T)
rt2=read.table('./uniCox_acrg.txt',header = T)
rt1$Dataset='TCGA'
rt2$Dataset='GEO'
rt=rbind(rt1,rt2)
rt$`-logP`=-log(rt$pvalue)
rt$logHR=log(rt$HR)
library(ggplot2)
rt$logHR=ifelse(rt$logHR>3,3,rt$logHR)
rt$logHR=ifelse(rt$logHR< -3,-3,rt$logHR)
ggplot(data=rt)+ geom_point(aes(y=Dataset,x=id,fill=logHR,size= `-logP`), color='black',shape=21,stroke=1.5)+
  scale_fill_gradientn(colours = c('#403D76','#E3B635','#C02E20'),limits=c(-3,3))+
  theme_bw()+ theme(axis.text.x=element_text(angle=90,hjust=1), panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  scale_size_area(breaks=c(1,2,3))
  