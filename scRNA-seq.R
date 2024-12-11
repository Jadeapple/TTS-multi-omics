##reanalysis of the single RNA-seq data downloaded from GSE189358

source activate seurat

library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)


samples=list.files("/mnt/sda/jwz/data/TTS/GSE189358/cellranger")
samples

names(samples) = c("Control","ISO")
samples


##批量读取数据并创建Seurat对象
scRNAlist <- list()
for(i in 1:length(samples)){
  counts <- Read10X(data.dir = samples[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3)
  #scRNAlist[[i]] <- CreateSeuratObject(counts)
  }

  
for(i in 1:length(scRNAlist)){
sc <- scRNAlist[[i]]
# 计算线粒体基因比例
sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
##计算核糖体基因比例
sc[["ribo_percent"]] <- PercentageFeatureSet(sc, pattern = "^Rp[sl]")
# 计算红细胞比例
HB_genes <- c("Hba-a1","Hba-a2","Hbegf","Hbb-y","Hbp1","Hbs1l","Hba-x")
HB_m <- match(HB_genes, rownames(sc@assays$RNA))
HB_genes <- rownames(sc@assays$RNA)[HB_m] 
HB_genes <- HB_genes[!is.na(HB_genes)] 
sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes) 
scRNAlist[[i]] <- sc
rm(sc)
}
scRNAlist[[1]]@meta.data[1:6,1:6]


pdf("./qc/qc-before.pdf")
VlnPlot(scRNAlist[[1]], features = c("nFeature_RNA", "nCount_RNA", "mt_percent","ribo_percent"), ncol = 4, pt.size = 0)
VlnPlot(scRNAlist[[2]], features = c("nFeature_RNA", "nCount_RNA", "mt_percent","ribo_percent"), ncol = 4, pt.size = 0)
dev.off()


##批量过滤细胞
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
	x <- subset(x, 
	subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & 
	mt_percent < 25
)})

pdf("./qc/qc-filter.pdf")
VlnPlot(scRNAlist[[1]], features = c("nFeature_RNA", "nCount_RNA", "mt_percent","ribo_percent"), ncol = 4, pt.size = 0)
VlnPlot(scRNAlist[[2]], features = c("nFeature_RNA", "nCount_RNA", "mt_percent","ribo_percent"), ncol = 4, pt.size = 0)
dev.off()


##利用merge()合并Seurat对象
scRNAlist <- merge(scRNAlist[[1]],y=c(scRNAlist[[2]]))
head(scRNAlist@meta.data,5)


##增加metadata分组信息
library("stringr")
phe=str_split(rownames(scRNAlist@meta.data),'_',simplify = T)
head(phe)
scRNAlist@meta.data$group=phe[,1]
head(scRNAlist@meta.data)


## 统计细胞数
table(scRNAlist[[]]$orig.ident)
dim(scRNAlist)   #查看基因数和细胞总数
table(scRNAlist@meta.data$orig.ident)  #查看每个样本的细胞数



##1.数据标准化，LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
scRNAlist <- NormalizeData(scRNAlist)
##2.鉴定高变基因
scRNAlist <- FindVariableFeatures(scRNAlist,selection.method = "vst",nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(scRNAlist), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(scRNAlist)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("VariableFeaturePlot.pdf")
plot1 + plot2
dev.off()
##3.ScaleData()
scRNAlist <- ScaleData(scRNAlist)

##降维可视化
##4. PCA降维
scRNAlist <- RunPCA(scRNAlist)        

##查看PCA结果
print(scRNAlist[["pca"]], dims = 1:5, nfeatures = 5)
##可视化PCA结果
##展示主成分基因
pdf("PCA-gene.pdf")
VizDimLoadings(scRNAlist, dims = 1:2, reduction = "pca")
dev.off()
##热图，第一主成分
pdf("PCA-heatmap.pdf")
DimHeatmap(scRNAlist, dims = 1, cells = 500, balanced = TRUE) 
dev.off()


##确定数据的维度
pdf("ElbowPlot-50.pdf")
ElbowPlot(scRNAlist, ndims = 50)
dev.off()

##保留20个主成分，累加贡献度达到90%
scRNAlist <- RunUMAP(scRNAlist,reduction = "pca", dims = 1:20, verbose = F)
scRNAlist <- RunTSNE(scRNAlist,reduction = "pca", dims = 1:20, verbose = F)


pdf("umap-before-integration-group.pdf")
DimPlot(scRNAlist,reduction = "umap", group.by = "group") + plot_annotation(title = "before integration")
dev.off()
pdf("tsne-before-integration-1.pdf")
DimPlot(scRNAlist,reduction = "tsne", group.by = "group") + plot_annotation(title = "before integration")
dev.off()



##harmony整合
scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]

# 聚类
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- scRNA_harmony %>% FindClusters(resolution = 0.1)
# umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)



pdf("umap-after-integration-group.pdf")
DimPlot(scRNA_harmony,reduction = "umap", group.by = "group") + plot_annotation(title = "after integration")
dev.off()
pdf("tsne-after-integration-group.pdf")
DimPlot(scRNA_harmony,reduction = "tsne", group.by = "group") + plot_annotation(title = "after integration")
dev.off()
pdf("harmony_umap-cluster.pdf")
DimPlot(scRNA_harmony,reduction = "umap",label = TRUE) 
dev.off()
pdf("harmony_tsne-cluster.pdf")
DimPlot(scRNA_harmony,reduction = "tsne",label = TRUE) 
dev.off()

# split.by 展示分组聚类
pdf("harmony_umap-split.pdf")
DimPlot(scRNA_harmony, reduction = "umap", split.by = "group",ncol=2,label = T)
dev.off()
pdf("harmony_tsne-split.pdf")
DimPlot(scRNA_harmony, reduction = "tsne", split.by = "group",ncol=2,label = T)
dev.off()


##寻找marker基因，注释细胞类型
markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_df = markers %>% group_by(cluster)
write.table(markers_df,file="all-markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
head(markers_df)
write.table(markers_df,file="markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)



#每个聚类前10个差异基因表达热图
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("marker-heatmap.pdf")
DoHeatmap(scRNA_harmony, features = top10$gene) + NoLegend()
dev.off()


##人工注释
#加上注释
new.cluster.ids <- c("Monocytes", "B cells", "Macrophages", "T cells", "Neutrophils", "NK cells", "T helper cells", "Dendritic Cells", "Monocytes", "Endothelial cells")

names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
#在metadata中，添加Celltype信息
scRNA_harmony$celltype <- Idents(scRNA_harmony)
table(scRNA_harmony$celltype,scRNA_harmony@meta.data$seurat_clusters)
##可视化
pdf("harmony-cluster-name-umap.pdf")
DimPlot(scRNA_harmony, reduction = "umap", group.by = "celltype", label = FALSE, pt.size = 0.5)
dev.off()

pdf("harmony-cluster-split.pdf")
DimPlot(scRNA_harmony, reduction = "umap", group.by = "celltype", split.by="group", ncol=1, label = FALSE, pt.size = 0.5)
dev.off()


##细胞比例计算及可视化
table(scRNA_harmony$orig.ident)  #查看各样本细胞数
table(Idents(scRNA_harmony), scRNA_harmony$orig.ident)  #各样本不同细胞群细胞数

table(scRNA_harmony$group)  #查看各组别细胞数
table(Idents(scRNA_harmony), scRNA_harmony$group)  #各组别不同细胞群细胞数

Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$group), margin = 2)#计算各组别不同细胞群比例

Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2=factor(Cellratio$Var2,levels=c("WT","BLM","AGGF1"))
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
pdf("cell-ratio.pdf")
ggplot(Cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Group',y = 'Cell Type Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()


