source activate signac

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(biovizBase)
library(GenomicRanges)
library(future)
library(dplyr)


setwd("/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/signac")

plan("multicore", workers = 4)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM


#Creating a common peak set
# read in peak sets
peaks.Control <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339A_B467_A-out/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.ISO1d <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339A_B486_A-out/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.ISO3d <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339B_B490_A-out/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.ISO7d <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339C_B484_A-out/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.Control <- makeGRangesFromDataFrame(peaks.Control)
gr.ISO1d <- makeGRangesFromDataFrame(peaks.ISO1d)
gr.ISO3d <- makeGRangesFromDataFrame(peaks.ISO3d)
gr.ISO7d <- makeGRangesFromDataFrame(peaks.ISO7d)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.Control, gr.ISO1d, gr.ISO3d, gr.ISO7d))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#Create Fragment objects
# load metadata
md.Control <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339A_B467_A-out/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.ISO1d <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339A_B486_A-out/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.ISO3d<- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339B_B490_A-out/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.ISO7d <- read.table(
  file = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339C_B484_A-out/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]



# create fragment objects
frags.Control <- CreateFragmentObject(
  path = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339A_B467_A-out/outs/fragments.tsv.gz",
  cells = rownames(md.Control)
)
frags.ISO1d <- CreateFragmentObject(
  path = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339A_B486_A-out/outs/fragments.tsv.gz",
  cells = rownames(md.ISO1d)
)
frags.ISO3d<- CreateFragmentObject(
  path = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339B_B490_A-out/outs/fragments.tsv.gz",
  cells = rownames(md.ISO3d)
)
frags.ISO7d <- CreateFragmentObject(
  path = "/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/231339C_B484_A-out/outs/fragments.tsv.gz",
  cells = rownames(md.ISO7d)
)

#Quantify peaks in each dataset
Control.counts <- FeatureMatrix(
  fragments = frags.Control,
  features = combined.peaks,
  cells = rownames(md.Control)
)

ISO1d.counts <- FeatureMatrix(
  fragments = frags.ISO1d,
  features = combined.peaks,
  cells = rownames(md.ISO1d)
)

ISO3d.counts <- FeatureMatrix(
  fragments = frags.ISO3d,
  features = combined.peaks,
  cells = rownames(md.ISO3d)
)

ISO7d.counts <- FeatureMatrix(
  fragments = frags.ISO7d,
  features = combined.peaks,
  cells = rownames(md.ISO7d)
)

#Create the objects
Control_assay <- CreateChromatinAssay(Control.counts,fragments = frags.Control,min.cells = 10, min.features = 200)
Control <- CreateSeuratObject(Control_assay, assay = "ATAC", meta.data=md.Control)

ISO1d_assay <- CreateChromatinAssay(ISO1d.counts,fragments = frags.ISO1d,min.cells = 10, min.features = 200)
ISO1d <- CreateSeuratObject(ISO1d_assay, assay = "ATAC", meta.data=md.ISO1d)

ISO3d_assay <- CreateChromatinAssay(ISO3d.counts,fragments = frags.ISO3d,min.cells = 10, min.features = 200)
ISO3d <- CreateSeuratObject(ISO3d_assay, assay = "ATAC", meta.data=md.ISO3d)

ISO7d_assay <- CreateChromatinAssay(ISO7d.counts,fragments = frags.ISO7d,min.cells = 10, min.features = 200)
ISO7d <- CreateSeuratObject(ISO7d_assay, assay = "ATAC", meta.data=md.ISO7d)

##每个样本分别质控
#添加基因注释
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- "UCSC"
#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(Control) <- annotations

#Computing QC Metrics
#核小体条带分布模式
Control <- NucleosomeSignal(object = Control)
Control$nucleosome_group <- ifelse(Control$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("Control-FragmentHistogram.pdf")
FragmentHistogram(object = Control, group.by = 'nucleosome_group',region = 'chr1-1-10000000')
dev.off()
table(Control@meta.data$nucleosome_group)

#转录起始位点富集分数
Control <- TSSEnrichment(Control, fast = FALSE)
Control$high.tss <- ifelse(Control$TSS.enrichment > 2, 'High', 'Low')
pdf("Control-TSSPlot.pdf")
TSSPlot(Control, group.by = 'high.tss') + NoLegend()
dev.off()
table(Control@meta.data$high.tss)

Control$pct_reads_in_peaks <- Control$peak_region_fragments / Control$passed_filters * 100
Control$blacklist_ratio <- Control$blacklist_region_fragments / Control$peak_region_fragments

#散点图
pdf("Control-DensityScatter.pdf")
DensityScatter(Control, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

###画qc指标的小提琴图
pdf("Control-qc-violin.pdf")
VlnPlot(
  object = Control,
  features = c('nCount_ATAC','pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal','blacklist_ratio'),
  pt.size = 0,
  ncol = 3
)
dev.off()
##过滤，根据质控指标过滤
Control <- subset(
  x = Control,
  subset = TSS.enrichment > 2 &
	nucleosome_signal < 4 &
	peak_region_fragments > 1000 &
    peak_region_fragments < 30000 &
	nCount_ATAC > 1000 &
    nCount_ATAC < 30000 &
    pct_reads_in_peaks > 5 &
    blacklist_ratio < 0.025
)
Control

##过滤后绘图
pdf("Control-FragmentHistogram-filter.pdf")
FragmentHistogram(object = test,region = 'chr1-1-10000000')
dev.off()
pdf("Control-TSSPlot-filter.pdf")
TSSPlot(test) + NoLegend()
dev.off()
pdf("Control-DensityScatter-filter.pdf")
DensityScatter(test, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
pdf("Control-qc-violin-filter.pdf")
VlnPlot(
  object = test,
  features = c('nCount_ATAC','pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal','blacklist_ratio'),
  pt.size = 0,
  ncol = 3
)
dev.off()


##合并整合后的质控小提琴图
pdf("integrated-qc-violin-1.pdf")
VlnPlot(
  object = integrated,
  group.by="dataset",
  cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#9467BDFF"),
  features = c('nCount_ATAC','pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal','blacklist_ratio'),
  pt.size = 0,
  ncol = 3
)
dev.off()



##简单合并
#Merge objects
# add information to identify dataset of origin
Control$dataset <- 'Control'
ISO1d$dataset <- 'ISO1d'
ISO3d$dataset <- 'ISO3d'
ISO7d$dataset <- 'ISO7d'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = Control,
  y = list(ISO1d, ISO3d, ISO7d),
  add.cell.ids = c("Control", "ISO1d", "ISO3d", "ISO7d")
)
combined[["ATAC"]]
head(combined@meta.data,5)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q0")
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')

pdf("combine-umap.pdf")
DimPlot(combined, group.by = 'dataset', pt.size = 0.1,cols=c("#1F77B4FF","#D62728FF","#FF7F0EFF","#2CA02CFF"))
dev.off()


CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)


# Clusters
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3,resolution = 0.1)
pdf("combined-umap-cluster.pdf")
DimPlot(object = combined, label = TRUE,cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF")) + NoLegend()
dev.off()
pdf("combined-umap-split.pdf")
DimPlot(combined, reduction = "umap", split.by = "dataset",ncol=2,,cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
dev.off()

save(Control,ISO1d,ISO3d,ISO7d,combined,file = "combined.Rdata")


##数据整合
# Rename
Control <- RenameCells(object = Control, add.cell.id = "Control")
ISO1d<- RenameCells(object = ISO1d, add.cell.id = "ISO1d")
ISO3d <- RenameCells(object = ISO3d, add.cell.id = "ISO3d")
ISO7d <- RenameCells(object = ISO7d, add.cell.id = "ISO7d")

# compute LSI,每个样本都需计算
Control <- FindTopFeatures(Control, min.cutoff = "q0")
Control <- RunTFIDF(Control)
Control <- RunSVD(Control)

ISO1d <- FindTopFeatures(ISO1d, min.cutoff = "q0")
ISO1d <- RunTFIDF(ISO1d)
ISO1d <- RunSVD(ISO1d)

ISO3d <- FindTopFeatures(ISO3d, min.cutoff = "q0")
ISO3d <- RunTFIDF(ISO3d)
ISO3d <- RunSVD(ISO3d)

ISO7d <- FindTopFeatures(ISO7d, min.cutoff = "q0")
ISO7d <- RunTFIDF(ISO7d)
ISO7d <- RunSVD(ISO7d)


##Integration
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(Control,ISO1d,ISO3d,ISO7d),
  anchor.features = rownames(Control),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

pdf("intergration-umap.pdf")
DimPlot(integrated, group.by = "dataset",cols=c("#1F77B4FF","#D62728FF","#FF7F0EFF","#2CA02CFF"))
dev.off()



# Clusters
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3,resolution = 0.1)
pdf("integration-umap-cluster.pdf")
DimPlot(object = integrated,cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
dev.off()
pdf("integration-umap-split.pdf")
DimPlot(integrated, reduction = "umap", split.by = "dataset",ncol=2,cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
dev.off()


#基因活性表达分析
#添加基因注释
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- "UCSC"
#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(integrated) <- annotations

# Create a gene activity matrix
gene.activities <- GeneActivity(integrated)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

# 查看基因活性矩阵
integrated[['RNA']]



#  visualize the activities of canonical marker genes
DefaultAssay(integrated) <- 'RNA'

##寻找marker
markers <- FindAllMarkers(integrated, only.pos = TRUE, logfc.threshold = 0.25,min.cells.feature=0)
markers_df = markers %>% group_by(cluster)
head(markers_df)
write.table(markers_df,file="integrated-all-markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.table(markers_df,file="integrated-markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)

save(integrated,file = "integrated.Rdata")
save(markers,file = "integrated-markers.Rdata")
##
#scRNAtoolVis可视化top3 marker
markers_top3 = markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
pdf("top3-marker.pdf")
jjDotPlot(object = integrated,
          gene = markers_top3$gene,
		  plot.margin=c(0.1,0.1,0.1,0.1))
dev.off()

pdf("top3-marker-1.pdf")
jjDotPlot(object = integrated,
          gene = markers_top3$gene,
		  plot.margin=c(0.1,0.1,0.1,0.1),
		  anno = T)
dev.off()

pdf("markerVolcano.pdf")
markerVolcano(markers = markers,
              topn = 2,
              labelCol = ggsci::pal_d3("category20")(13))
dev.off()

# no facet group
pdf("marker-1.pdf")
featureCornerAxes(object = integrated,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("C1qb","Cd163", "Myh6","Actc1"),
                  aspect.ratio = 1,
				  nLayout=2)
dev.off()
pdf("marker-2.pdf")
featureCornerAxes(object = integrated,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("Fibin","Cytl1", "Sox17","Krt19"),
                  aspect.ratio = 1,
				  nLayout=2)
dev.off()
pdf("marker-3.pdf")
featureCornerAxes(object = integrated,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("Rgs4","Gzmb", "S100a9","Iglc2"),
                  aspect.ratio = 1,
				  nLayout=2)
dev.off()




pdf("marker-FeaturePlot-1.pdf")
FeaturePlot(
  object = integrated,
  features = c('Fibin', 'Camp', 'C1qb', 'Cd163', 'Ephb4', 'Nos3','Gpihbp1','Actc1','Myh6','Krt19','Krt14','Rgs5','Trpc3','Cxcr6','Cd3d','S100a9','Retnlg','Cd19','Cd79a'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()
VlnPlot(integrated,
        features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
        flip = T,stack = T )+scale_fill_d3()


##EOF
##利用scRNA-seq数据注释scATAC-seq
load("/mnt/sda/jwz/data/TTS/GSE189358/seurat/scdata-harmony-celltype.Rdata")
DefaultAssay(scRNA_harmony)
DefaultAssay(integrated)

# Identify anchors
scRNA_harmony <- FindVariableFeatures(scRNA_harmony,selection.method = "vst",nfeatures = 2000)
gene.activities <- GeneActivity(integrated,features = VariableFeatures(scRNA_harmony))
integrated[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(
  object = integrated,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)
integrated <- ScaleData(integrated, features = rownames(integrated))

transfer.anchors <- FindTransferAnchors(
  reference = scRNA_harmony, # 参考集
  query = integrated, # 目标集
  reduction = 'cca', # 算法
  reference.assay = 'RNA',
  query.assay = 'ACTIVITY',
  features = VariableFeatures(scRNA_harmony)
)
head(transfer.anchors@anchors)

#通过标签转移注释 scATAC-seq 细胞
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = scRNA_harmony$celltype,
    weight.reduction = integrated[["integrated_lsi"]], dims = 2:30)
predicted.labels[1:3,]
 
#加入metadata
integrated <- AddMetaData(integrated, metadata = predicted.labels)
#完成注释转移后，ATAC-seq的细胞获得了预测的细胞类型标签，这些标签源自scRNA-seq数据集，并存储在predicted.id字段中。
##可视化
pdf("umap-cluster-name-1.pdf")
DimPlot(object = integrated, label = TRUE,reduction = "umap",group.by = 'predicted.id') + NoLegend()
dev.off()
pdf("umap-split-name-1.pdf")
DimPlot(integrated, reduction = "umap", split.by = "dataset",ncol=2)
dev.off()
##EOF

##人工注释
##人工注释
#加上注释

new.cluster.ids <- c("Monocytes", "Fibroblasts", "Macrophages", "Endocardial cells", "Endothelial cells", "Cardiomyocytes", "Epicardial cells", "Pericytes", "T cells", "Neutrophils","B cells","Astrocytes","Adipocytes")

names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)
#在metadata中，添加Celltype信息
integrated$celltype <- Idents(integrated)
table(integrated$celltype,integrated@meta.data$seurat_clusters)
##可视化
pdf("integrated-cluster-celltype-umap.pdf")
DimPlot(integrated, reduction = "umap", group.by = "celltype", label = FALSE,cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
#DimPlot(integrated, reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
dev.off()

pdf("integrated-celltype-umap-split-1.pdf")
DimPlot(integrated, reduction = "umap", group.by = "celltype", split.by="dataset", ncol=2, label = FALSE,cols=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))+ NoLegend()
dev.off()


#保存
save(integrated,file = "integrated-celltype.Rdata")

##细胞比例计算及可视化
table(integrated$orig.ident)  #查看各样本细胞数
table(Idents(integrated), integrated$orig.ident)  #各样本不同细胞群细胞数

table(integrated$group)  #查看各组别细胞数
table(Idents(integrated), integrated$group)  #各组别不同细胞群细胞数

Cellratio <- prop.table(table(Idents(integrated), integrated$dataset), margin = 2)#计算各组别不同细胞群比例

Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
colnames(Cellratio)=c("Celltype","Group","Freq")
library(ggplot2)
pdf("cell-ratio.pdf")
ggplot(Cellratio) + geom_bar(aes(x =Group, y= Freq, fill = Celltype),stat = "identity",width = 0.5,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Group',y = 'Cell Type Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

Cellratio <- prop.table(table(Idents(integrated), integrated$orig.ident), margin = 2)#计算各样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
pdf("cell-ratio-1.pdf")
ggplot(Cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5)+ 
  theme_classic() +
  labs(x='Sample',y = 'Cell Type Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

##横向柱状图
pdf("cell-ratio-1.pdf")
ggplot(Cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Group',y = 'Cell Type Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

save(integrated,file = "integrated-celltype.Rdata")

##marker可视化
gene_marker <- c("Serpine1","Fibin","C1qb","Eltd1","Flt4","Myh6","Krt19","Rgs5","Cd3d","S100a9","Iglc2","Bpifa1","Pck1","Usp9y")
gene_marker <- unique(gene_marker)
pdf("marker.pdf")
DotPlot(integrated, features = gene_marker,assay='RNA') + coord_flip()
dev.off()
pdf("marker-VlnPlot.pdf")
VlnPlot(integrated, features = gene_marker,stack = T,pt.size = 0,flip=TRUE)+
  NoLegend()+
  theme(axis.title = element_blank())
dev.off()

##可视化maker或感兴趣基因或某个区域的ATAC-seq信号
DefaultAssay(integrated) <- 'ATAC'
pdf("Fosl1-coverage.pdf")
CoveragePlot(
  object = integrated,
  region = "Fosl1",
  extend.upstream = 5000,
  extend.downstream = 5000
)
dev.off()



#EOF
##scATAC-seq提取Monocytes进行亚群分析
Monocytes = Monocytes[,Monocytes@meta.data$seurat_clusters %in% c(0)]
#标准化及线性降维
Monocytes <- RunTFIDF(Monocytes)
Monocytes <- FindTopFeatures(Monocytes, min.cutoff = 'q0')
Monocytes <- RunSVD(object = Monocytes)
# create a new UMAP using the Monocytes embeddings
Monocytes <- RunUMAP(Monocytes, reduction = "lsi", dims = 2:30)

# Clusters
Monocytes <- FindNeighbors(object = Monocytes, reduction = 'lsi', dims = 2:30)
Monocytes <- FindClusters(object = Monocytes, verbose = FALSE, algorithm = 3,resolution = 0.1)

table(Monocytes$seurat_clusters) 


pdf("Monocytes-cluster.pdf")
DimPlot(Monocytes, reduction = 'umap')
dev.off()
pdf("Monocytes-cluster-sample.pdf")
DimPlot(Monocytes, reduction = "umap", split.by="dataset", ncol=2, label = FALSE)
dev.off()

##marker
DefaultAssay(Monocytes) <- 'RNA'
markers <- FindAllMarkers(Monocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_df = markers %>% group_by(cluster)
head(markers_df)
write.table(markers_df,file="Monocytes-all-markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
head(markers_df)
write.table(markers_df,file="Monocytes-markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)
#EOF


DefaultAssay(integrated) <- 'ATAC'
##call peak
peaks <- CallPeaks(
  object = integrated,
  group.by = "celltype"
)
head(peaks)
peak_new <- as.data.frame(peaks)   ##Convert GRanges to a data frame
write.table(peak_new,file="signac_peak.txt",sep = '\t', quote = FALSE, row.names = FALSE)

##可视化maker或感兴趣基因或某个区域的ATAC-seq信号
DefaultAssay(integrated) <- 'ATAC'
pdf("Fosl1-coverage.pdf")
CoveragePlot(
  object = integrated,
  region = "Fosl1",
  extend.upstream = 5000,
  extend.downstream = 5000
)
dev.off()



##上面生成的peaks是一个GRanges object,所有细胞类型的peak在一起
##combine.peaks =FALSE 每个细胞类型的peaks分别GRanges object，构成一个GRangesList object
peaks <- CallPeaks(
  object = integrated,
  group.by = "celltype",
  outdir="peaks",
  combine.peaks =FALSE,
  extsize = 150,
  shift =75
)
head(peaks[[1]])
for (i in seq_along(peaks)) {
	gr_df <- as.data.frame(peaks[[i]])
	gr_df_flt <- gr_df[grepl("^chr", gr_df$seqnames), ]  ##过滤非chr peaks
    file_name <- paste0(peaks[[i]]$ident[1], ".txt")
    write.table(gr_df_flt, file = file_name,sep="\t",quote = FALSE,  row.names = FALSE)
}



##细胞类型特异性DAR分析
DefaultAssay(integrated) <- 'ATAC'
#cellType.DARs <- FindAllMarkers(integrated, test.use = 'LR', logfc.threshold=0.25, min.pct = 0.05, latent.vars = "nCount_ATAC")
cellType.DARs <- FindAllMarkers(integrated, test.use = 'wilcox', logfc.threshold=0.25, min.pct = 0.05, latent.vars = "nCount_ATAC")
markers_df = cellType.DARs %>% group_by(cluster)
write.table(markers_df,file="cellType.DARs.txt",quote=F,sep="\t",row.names=F,col.names=T)

##DAR可视化
#每个聚类前10个DAR热图(如果小于10，则绘制所有标记)
top10 <- cellType.DARs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("top10DAR-heatmap.pdf")
DoHeatmap(integrated, features = top10$gene,slot="data") + NoLegend()
dev.off()

##可视化每个细胞类型的top1 DAR
top1 <- cellType.DARs %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
pdf("Adipocytes_da_peaks.pdf")
CoveragePlot(
  object = integrated,
  region = top1$gene[13],
  extend.upstream = 10000,
  extend.downstream = 10000
)
dev.off()


pdf("0-Tmem247.pdf")
p <- CoveragePlot(
  object = integrated,
  region = "Tmem247",
  extend.upstream = 5000,
  extend.downstream = 5000
)
p & scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
dev.off()

##每种细胞类型的motif
filtered_da_peaks <- cellType.DARs[grepl("^chr", row.names(cellType.DARs)), ]
Monocytes.da.peak <- filtered_da_peaks[filtered_da_peaks$cluster=="Monocytes", ]

##须先添加motif矩阵
enriched.motifs <- FindMotifs(
  object = integrated,
  features = head(rownames(Monocytes.da.peak),1000)
)

enriched.motifs <- enriched.motifs[order(enriched.motifs[, 7], -enriched.motifs[, 6]), ]
head(enriched.motifs)
write.table(enriched.motifs,file="Monocytes.enriched.motifs",sep="\t",quote = FALSE, row.names = FALSE)

pdf("Monocytes-motif.pdf")
MotifPlot(
  object = integrated,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

##Motif footprinting
# gather the footprinting information for sets of motifs
library(motifmatchr)
integrated <- Footprint(
  object = integrated,
  motif.name = c("KLF15", "BHLHA15", "SPIB","GATA5","SOX8","MEF2B","TEAD3","Ebf2","RUNX1","EHF","POU5F1B","FOXI1","PPARA::RXRA"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)

integrated <- Footprint(
  object = integrated,
  motif.name = c("KLF15", "Tcf12", "SPIB","GATA5","SOX8","MEF2B","TEAD3","Ebf2","RUNX1","EHF","POU5F1B","FOXI1"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# plot the footprint data for each group of cells
pdf("motif-footprint.pdf")
p2 <- PlotFootprint(integrated, features = c("KLF15", "Tcf12", "SPIB","GATA5","SOX8","MEF2B","TEAD3","Ebf2","RUNX1","EHF","POU5F1B","FOXI1"))
p2 + patchwork::plot_layout(ncol = 1)
dev.off()

pdf("KLF15-footprint-1.pdf")
p<-PlotFootprint(integrated, features = c("KLF15"))
p & scale_colour_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
dev.off()




#两种细胞类型差异可及性分析及motif分析
# Find differentially accessible peaks between cell types
# change back to working with peaks instead of gene activities
DefaultAssay(integrated) <- 'ATAC'
da_peaks <- FindMarkers(
  object = integrated,
  ident.1 = "Monocytes",
  ident.2 = "Macrophages",
  test.use = 'LR',
  min.pct = 0.05,   
  latent.vars = 'nCount_ATAC'
)
head(da_peaks)


#Filter rows where row names start with 'chr'
##da_peaks中含有非chr开头的peak，如JH584304.1-9565-10609，在使用 FindMotifs()时会报错，故过滤。
filtered_da_peaks <- da_peaks[grepl("^chr", row.names(da_peaks)), ]
write.table(filtered_da_peaks,file="Monocytes-Macrophages-da_peaks.txt",sep="\t",quote = FALSE)
# get top differentially accessible peaks，没结果？
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
top.da.peak <- rownames(filtered_da_peaks[filtered_da_peaks$p_val < 0.005, ])

pdf("Monocytes-Macrophages-da_peaks-VlnPlot.pdf")
VlnPlot(
  object = integrated,
  features = rownames(filtered_da_peaks)[1],
  pt.size = 0,
  idents = c("Monocytes","Macrophages"),
  cols=c("#1F77B4FF","#2CA02CFF")
)
dev.off()
pdf("Monocytes-Macrophages-da_peaks-FeaturePlot.pdf")
FeaturePlot(
  object = integrated,
  features = rownames(filtered_da_peaks)[1],
  pt.size = 0.1
)
dev.off()

open_Monocytes <- rownames(filtered_da_peaks[filtered_da_peaks$avg_log2FC > 2, ])
open_Macrophages <- rownames(filtered_da_peaks[filtered_da_peaks$avg_log2FC < -2, ])
##关联peak最近的基因
closest_genes_Monocytes <- ClosestFeature(integrated, regions = open_Monocytes)
closest_genes_Macrophages <- ClosestFeature(integrated, regions = open_Macrophages)
head(closest_genes_Monocytes)
head(closest_genes_Macrophages)
write.table(closest_genes_Macrophages,file="closest_genes_Macrophages.txt",sep = '\t', quote = FALSE, row.names = FALSE)

pdf("da_peaks-1.pdf")
CoveragePlot(
  object = integrated,
  region = rownames(filtered_da_peaks)[1],
  extend.upstream = 10000,
  extend.downstream = 10000
)+scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF","#AEC7E8FF","#FFBB78FF","#98DF8AFF"))
dev.off()


##差异可及性区域的motif分析
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

#Adding motif information to the Seurat object
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species="Mus musculus",tax_group = 'vertebrates', all_versions = FALSE)
)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
print(pfm)
# add motif information
integrated <- AddMotifs(
  object = integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# 使用FindMotifs函数进行motif富集分析
enriched.motifs <- FindMotifs(
  object = integrated,
  features = head(rownames(filtered_da_peaks), 1000)
)

# 查看motif富集分析的结果
# sort by p-value and fold change
enriched.motifs <- enriched.motifs[order(enriched.motifs[, 7], -enriched.motifs[, 6]), ]
head(enriched.motifs)

#使用MotifPlot函数绘制motif的位置权重矩阵，可视化不同的motif序列
library(ggseqlogo)
pdf("Monocytes-Macrophages-da_peaks-motif.pdf")
MotifPlot(
  object = integrated,
  motifs = head(rownames(enriched.motifs))
)
dev.off()



##method2:chromVAR
#运行chromVAR计算每个细胞的基序活性得分(motif activity score)，这样我们可以查看每个细胞的motif activities，并且还提供了一种识别细胞类型之间差异激活的基序的替代方法。ChromVAR可识别与细胞之间染色质可及性变异相关的基序
# 使用RunChromVAR函数计算所有细胞中的motif activities
integrated <- RunChromVAR(
  object = integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

##<<EOF
##注意，一下报错和警告是由于chromosome names in the peak matrix that are not present in the BSgenome object, so the sequence for those regions of the genome cannot be retrieved. 
#Error in .getOneSeqFromBSgenomeMultipleSequences(x, name, start, NA, width,  : 
#  sequence GL456210.1 not found
#In addition: Warning message:
#In .merge_two_Seqinfo_objects(x, y) :
#  Each of the 2 combined objects has sequence levels not in the other:
#  - in 'x': chrM, chr1_GL455991_alt,
##EOF

#The best solution is for you to remove these sequences before creating your Seurat object, eg:
features.keep <- as.character(seqnames(granges(integrated))) %in% standardChromosomes(granges(object))
integrated <- integrated[features.keep, ] # if you have multiple assays you'll need to adjust this to keep features from the different assays

DefaultAssay(integrated) <- 'chromvar'

# look at the activity of top1 motif SPIB, the top hit from the overrepresentation testing
FeaturePlot(
  object = integrated,
  features = rownames(enriched.motifs)[[1]],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

##直接测试不同细胞类型之间motif的差异活性得分，这和对不同细胞类型之间的差异可及性peaks进行富集测试的结果相类似。
differential.activity <- FindMarkers(
  object = integrated,
  ident.1 = 'Monocytes',
  ident.2 = 'Macrophages',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

# 使用MotifPlot函数对富集到的motif进行可视化
MotifPlot(
  object = integrated,
  motifs = head(rownames(differential.activity)),
  assay = 'ATAC'
)

