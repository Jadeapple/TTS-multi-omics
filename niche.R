##细胞类型生态位分析

library(Seurat)
library(compositions)
library (tidyverse)
library (clustree)
library (uwot)
library (scran)
library(cluster)
library (RColorBrewer)


#cd /mnt/sda/jwz/data/TTS/spatial/cell2location/celltype_niche
#sed 's/q05cell_abundance_w_sf_//g' st_cell2location_ratio.txt >st_cell2location_celltype.txt

##读入spot X cell2location celltype文件
integrated_compositions=read.table("st_cell2location_celltype.txt",head=T,sep="\t",row.names=1)
head(integrated_compositions)
integrated_compositions[1:5,1:5]
max(integrated_compositions)
min(integrated_compositions)
##ILR转换
baseILR <- ilrBase(x = integrated_compositions , method = "basic")
cell_ilr <- as.matrix(ilr(integrated_compositions, baseILR)) 
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

# Make community graph
k_vect <- c(10,20,30,40,50)
k_vect <- set_names(k_vect, paste0("k_", k_vect))
cluster_info <- map(k_vect, function(k) {
print(k)
print("Generating SNIN")
snn_graph <- scran::buildSNNGraph(x = t(cell_ilr %>% as.data.frame() %>% as.matrix()),k=k)
print ("Louvain clustering")
clust.louvain <- igraph::cluster_louvain(snn_graph)
clusters <- tibble(cluster = clust.louvain$membership, spot_id = rownames(cell_ilr))
})

cluster_info <- cluster_info %>% enframe() %>% unnest() %>% pivot_wider(names_from = name, values_from = cluster)

#查看cluster数目
table(cluster_info$k_50)


# Create UMAP and plot the compositlons
#n_neighbors 数值与k_50数值保持一致
comp_umap <- umap(cell_ilr, n_neighbors = 50, n_epochs = 1000,metric = "cosine") %>% as.data.frame() %>% mutate(row_id= rownames(cell_ilr) )
head(comp_umap)
comp_umap <- comp_umap[, 1:3]
comp_umap <- comp_umap %>% left_join(cluster_info, by = c("row_id" = "spot_id") )
head(comp_umap)
dim(comp_umap)
dim(cluster_info)


# Make the niche annotation meta-data
cluster_info <- comp_umap %>% dplyr::select(c("row_id", "k_50") ) %>% dplyr::rename("niche" = k_50) %>% dplyr::mutate(ct_niche = paste0("niche_", niche))
head(cluster_info)
write.table(cluster_info,file="cluster_info.txt",sep="\t",quote = F, row.names = F)


##合并comp_umap和cluster_info，使用ggplot2绘制umap图
umap <- merge(comp_umap,cluster_info,by="row_id")
mycolor=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF")
p <- ggplot(umap,aes(x= V1 , y = V2 ,color = ct_niche)) +  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = mycolor) +
	theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white")) +
	theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=5)))
	##添加箭头
	#geom_segment(aes(x = min(umap$V1) , y = min(umap$V2) ,
    #   xend = min(umap$V1) +3, yend = min(umap$V2) ),
    #    colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
    #geom_segment(aes(x = min(umap$V1)  , y = min(umap$V2)  ,
    #    xend = min(umap$V1) , yend = min(umap$V2) + 3),
    #    colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
    #annotate("text", x = min(umap$V1) +1.5, y = min(umap$V2) -1, label = "UMAP1",
    #    color="black",size = 3, fontface="bold" ) + 
    #annotate("text", x = min(umap$V1) -1, y = min(umap$V2) + 1.5, label = "UMAP2",
    #    color="black",size = 3, fontface="bold" ,angle=90) 
pdf("niche-umap.pdf")
p
dev.off()

##保存umap信息
write.table(umap,file="niche_umap.txt",sep="\t",quote = F, row.names = F)


# What are the cells that define the niches?
niche_summary_pat <- integrated_compositions %>% as.data.frame() %>% rownames_to_column("row_id") %>% pivot_longer(-row_id, values_to =
"ct_prop", names_to = "cell_type") %>% left_join(cluster_info) %>% mutate(orig.ident = strsplit(row_id,"[..]") %>% map_chr(., ~ .x[1])) %>% group_by(orig.ident,ct_niche,cell_type) %>% summarize(median_ct_prop=median(ct_prop))

niche_summary <- niche_summary_pat %>% ungroup() %>% group_by(ct_niche, cell_type) %>% summarise(patient_median_ct_prop = median(median_ct_prop))

write.table(niche_summary_pat, file ="./niche_summary_pat.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# Show the compositions of cells of each niche 
pdf("niche_summary_pat.pdf")
niche_summary_pat %>%
  ggplot(aes(x = ct_niche, y  = median_ct_prop)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(.~cell_type, ncol = 3,scales = "free_y")
dev.off()



# Data manipulation to have clustered data
niche_summary_mat <- niche_summary %>% pivot_wider(values_from = patient_median_ct_prop,names_from=cell_type,values_fill=0) %>% column_to_rownames("ct_niche") %>% as.matrix()
niche_order <- hclust(dist(niche_summary_mat))
niche_order <- niche_order$labels[niche_order$order]
ct_order <- hclust(dist(t(niche_summary_mat)))
ct_order <- ct_order$labels[ct_order$order]

# Find characteristic cell types of each niche
# We have per patient the proportion of each cell-type in each niche
# We have per patient the proportion of each cell-type in each niche
run_wilcox_up <- function(prop_data) {
prop_data_group <- prop_data[["ct_niche"]] %>% unique() %>% set_names()
	map(prop_data_group, function(g) {
		test_data <- prop_data %>% mutate(test_group = ifelse(ct_niche == g, "target", "rest")) %>% mutate( test_group = factor(test_group, levels=c("target", "rest")))
		wilcox.test(median_ct_prop ~ test_group, data = test_data, alternative = "greater") %>% broom::tidy()
}) %>% enframe("ct_niche") %>% unnest()
}

wilcoxon_res <- niche_summary_pat %>%
ungroup() %>%
group_by(cell_type) %>% 
nest() %>%
mutate(wres = map(data, run_wilcox_up)) %>%
dplyr::select(wres) %>% 
unnest() %>% 
ungroup() %>%
mutate(p_corr = p.adjust(p.value) )

wilcoxon_res <- wilcoxon_res %>% mutate(significant = ifelse(p_corr <= 0.05, "*", ""))

write.table(wilcoxon_res, file ="./wilcoxon_res_cell_niches.txt", col.names = T, row.names = F, quote = F, sep = "\t")

##niche细胞组成热图
mean_ct_prop_plt <- niche_summary %>%
left_join(wilcoxon_res, by = c("ct_niche","cell_type")) %>%
mutate(cell_type = factor(cell_type, levels = ct_order),
	ct_niche = factor(ct_niche, levels = niche_order)) %>%
ungroup() %>%
group_by(cell_type) %>%
mutate(scaled_comp=(patient_median_ct_prop - mean(patient_median_ct_prop))/sd(patient_median_ct_prop)) %>%
ungroup() %>%
ggplot(aes (x = ct_niche, y = cell_type, fill = scaled_comp) )+
geom_tile() +
geom_text(aes(label = significant)) +
theme_bw()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
	legend.position = "bottom",
	plot.margin = unit(c(1, 1,1, 1), "cm"),
	axis.text.y = element_text(size=12))+
scale_fill_gradient2(high ="red", mid ='white', low = "blue")
pdf("niches_celltype.pdf")
mean_ct_prop_plt
dev.off()


##每个样本niches占比barplot
awk 'NR>1' cluster_info.txt | sed 's/-/\t/g' | awk '{print $3"-"$4"\t"$6}' | sed '1i sample\tniche' >sample_cluster_info.txt
data <- read.table("sample_cluster_info.txt",head=T)
table(data$niche,data$sample)
nicheratio <- prop.table(table(data$niche,data$sample), margin = 2)
nicheratio
nicheratio <- as.data.frame(nicheratio)
colnames(nicheratio)=c("Niche","Sample","Freq")
library(ggplot2)
pdf("sample-mniche-ratio.pdf")
ggplot(nicheratio) + geom_bar(aes(x =Sample, y= Freq, fill = Niche),stat = "identity",width = 0.7,size = 0.5)+ 
  scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF"))+
  theme_classic() +
  labs(x='Sample',y = 'Niche Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()



##EOF
# Finally. describe the proportions of those niches in alt the data
cluster_counts <- cluster_info %>%
dplyr::select_at(c("row_id", "ct_niche")) %>%
group_by(ct_niche) %>%
summarise(nspots = length(ct_niche)) %>%
mutate(prop_spots = nspots/sum(nspots))

write.table(cluster_counts, file = "/niche_prop_summary.txt",quote = F, sep = "\t")

##绘图
barplts <- cluster_counts %>%
mutate(ct_niche = factor(ct_niche, levels = niche_order)) %>%
ggplot (aes (y = ct_niche, x = prop_spots) )+
geom_bar (stat = "identity")+
theme_classic() + ylab ("") +
theme(axis.text.y = element_blank() ,
	plot.margin = unit(c(0, 0,0,0),"cm")
	axis.text.x = element_text(size=12))
	
niche_summary_plt <- cowplot::plot_grid(mean_ct_prop_plt, barplts, align = "hv", axis = "tb")
pdf("characteristic_ct_niches.pdf", width = 12, height = 7)
plot(niche_summary_plt)
dev.off()
##EOF


##niche空间可视化，python,见scanpy.py
##将niche添加到adata.obs
adata_spatial= sc.read_h5ad("/mnt/sda/jwz/data/TTS/spatial/cell2location/adata_spatial_cell2location.h5ad")
#awk '{print $3}' cluster_info.txt >niche.txt
a=pd.read_csv("niche.txt",sep='\t', header=0, index_col=None)
a1=np.array(a["ct_niche"])
adata_spatial.obs['niche']=pd.Categorical(a1)
adata_spatial.obs['niche']

##选择芯片可视化
from cell2location.utils import select_slide
slide = select_slide(adata_spatial, 'ISO1d-2')
sc.pl.spatial(slide, cmap='viridis',
                  # selected cell types
                  color="niche",
                  ncols=5, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  #vmin=0, vmax='p99.2',
                  show=True
                 )	
plt.savefig("ISO1d-2-niche.pdf") 

##循环绘图
sample=adata_spatial.obs['sample'].unique().tolist()
fig, axes = plt.subplots(2, 4, figsize=(18, 12), tight_layout=True,dpi=600) # nrows和ncols取决于想要的画图个数
x=0;y=0
for i in sample: #样本循环
	slide = select_slide(adata_spatial, i)
	fig = sc.pl.spatial(slide, cmap='viridis',color="niche",size=1.3,img_key='hires', ax=axes[y][x],show=True)
	x = x +1 if x <3 else 0
	y = y +1 if y <1 and x ==0 else y

plt.savefig('sample-niche-spatial.pdf')

adata_spatial.write("adata_spatial_cell2location_niche.h5ad")



##分子生态位细胞组成分析
sed '1d' barcode_niche.txt | sed '1i row_id\tct_niche' >molecular_cluster_info.txt
cluster_info=read.table("molecular_cluster_info.txt",head=T,sep="\t")

# What are the cells that define the niches?
niche_summary_pat <- integrated_compositions %>% as.data.frame() %>% rownames_to_column("row_id") %>% pivot_longer(-row_id, values_to =
"ct_prop", names_to = "cell_type") %>% left_join(cluster_info) %>% mutate(orig.ident = strsplit(row_id,"[..]") %>% map_chr(., ~ .x[1])) %>% group_by(orig.ident,ct_niche,cell_type) %>% summarize(median_ct_prop=median(ct_prop))

niche_summary <- niche_summary_pat %>% ungroup() %>% group_by(ct_niche, cell_type) %>% summarise(patient_median_ct_prop = median(median_ct_prop))

write.table(niche_summary_pat, file ="./molecular_niche_summary_pat.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# Show the compositions of cells of each niche 
pdf("molecular_niche_summary_pat.pdf")
niche_summary_pat %>%
  ggplot(aes(x = ct_niche, y  = median_ct_prop)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(.~cell_type, ncol = 3,scales = "free_y")
dev.off()



# Data manipulation to have clustered data
niche_summary_mat <- niche_summary %>% pivot_wider(values_from = patient_median_ct_prop,names_from=cell_type,values_fill=0) %>% column_to_rownames("ct_niche") %>% as.matrix()
niche_order <- hclust(dist(niche_summary_mat))
niche_order <- niche_order$labels[niche_order$order]
ct_order <- hclust(dist(t(niche_summary_mat)))
ct_order <- ct_order$labels[ct_order$order]

# Find characteristic cell types of each niche
# We have per patient the proportion of each cell-type in each niche
# We have per patient the proportion of each cell-type in each niche
run_wilcox_up <- function(prop_data) {
prop_data_group <- prop_data[["ct_niche"]] %>% unique() %>% set_names()
	map(prop_data_group, function(g) {
		test_data <- prop_data %>% mutate(test_group = ifelse(ct_niche == g, "target", "rest")) %>% mutate( test_group = factor(test_group, levels=c("target", "rest")))
		wilcox.test(median_ct_prop ~ test_group, data = test_data, alternative = "greater") %>% broom::tidy()
}) %>% enframe("ct_niche") %>% unnest()
}

wilcoxon_res <- niche_summary_pat %>%
ungroup() %>%
group_by(cell_type) %>% 
nest() %>%
mutate(wres = map(data, run_wilcox_up)) %>%
dplyr::select(wres) %>% 
unnest() %>% 
ungroup() %>%
mutate(p_corr = p.adjust(p.value) )

wilcoxon_res <- wilcoxon_res %>% mutate(significant = ifelse(p_corr <= 0.05, "*", ""))

write.table(wilcoxon_res, file ="./wilcoxon_res_molecular_cell_niches.txt", col.names = T, row.names = F, quote = F, sep = "\t")

##niche细胞组成热图
mean_ct_prop_plt <- niche_summary %>%
left_join(wilcoxon_res, by = c("ct_niche","cell_type")) %>%
mutate(cell_type = factor(cell_type, levels = ct_order),
	ct_niche = factor(ct_niche, levels = niche_order)) %>%
ungroup() %>%
group_by(cell_type) %>%
mutate(scaled_comp=(patient_median_ct_prop - mean(patient_median_ct_prop))/sd(patient_median_ct_prop)) %>%
ungroup() %>%
ggplot(aes (x = ct_niche, y = cell_type, fill = scaled_comp) )+
xlab("molecular_niche")+
geom_tile() +
geom_text(aes(label = significant)) +
theme_bw()+
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
	legend.position = "bottom",
	plot.margin = unit(c(1, 1,1, 1), "cm"),
	axis.text.y = element_text(size=12))+
scale_fill_gradient2(high ="red", mid ='white', low = "blue")
pdf("molecular_niches_celltype.pdf")
mean_ct_prop_plt
dev.off()
