##单细胞和空转联合分析

##cell2location
#线粒体编码基因与空间定位无关,去除
# find mitochondria-encoded (MT) genes
adata_spatial.var['MT_gene'] = [gene.startswith('mt-') for gene in adata_spatial.var.index]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_spatial.obsm['MT'] = adata_spatial[:, adata_spatial.var['MT_gene'].values].X.toarray()
adata_spatial = adata_spatial[:, ~adata_spatial.var['MT_gene'].values]
adata_spatial

##目录操作
import os
os.getcwd()
os.chdir("/mnt/sda/jwz/data/TTS/spatial/cell2location")


#2.1加载单细胞数据
adata_sc = sc.read_text("/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/signac/new/integrated_counts.csv",delimiter='\t',first_column_names=True)
adata_sc.X = sp.csr_matrix(adata_sc.X)
metadata = pd.read_csv("/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/signac/new/integrated_metadata.csv",sep='\t', header=0, index_col=0)
adata_sc.obs = metadata
r_embedding = pd.read_csv("/mnt/sda/jwz/data/TTS/spatial/scATAC-seq/signac/new/integrated_embedding.csv",sep='\t', header=0, index_col=0)
adata_sc.obsm["X_umap"] = r_embedding.values
adata_sc.write("adata_sc.h5ad")
#读入
adata_sc = sc.read_h5ad('adata_sc.h5ad')
##umap绘图
sc.pl.umap(adata_sc, color='celltype',save="adata_sc-umap.pdf")

#2.3 基因筛选
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
plt.savefig("adata_sc_filter_gene.pdf")
# filter the object
adata_sc = adata_sc[:, selected].copy()

#2.4 参考细胞类型特征估计(NB回归)
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_sc,
                         batch_key='dataset', # 10X reaction / sample / batch
                         labels_key='celltype' # cell type, covariate used for constructing signatures
                         #categorical_covariate_keys=['Method'] # multiplicative technical effects (platform, 3' vs 5', donor effect)
                       )
# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_sc)

# view anndata_setup as a sanity check
mod.view_anndata_setup()	

#训练模型以估计参考细胞类型的特征。
mod.train(max_epochs=250)

#确定模型是否需要更多训练。
#在这里，我们绘制训练过程中的 ELBO 损失历史，从图中去掉前 20 个 epochs。这个图应该呈下降趋势，并在训练结束时趋于平稳。如果仍在下降，请增加 max_epochs。
#没有曲线，不知何故？
mod.plot_history(20)
plt.savefig("mod.plot_history.pdf")

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_sc = mod.export_posterior(
    adata_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

#可选用户可以直接计算后验分布的5%、50%和95%分位数，而不是使用1000个样本（或其他任何分位数）从分布中抽取。这样可以加快在大型数据集上的应用并减少内存需求，但无法用这种方式计算后验均值和标准差。
adata_sc = mod.export_posterior(
    adata_sc, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'use_gpu': True}
)

##QC plot
mod.plot_QC()
plt.savefig("mod.plot_QC1.pdf")
# Save model
mod.save("./mod.train5000/", overwrite=True)

# Save anndata object with results
adata_sc.write("./sc_cell2location.h5ad")	

#模型和输出的h5ad，可以像这样加载:
adata_sc = sc.read_h5ad("./sc_cell2location.h5ad")
mod = cell2location.models.RegressionModel.load("./mod.train5000/", adata_sc)


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
    inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
##good
inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_sc.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


#后续Spatial Mapping仅需要inf_aver文件作为输入，在此保存和读取
inf_aver.to_csv('inf_aver.csv',index=True)
inf_aver = pd.read_csv('inf_aver.csv',index_col=0)


#Cell2location: spatial mapping
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_spatial.var_names, inf_aver.index)
adata_spatial = adata_spatial[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model,need unnormalized count data
cell2location.models.Cell2location.setup_anndata(adata=adata_spatial, batch_key="library_id")
		   
# create and train the model
mod = cell2location.models.Cell2location(
    adata_spatial, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

#Training cell2location
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          #use_gpu=True,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.savefig("full_data_training.pdf")  
plt.legend(labels=['full data training'])

#Exporting estimated posterior distributions of cell abundance and saving results
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_spatial = mod.export_posterior(
    adata_spatial, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save("./spatial_cell2location", overwrite=True)

##QC
mod.plot_QC()
plt.savefig("cell2location_qc.pdf") 

#导出预测的细胞丰度结果为csv文件：作者建议使用后验分布的5%分位数的结果，表示模型具有高置信度的细胞丰度值(即“至少存在这个数量”)。结果在adata_st.obsm下
pd.DataFrame(adata_spatial.obsm['q05_cell_abundance_w_sf']).to_csv("st_cell2location_ratio.csv")
pd.DataFrame(adata_spatial.obsm['q05_cell_abundance_w_sf']).to_csv("st_cell2location_ratio.txt",sep="\t")

# Save anndata object with results
adata_spatial.write("adata_spatial_cell2location.h5ad")

#The model and output h5ad can be loaded later like this:
adata_spatial= sc.read_h5ad("adata_spatial_cell2location.h5ad")
mod = cell2location.models.Cell2location.load("./spatial_cell2location", adata_spatial)

##Visualising cell abundance in spatial coordinates
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_spatial.obs[adata_spatial.uns['mod']['factor_names']] = adata_spatial.obsm['q05_cell_abundance_w_sf']

# select one slide
from cell2location.utils import select_slide
adata_spatial.obs['sample']
slide = select_slide(adata_spatial, 'ISO3d-2')   ##根据sample name选择
# plot in spatial coordinates
sc.pl.spatial(slide, cmap='viridis',
                  # show first 8 cell types
                  color=['Monocytes', 'Macrophages', 'Cardiomyocytes', 'Arterial endothelial cells'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
				 

sc.pl.spatial(slide, cmap='viridis',
                  # selected cell types
                  color=['Adipocytes', 'Arterial endothelial cells', 'Astrocytes', 'B cells', 'Cardiomyocytes', 'Epithelial cell', 'Fibroblasts', 'Macrophages', 'Monocytes', 'Neutrophils', 'Pericytes', 'T cells', 'Venous endothelial cells'],
                  ncols=5, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2',
                  show=True
                 )				 
plt.savefig("ISO3d-2-cell2location.pdf") 


##画单个细胞的分布随时间变化
fig, axes = plt.subplots(1, 8, figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="Cardiomyocytes",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin",vmax=25)

plt.savefig('Cardiomyocytes-spatial.pdf')




##用Leiden聚类识别离散组织区域
##通过使用由cell2location估计的细胞丰度对位置进行聚类，我们识别在细胞组成方面存在差异的组织区域。		
# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_spatial, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 10)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_spatial, resolution=0.3)

# add region as categorical variable
adata_spatial.obs["region_cluster"] = adata_spatial.obs["leiden"].astype("category")		 

## compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_spatial, min_dist = 0.3, spread = 1)

# show regions in UMAP coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_spatial, color=['region_cluster'], size=30,
               ncols = 2, legend_loc='on data',
               legend_fontsize=15)

plt.savefig("cell2location_region_cluster.pdf")

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_spatial, color=['sample'], size=30,
               ncols = 2)
plt.savefig("cell2location_sample.pdf")			   
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata_spatial, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5)