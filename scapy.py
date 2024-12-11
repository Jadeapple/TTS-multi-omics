source activate scanpy

import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as an
#import scanorama
#import scipy.sparse as sp
#from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


##多样本整合分析
##
sample = pd.read_csv('sample.csv', sep =',')
sample

samples = {}
for i in range(sample.shape[0]):
    key = sample.iloc[i, 0]
    samples[key] = sample.iloc[i, 1]

samples    


adatas = {}
for sample_id,filename in samples.items():
    sample_adata = sc.read_visium(filename)
    sample_adata.var_names_make_unique()
    sample_adata.var["mt"] = sample_adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(sample_adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adatas[sample_id] = sample_adata
    sample_adata.uns['spatial'][sample_id] = sample_adata.uns['spatial'][list(sample_adata.uns['spatial'].keys())[0]]
    del sample_adata.uns['spatial'][list(sample_adata.uns['spatial'].keys())[0]]
    
    
#＃合并
adata_spatial = sc.concat(adatas, label= "sample", uns_merge="unique")
adata_spatial.obs_names_make_unique()

##输出原始count数据，cell2location需要原始count矩阵
adata_spatial.write("adata_spatial_original.h5ad")

##filter data
#sc.pp.filter_cells(adata, min_genes=100)
#sc.pp.filter_genes(adata_spatial, min_cells=3)
# Normalize data
sc.pp.normalize_total(adata_spatial)
# Log transformation
sc.pp.log1p(adata_spatial)
# Store raw data
#adata_spatial.raw = adata_spatial
# Extract top highly variable genes
#sc.pp.highly_variable_genes(adata_spatial,min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata_spatial, flavor="seurat", n_top_genes=2000)
##取高变基因
##adata_spatial=adata_spatial[:,adata_spatial.var.highly_variable]
#scale the data
sc.pp.scale(adata_spatial, max_value=10)
# Run dimensionality reduction
sc.pp.pca(adata_spatial)

##harmony整合
###Run integration with harmony
sce.pp.harmony_integrate(adata_spatial, 'sample')

sc.pp.neighbors(adata_spatial, n_pcs =30, use_rep = 'X_pca_harmony')
sc.tl.umap(adata_spatial)
sc.tl.leiden(adata_spatial, resolution=0.5, flavor="igraph", directed=False, n_iterations=2)
#sc.tl.leiden(adata_spatial, resolution=0.7,key_added="clusters") 


##修改leiden cluster name为Clus.0 Clust.1等
##adata_spatial.obs['leiden'] = ['Clust. {0}'.format(i) for i in adata_spatial.obs['leiden']]

##umap绘图
sc.pl.umap(adata_spatial, color=["sample", "leiden"])
plt.savefig('sample_cluster_umap.pdf')
###
##循环绘图
width = 4 * len(list(samples.keys()))
fig, axes = plt.subplots(1, len(list(samples.keys ())), figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="leiden",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")

plt.savefig('sample-cluster-spatial.pdf')



##保存数据
adata_spatial.write("adata_spatial.h5ad")
adata_spatial = sc.read_h5ad('adata_spatial.h5ad')


#查看某些基因在芯片上的表达
sc.pl.spatial(adata, img_key="hires", color=["Vim", "Mgp"], alpha=0.7)
plt.savefig("B468.marker.Vim.pdf")

#width = 4 * len(list(samples.keys()))
width = 4 * 8
fig, axes = plt.subplots(1, 8, figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="Fosl1",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")
plt.savefig('Fosl1-spatial.pdf')

##vmin,vmax设置colorbar范围
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="Cxcl3",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin",vmax=10)



width = 4 * 8
fig, axes = plt.subplots(1, 8, figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="Nrxn3",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")
plt.savefig('Nrxn3-spatial.pdf')


##groups=["niche_8"]指定可视化某个cluster，##循环注意缩进
niche=['niche_1', 'niche_2', 'niche_3', 'niche_4', 'niche_5', 'niche_6', 'niche_7','niche_8']
for i in niche:
	fig, axes = plt.subplots(1, 8, figsize=(width, 3))
	for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
		sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="ct_niche", groups=[i],size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")

    plt.savefig(i+'-spatial.pdf')

##基因集打分
genelist=['Cacna1d','Cntn5','Gja5','Irx3','Kcnj3','Hcn1','Hcn4','Myl4','Tbx3']

##从文件读取基因列表
genes = pd.read_csv('/mnt/sda/jwz/genome/mouse/mm39/collagen.txt',header=0)
print(genes)
genelist=genes['name'].values.tolist()

sc.tl.score_genes(adata_spatial,genelist)
fig, axes = plt.subplots(1, 8, figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="score",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")
	
plt.savefig('CCS-gene-score-spatial.pdf')	




##总共获得12个分子生态位
##添加分子生态位绘图
adata_spatial.obs['niche']=adata_spatial.obs['leiden']
#niche=["niche_0","niche_1","niche_2","niche_3","niche_4","niche_5","niche_6","niche_7","niche_8","niche_9","niche_10","niche_11"]
#adata_spatial.rename_categories('niche',niche)
adata_spatial.obs['niche'] = ['niche{0}'.format(i) for i in adata_spatial.obs['leiden']]
adata_spatial.obs['niche']

##添加分子生态位绘图
sc.pl.umap(adata_spatial, color=["sample", "niche"])
plt.savefig('sample_niche_umap.pdf')

width = 4 * len(list(samples.keys()))  ##len(list(samples.keys()))为样本数目
fig, axes = plt.subplots(1, len(list(samples.keys ())), figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="niche",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")

plt.savefig('sample-niche-spatial.pdf')


##添加细胞类型生态位绘图
a=pd.read_csv("/mnt/sda/jwz/data/TTS/spatial/cell2location/celltype_niche/k_50/niche.txt",sep='\t', header=0, index_col=None)
a1=np.array(a["ct_niche"])
adata_spatial.obs['ct_niche']=pd.Categorical(a1)
adata_spatial.obs['ct_niche']


width = 4 * len(list(samples.keys()))  ##len(list(samples.keys()))为样本数目
fig, axes = plt.subplots(1, len(list(samples.keys ())), figsize=(width, 3))
for x,y in zip(axes, sorted(adata_spatial.uns['spatial'].keys() )):
    sc.pl.spatial(adata_spatial[adata_spatial.obs['sample']==y], frameon=False,color="ct_niche",size=1.3, library_id=y, title=y,ax=x,show=False, legend_loc="right margin")

plt.savefig('sample-ct_niche-spatial.pdf')

adata_spatial.write("adata_spatial_niche.h5ad")



##PROGENy通路活性分析
source activate decoupler
import decoupler as dc
adata_spatial = sc.read_h5ad('/mnt/sda/jwz/data/TTS/spatial/spaceranger/scanpy/adata_spatial_niche.h5ad')
model_500=pd.read_table("/mnt/sda/jwz/data/TTS/spatial/cell2location/celltype_niche/k_50/progeny_mouse_model_500.txt",index_col=None)
model_500

dc.run_mlm(
    mat=adata_spatial,
    net=model_500,
    source='pathway',
    target='gene',
    weight='weight',
    verbose=True,
    use_raw=False
)

# Store in new obsm keys
adata_spatial.obsm['progeny_mlm_estimate'] = adata_spatial.obsm['mlm_estimate'].copy()
adata_spatial.obsm['progeny_mlm_pvals'] = adata_spatial.obsm['mlm_pvals'].copy()
adata_spatial.obsm['progeny_mlm_estimate']

##visualization
acts = dc.get_acts(adata_spatial, obsm_key='progeny_mlm_estimate')
acts

##visualize which pathways are more active in each niche
sc.pl.matrixplot(acts, var_names=acts.var_names, groupby='niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis')
 
plt.savefig("niche-progeny.pdf") 



##寻找每个细胞niche的marker基因
sc.tl.rank_genes_groups(adata_spatial, 'ct_niche', method='t-test')

##查看某个niche的marker
sc.get.rank_genes_groups_df(adata_spatial,group="niche_8")

#可视化top marker
sc.pl.rank_genes_groups(adata_spatial, n_genes=25, sharey=False)
plt.savefig("ct_niche-marker25.pdf") 

sc.pl.rank_genes_groups_dotplot(adata_spatial, n_genes=4)
plt.savefig("ct_niche-marker-dotplot.pdf") 

sc.pl.rank_genes_groups_matrixplot(adata_spatial, n_genes=4, standard_scale='var', cmap='Blues')
#sc.pl.rank_genes_groups_matrixplot(adata_spatial, n_genes=4, standard_scale='var', cmap='bwr')
plt.savefig("ct_niche-marker-matrixplot.pdf") 

sc.pl.violin(adata_spatial, ['Vim', 'Eef1a1', 'Sparc'], groupby='ct_niche')
plt.savefig("ct_niche8-marker-violin.pdf")


##输出
df  = sc.get.rank_genes_groups_df(adata_spatial, group="niche_1")
df = df.sort_values(by="scores", ascending=False)
df.to_csv('niche_1_marker.txt', index=0, sep='\t')

group=['niche_1','niche_2','niche_3','niche_4','niche_5','niche_6','niche_7','niche_8']
for i in group:
    df  = sc.get.rank_genes_groups_df(adata_spatial, group=i)
    df = df.sort_values(by="scores", ascending=False)
    df.to_csv(i+'_marker.txt', index=0, sep='\t')




#保存
adata_spatial.write("adata_spatial.h5ad")




