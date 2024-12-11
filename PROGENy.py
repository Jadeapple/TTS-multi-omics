##1. PROGENy通路活性分析
source activate decoupler

import scanpy as sc
import decoupler as dc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##plotting options，可选
sc.setting.set_figure_params(dpi=200,frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4,4))
adata = sc.read_h5ad('adata_spatial.h5ad')

progeny = dc.get_progeny(organism='mouse', top=500) ##python报错，获得困难，human可以直接获得

##EOF
##通过R获取mouse top500 significant genes per pathway
##retrieve all genes from each pathway in mouse
library(progeny)
library(dplyr)
model <- progeny::model_mouse_full
head(model)
# Get top 500 significant genes per pathway
model_500 <- model %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 500)
write.table(model_500,file="progeny_mouse_model_500.txt",sep = '\t', quote = FALSE, row.names = FALSE)
##EOF

model_500=pd.read_table("progeny_mouse_model_500.txt",index_col=None)
model_500

adata = sc.read_h5ad('adata_spatial.h5ad')
dc.run_mlm(
    mat=adata,
    net=model_500,
    source='pathway',
    target='gene',
    weight='weight',
    verbose=True,
    use_raw=False
)

# Store in new obsm keys
adata.obsm['progeny_mlm_estimate'] = adata.obsm['mlm_estimate'].copy()
adata.obsm['progeny_mlm_pvals'] = adata.obsm['mlm_pvals'].copy()
adata.obsm['progeny_mlm_estimate']

##visualization
acts = dc.get_acts(adata, obsm_key='progeny_mlm_estimate')
acts

ad = acts[acts.obs.library_id == "ISO7d-2", :].copy()
sc.pl.spatial(
    ad,
    color=['Androgen','EGFR','Estrogen','Hypoxia','JAK-STAT','MAPK','NFkB','p53', 'PI3K', 'TGFb','TNFa','Trail','VEGF','WNT'],
    cmap='viridis',
    ncols=5,
    size=1.3,
    vcenter=0,
    frameon=False,
    library_id="ISO7d-2"
)
plt.savefig("ISO7d-2-progeny.pdf")




##niche通路活性分析
adata_spatial= sc.read_h5ad("/mnt/sda/jwz/data/TTS/spatial/spaceranger/scanpy/old/adata_spatial.h5ad")
##将niche添加到adata.obs
#awk '{print $3}' cluster_info.txt >niche.txt
a=pd.read_csv("/mnt/sda/jwz/data/TTS/spatial/cell2location/celltype_niche/k_50/niche.txt",sep='\t', header=0, index_col=None)
a1=np.array(a["ct_niche"])
adata_spatial.obs['niche']=pd.Categorical(a1)
adata_spatial.obs['niche']


model_500=pd.read_table("progeny_mouse_model_500.txt",index_col=None)
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



##2. 转录因子活性分析 Transcription factor activity inference
#net = dc.get_collectri(organism='human', split_complexes=False)
#pd.DataFrame(net).to_csv("/mnt/sda/jwz/genome/decoupler/TF_net_human.csv",sep="\t") 
net = dc.get_collectri(organism='mouse', split_complexes=False)
pd.DataFrame(net).to_csv("/mnt/sda/jwz/genome/decoupler/TF_net_mouse.csv",sep="\t")   

adata_spatial = sc.read_h5ad('adata_spatial_niche.h5ad')
net=pd.read_csv("/mnt/sda/jwz/genome/decoupler/net_mouse.csv",sep='\t', header=0, index_col=0)
dc.run_ulm(
    mat=adata_spatial,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw=False
)
adata_spatial.obsm['ulm_estimate']

adata_spatial.obsm['collectri_ulm_estimate'] = adata_spatial.obsm['ulm_estimate'].copy()
adata_spatial.obsm['collectri_ulm_pvals'] = adata_spatial.obsm['ulm_pvals'].copy()
adata_spatial

acts = dc.get_acts(adata_spatial, obsm_key='collectri_ulm_estimate')
acts
# identify which are the top TF per cluster. 
##2.1. 分子生态位TF
df = dc.rank_sources_groups(acts, groupby='niche', reference='rest', method='t-test_overestim_var')
df
pd.DataFrame(df).to_csv("molecular_niche_TF_activity.txt",sep="\t") 

##extract the top 3 markers per cluster
n_markers = 3
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
source_markers

#plot the top 3 markers
sc.pl.matrixplot(acts, source_markers, 'niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis')
plt.savefig("molecular_niche_TF_top3.pdf")

#visualize TF on the tissue some marker TFs:
sc.pl.spatial(
    acts,
    library_id='Control-1',
    color=['Jun', 'Smad4', 'Ets1', 'Sp1','Stat3','Atf2','Rela','Fos'],
    cmap='viridis',
    size=1.5,
    vcenter=0,
    frameon=False
)
plt.savefig("Control-1_molecular_niche_TF_spatial.pdf")
##循环可视化
for sample in ['Control-1', 'Control-2', 'ISO1d-1', 'ISO1d-2', 'ISO3d-1', 'ISO3d-2', 'ISO7d-1', 'ISO7d-2']:
    sc.pl.spatial(acts[acts.obs['sample']==sample], frameon=False,color=['Jun', 'Smad4', 'Ets1', 'Sp1','Stat3','Atf2','Rela','Fos'],size=1.3, library_id=sample,cmap='viridis',vcenter=0)
    file=sample+'_molecular_niche_TF_spatial.pdf'
    plt.savefig(file)

for sample in ['Control-1', 'Control-2', 'ISO1d-1', 'ISO1d-2', 'ISO3d-1', 'ISO3d-2', 'ISO7d-1', 'ISO7d-2']:
    sc.pl.spatial(acts[acts.obs['sample']==sample], frameon=False,color=['Hsf4', 'Pitx1', 'Tlx1', 'Hsf2'],size=1.3, library_id=sample,cmap='viridis',vcenter=0)
    file=sample+'_molecular_niche7_TF_spatial.pdf'
    plt.savefig(file)

term=['Jun', 'Smad4', 'Ets1', 'Sp1','Stat3']               
for i in term:
    sc.pl.violin(acts,keys=i,groupby='niche',rotation=90)
    file='molecular_niche_TF_'+i+'.pdf'
    plt.savefig(file) 


##2.2 细胞类型生态位TF    
df = dc.rank_sources_groups(acts, groupby='ct_niche', reference='rest', method='t-test_overestim_var')
df
pd.DataFrame(df).to_csv("ct_niche_TF_activity.txt",sep="\t") 

##extract the top 3 markers per cluster
n_markers = 3
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
source_markers

#plot the top 3 markers
sc.pl.matrixplot(acts, source_markers, 'ct_niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis')
plt.savefig("ct_niche_TF_top3.pdf")

#visualize TF on the tissue some marker TFs:
##循环可视化
for sample in ['Control-1', 'Control-2', 'ISO1d-1', 'ISO1d-2', 'ISO3d-1', 'ISO3d-2', 'ISO7d-1', 'ISO7d-2']:
    sc.pl.spatial(acts[acts.obs['sample']==sample], frameon=False,color=['Jun', 'Sp1', 'Smad3', 'Sp3','Ets1','Stat3','Smad4','Fosl1'],size=1.3, library_id=sample,cmap='viridis',vcenter=0)
    file=sample+'_ct_niche_TF_spatial.pdf'
    plt.savefig(file) 

adata_spatial.write("./TF_activity/adata_spatial_TF_activity.h5ad") 

term=['Jun', 'Sp1', 'Smad3', 'Sp3','Ets1']               
for i in term:
    sc.pl.violin(acts,keys=i,groupby='ct_niche',rotation=90)
    file='ct_niche_TF_'+i+'.pdf'
    plt.savefig(file) 


##3. MSigDB gene sets
#The Molecular Signatures Database (MSigDB) is a resource containing a collection of gene sets annotated to different biological processes.
msigdb = dc.get_resource('MSigDB'，organism='mouse') 
##或者读取人的数据，然后转换成小鼠数据
msigdb = dc.get_resource('MSigDB')  ##人数据
msigdb
pd.DataFrame(msigdb).to_csv("/mnt/sda/jwz/genome/decoupler/msigdb_human.csv",sep="\t")   

# Filter by hallmark
msigdb = msigdb[msigdb['collection']=='hallmark']

# Remove duplicated entries
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
msigdb

# transform the obtained resource into mouse genes. Organisms can be defined by their common name, latin name or NCBI Taxonomy identifier.
# Translate targets
mouse_msigdb = dc.translate_net(msigdb, target_organism = 'mouse', unique_by = ('geneset', 'genesymbol'))
mouse_msigdb

msigdb=mouse_msigdb

# Rename
msigdb.loc[:, 'geneset'] = [name.split('HALLMARK_')[1] for name in msigdb['geneset']]
msigdb

dc.run_ora(
    mat=adata_spatial,
    net=msigdb,
    source='geneset',
    target='genesymbol',
    verbose=True,
    use_raw=False
)

# Store in a different key
adata_spatial.obsm['msigdb_ora_estimate'] = adata_spatial.obsm['ora_estimate'].copy()
adata_spatial.obsm['msigdb_ora_pvals'] = adata_spatial.obsm['ora_pvals'].copy()
##查看
adata_spatial.obsm['msigdb_ora_estimate']
adata_spatial.obsm['msigdb_ora_estimate'].iloc[:, 0:5]

#visualization
#extract the activities from the adata object.
acts = dc.get_acts(adata_spatial, obsm_key='msigdb_ora_estimate')

# We need to remove inf and set them to the maximum value observed
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

acts

#identify which are the top terms per cluster
##3.1分子生态位
df = dc.rank_sources_groups(acts, groupby='niche', reference='rest', method='t-test_overestim_var')
df
pd.DataFrame(df).to_csv("molecular_niche_msigdb.txt",sep="\t") 

n_top = 3
term_markers = df.groupby('group').head(n_top).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
term_markers


sc.pl.matrixplot(acts, term_markers, 'niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis')   
plt.savefig("molecular_niche_msigdb_top3.pdf") 

sc.pl.matrixplot(acts, term_markers, 'niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis',swap_axes=True)   ##swap_axes=True用于旋转坐标轴
plt.savefig("molecular_niche_msigdb_top3_1.pdf") 

##小提琴图
sc.pl.violin(
    acts,
    keys='ANGIOGENESIS',
    groupby='niche',
    rotation=90
) 
plt.savefig("molecular_niche_msigdb_ANGIOGENESIS.pdf")

sc.pl.violin(
    acts,
    keys='EPITHELIAL_MESENCHYMAL_TRANSITION',
    groupby='niche',
    rotation=90
) 
plt.savefig("molecular_niche_msigdb_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf")  

##循环绘图
term=['INFLAMMATORY_RESPONSE','KRAS_SIGNALING_UP','TNFA_SIGNALING_VIA_NFKB']
               
for i in term:
    sc.pl.violin(acts,keys=i,groupby='niche',rotation=90)
    file='molecular_niche_msigdb_'+i+'.pdf'
    plt.savefig(file)                 
    
##空间循环可视化
for sample in ['Control-1', 'Control-2', 'ISO1d-1', 'ISO1d-2', 'ISO3d-1', 'ISO3d-2', 'ISO7d-1', 'ISO7d-2']:
    sc.pl.spatial(acts[acts.obs['sample']==sample], frameon=False,color=['ANGIOGENESIS', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INFLAMMATORY_RESPONSE','KRAS_SIGNALING_UP'],size=1.3, library_id=sample,cmap='viridis',vcenter=0)
    file=sample+'_molecular_niche_msigdb_spatial.pdf'
    plt.savefig(file) 

 

###3.2细胞类型生态位
df = dc.rank_sources_groups(acts, groupby='ct_niche', reference='rest', method='t-test_overestim_var')
df
pd.DataFrame(df).to_csv("ct_niche_msigdb.txt",sep="\t") 

n_top = 3
term_markers = df.groupby('group').head(n_top).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
term_markers


sc.pl.matrixplot(acts, term_markers, 'ct_niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis')   
plt.savefig("ct_niche_msigdb_top3.pdf") 

sc.pl.matrixplot(acts, term_markers, 'ct_niche', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='viridis',swap_axes=True)   ##swap_axes=True用于旋转坐标轴
plt.savefig("ct_niche_msigdb_top3_1.pdf") 

##小提琴图
sc.pl.violin(
    acts,
    keys='ANGIOGENESIS',
    groupby='niche',
    rotation=90
) 
plt.savefig("molecular_niche_msigdb_ANGIOGENESIS.pdf")

sc.pl.violin(
    acts,
    keys='EPITHELIAL_MESENCHYMAL_TRANSITION',
    groupby='niche',
    rotation=90
) 
plt.savefig("molecular_niche_msigdb_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf")  

##循环绘图
term=['EPITHELIAL_MESENCHYMAL_TRANSITION','APICAL_JUNCTION','KRAS_SIGNALING_UP','ANGIOGENESIS','TNFA_SIGNALING_VIA_NFKB']
               
for i in term:
    sc.pl.violin(acts,keys=i,groupby='ct_niche',rotation=90)
    file='ct_niche_msigdb_'+i+'.pdf'
    plt.savefig(file)                 
    
##空间循环可视化
for sample in ['Control-1', 'Control-2', 'ISO1d-1', 'ISO1d-2', 'ISO3d-1', 'ISO3d-2', 'ISO7d-1', 'ISO7d-2']:
    sc.pl.spatial(acts[acts.obs['sample']==sample], frameon=False,color=['EPITHELIAL_MESENCHYMAL_TRANSITION','APICAL_JUNCTION','KRAS_SIGNALING_UP','ANGIOGENESIS','TNFA_SIGNALING_VIA_NFKB'],size=1.3, library_id=sample,cmap='viridis',vcenter=0,ncols=5)
    file=sample+'_ct_niche_msigdb_spatial.pdf'
    plt.savefig(file) 

adata_spatial.write(".adata_spatial_msigdb.h5ad")  