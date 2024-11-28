import numpy as np
import scanpy as sc
import scipy.sparse as sp
import loompy as lp
import matplotlib.pyplot as plt
import pandas as pd
import json
import zlib
import base64
import umap
from MulticoreTSNE import MulticoreTSNE as TSNE

adata = sc.read_h5ad("20241018_AG_symbs.h5ad")

sc.pp.filter_genes(adata, min_cells=5 )

expression_matrix_dense = adata.X.transpose().toarray()
row_attrs = { 
  "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
  "CellID":  np.array(adata.obs.index) ,
  "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
  "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( "20241017_AG_scenic.loom", expression_matrix_dense, row_attrs, col_attrs)

# pySCENIC output
lf = lp.connect( "20241017_AG_pyscenic_output.loom", mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )

# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )

auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
  tmp = x.get('regulon').replace("(","_(")
  x.update( {'regulon': tmp} )

adata = sc.read_h5ad("20241017_AG_symbs.h5ad")  
sc.pp.filter_genes(adata, min_cells=3 )
# simply compute the number of genes per cell (computers 'n_genes' column)
sc.pp.filter_cells(adata, min_genes=0)
# mito and genes/counts cuts
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
adata.obs['percent_mito'] = np.sum(
  adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
# save a copy of the raw data
adata.raw = adata

# Total-count normalize (library-size correct) to 10,000 reads/cell
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

# log transform the data.
sc.pp.log1p(adata)

# identify highly variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

canonical_genes = ["INS", "GHRL", "SST", "PPY", "GCG","CD4","TRAC",
                   "CD8A","CD8B","CD20","CD79A","CD79B","SOX10","CRYAB","KRT19","MUC5B","SPP1","MMP7",
                   "CPA1","REG1A",
                   "TPSAB1","TPSB2","CPA3",
                   "PECAM1","VWF","CD34","CD93","PLVAP",
                   "CD163","CD68","CD14","S100A8","LYZ",
                   "COL1A1","PDGFRB","COL1A2","FN1","RGS5"
]


genes_to_keep_bool = np.logical_or(adata.var.highly_variable, adata.var_names.isin(canonical_genes))

# Subset the AnnData object
adata = adata[:, genes_to_keep_bool]
adata = adata.copy()

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
# compute UMAP
sc.tl.umap(adata)
# tSNE
X_dense = adata.X.toarray()

# Run t-SNE on the dense matrix
tsne = TSNE(n_jobs=16, random_state=42)
adata.obsm['X_tsne'] = tsne.fit_transform(X_dense)

dr_tsne_df = pd.DataFrame(dr_tsne, index=adata.obs.index, columns=['X', 'Y'])
dr_umap_df = pd.DataFrame(dr_umap, index=adata.obs.index, columns=['X', 'Y'])

tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])
Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
  pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
  pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
  dr_tsne_df['X'] ,
  dr_umap_df['X']
], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
  pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
  pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
  dr_tsne_df['Y'] ,
  dr_umap_df['Y']
], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3','4']

metaJson = {}

metaJson['embeddings'] = [
  {
    "id": -1,
    "name": f"Scanpy t-SNE (highly variable genes)"
  },
  {
    "id": 1,
    "name": f"Scanpy UMAP  (highly variable genes)"
  },
  {
    "id": 2,
    "name": "Scanpy PC1/PC2"
  },
  {
    "id": 3,
    "name": "SCENIC AUC t-SNE"
  },
  {
    "id": 4,
    "name": "SCENIC AUC UMAP"
  },
]

metaJson["clusterings"] = [{
  "id": 0,
  "group": "Scanpy",
  "name": "Scanpy louvain default resolution",
  "clusters": [],
}]

metaJson["metrics"] = [
  {
    "name": "nUMI"
  }, {
    "name": "nGene"
  }, {
    "name": "Percent_mito"
  }
]

sc.tl.louvain(adata,resolution=0.4)

metaJson["annotations"] = [
  {
    "name": "Louvain_clusters_Scanpy",
    "values": list(set( adata.obs['louvain'].astype(str) ))
  },
  {
    "name": "celltype",
    "values": list(set(adata.obs['celltype'].values))
  },
   {
    "name": "Donor",
    "values": list(set( adata.obs['Donor'].astype(str) ))
  }, 
]

# SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
  clustDict = {}
  clustDict['id'] = i
  clustDict['description'] = f'Unannotated Cluster {i + 1}'
  metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()
clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)

def dfToNamedMatrix(df):
  arr_ip = [tuple(i) for i in df.values]
  dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
  arr = np.array(arr_ip, dtype=dtyp)
  return arr

col_attrs = {
  "CellID": np.array(adata.obs.index),
  "nUMI": np.array(adata.obs['n_counts'].values),
  "nGene": np.array(adata.obs['n_genes'].values),
  "Louvain_clusters_Scanpy": np.array( adata.obs['louvain'].values ),
  "celltype":  np.array( adata.obs['celltype'].values),
  "Donor": np.array( adata.obs['Donor'].values),
  "Percent_mito": np.array(adata.obs['percent_mito'].values),
  "Embedding": dfToNamedMatrix(tsneDF),
  "Embeddings_X": dfToNamedMatrix(Embeddings_X),
  "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
  "RegulonsAUC": dfToNamedMatrix(auc_mtx),
  "Clusterings": dfToNamedMatrix(clusterings),
  "ClusterID": np.array(adata.obs['louvain'].values)
  }
  
row_attrs = {
  "Gene": lf.ra.Gene,
  "Regulons": regulons,
}

attrs = {
  "title": "sampleTitle",
  "MetaData": json.dumps(metaJson),
  "Genome": 'hg38',
  "SCopeTreeL1": "",
  "SCopeTreeL2": "",
  "SCopeTreeL3": ""
}

# compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

lp.create(
  filename = "20241018_AG_pyscenic_final.loom" ,
  layers=lf[:,:],
  row_attrs=row_attrs, 
  col_attrs=col_attrs, 
  file_attrs=attrs
)
lf.close()

lf = lp.connect( "20241018_AG_pyscenic_final.loom", mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
  regulons[i] =  list(r[r==1].index.values)

# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
  [
    pd.DataFrame( lf.ca.Louvain_clusters_Scanpy, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.celltype, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.Donor, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.Percent_mito, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
  ],
  axis=1
)

cellAnnot.columns = [
  'Louvain_clusters_Scanpy',
  'celltype',
  'Donor',
  'Percent_mito',
  'nGene',
  'nUMI']

# capture embeddings:
dr = [
  pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
]

dr_names = [
  meta['embeddings'][0]['name'].replace(" ","_")
]

# add other embeddings
drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )

for i in range( len(drx.columns) ):
  dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
  dr_names.append( meta['embeddings'][i+1]['name'].replace(" ","_").replace('/','-') )

# rename columns:
for i,x in enumerate( dr ):
  x.columns = ['X','Y']

lf.close()

for i, x in enumerate(dr):
  adata.obsm['X_' + dr_names[i]] = x.values

# Assuming auc_mtx is the matrix of AUC scores (from the AUCell step)
adata.obsm['RegulonsAUC'] = auc_mtx

adata.write_h5ad("20241018_adata_AG_SCENIC.h5ad")
