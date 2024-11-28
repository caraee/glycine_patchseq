library("celda")
library("Seurat")
load("../human_islets_scRNAseq/20241015_patchSeq_integrated_patched-only_decon.RData")
load("20241017_human_seurat_object_GFP.RData")

human_seurobj<-NormalizeData(human_seurobj)
sce <- as.SingleCellExperiment(human_seurobj)
sce <- decontX(sce,
               assayName = "counts")

human_seurobj <- as.Seurat(sce)
human_seurobj[["RNA"]] <- CreateAssayObject(counts = round(decontXcounts(sce)))
human_seurobj

human_islets<-subset(human_seurobj,User=="XQ"|User=="AG"|Group=="Human Primary")
human_islets<-subset(human_islets,nFeature_RNA>0)

human_islets<-NormalizeData(human_islets, verbose = FALSE)
human_islets <- FindVariableFeatures(human_islets, 
                                     selection.method = "vst", nfeatures = 500,
                                     verbose = FALSE)
pclamp_patched_filtered <- subset(x = pclamp_patched_all, 
                                  nFeature_RNA < 12500 & nFeature_RNA > 200 &
                                    percent.mt < summary(pclamp_patched_all$percent.mt)[5])


anchors <- FindTransferAnchors(
  reference = pclamp_patched_filtered,
  normalization.method = "LogNormalize",
  reduction = "pcaproject",
  reference.assay="RNA",
  query = human_islets,
  dims = 1:30,
  reference.reduction = "pca",
  verbose = T
)

human_islets <- MapQuery(anchorset = anchors, reference = pclamp_patched_filtered, 
                         query = human_islets,
                         refdata = list(celltype = "celltype"), 
                         reference.reduction = "ref.integrated.rpca", 
                         reduction.model = "umap")

DimPlot(pclamp_patched_filtered,
        raster=F,
        reduction = "umap",
        group.by = "celltype",
        order=T)

combined<-merge(pclamp_patched_filtered,human_islets)
combined<-JoinLayers(combined)
combined<-NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, 
                                 selection.method = "vst", nfeatures = 500,
                                 verbose = FALSE)
genes<-ensembldb::select(EnsDb.Hsapiens.v86, 
                         key=rownames(combined),
                         columns=c("SYMBOL"), keytype="GENEID")

circadian<-c(
  "CLOCK",
  "BMAL1",
  "PER1",
  "PER2",
  "CRY1",
  "CRY2",
  "NR1D1", #REV-ERBɑ
  "RORA"
)

circadian_ensg<-ensembldb::select(EnsDb.Hsapiens.v86, 
                                  key=circadian,
                                  columns=c("GENEID"), 
                                  keytype="SYMBOL")

circadian_ensg<-circadian_ensg$GENEID
circadian<-list(circadian_ensg)

combined<-AddModuleScore(
  combined,
  features=circadian,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "circadian_features",
  seed = 42,
  search = FALSE,
  slot = "data"
)

combined$Study[is.na(combined$Study)]<-"PEM_smartseq3"

mt_genes<-genes%>%dplyr::filter(str_detect(SYMBOL,"^MT-"))%>%
  dplyr::select(GENEID)
mt_genes<-mt_genes$GENEID
combined[["percent.mt"]] <- PercentageFeatureSet(object = combined,
                                                 features=mt_genes) 
combined <- ScaleData(combined, verbose = T,
                      vars.to.regress = c("Study","nFeature_RNA",
                                          "percent.mt","circadian_features1"))
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- FindClusters(combined, resolution = 0.8, 
                         cluster.name = "unintegrated_clusters")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", 
                    reduction.name = "umap.unintegrated")

DimPlot(combined, reduction = "umap.unintegrated", cols=cluster_cols2,
        group.by = "Study")&
  theme(axis.text = element_text(family = "Arial", color="black", size=12))&
  theme(plot.title = element_blank())

DimPlot(combined, reduction = "umap.unintegrated", cols=cluster_cols2,
        group.by = "celltype")&
  theme(axis.text = element_text(family = "Arial", color="black", size=12))&
  theme(plot.title = element_blank())

unique(human_islets$predicted.celltype)

cluster_cols2<-c("#343d00","#834efd","#d7ff44","#011b9e",
                 "#9acb00","#e953ff","#00df71","#ac00ba",
                 "#438500","#dc0086","#008e5a","#ff0072",
                 "#02edf7","#ff3401","#70a0ff","#ffdb44",
                 "#001756","#deffb0","#3b003a","#91ffe3",
                 "#790034","#a5dcff","#892f00","#a78cff",
                 "#827a00","#ff99ce","#002d00","#ff5f57",
                 "#029cb4","#ff9a62","#000f14","#ffecc1",
                 "#3a0019","#f9d2ff","#00536d")

FeaturePlot(object = combined, features = c("ENSG00000157005"), #SST
            raster=F,
            reduction = "umap.unintegrated")
sst<-WhichCells(object=combined,
                expression = `ENSG00000157005`> 7)

FeaturePlot(object = combined, features = c("ENSG00000108849"), #PPY
            raster=F,
            reduction = "umap.unintegrated")
pp<-WhichCells(object=human_islets,
               expression = `ENSG00000108849`> 7)

FeaturePlot(object = combined, features = c("ENSG00000157017"), #GHRL
            raster=F,
            reduction = "umap.unintegrated")
ghrl<-WhichCells(object=combined,
                 expression = `ENSG00000157017`> 4)

FeaturePlot(object = combined, features = c("ENSG00000091704"), #CPA1
            raster=F,
            reduction = "umap.unintegrated")
cpa1<-WhichCells(object=combined,
                 expression = `ENSG00000091704` > 5.5,
                 slot="data")

FeaturePlot(object = combined, features = c("ENSG00000171345"), #KRT19
            raster=F,
            reduction = "umap.unintegrated")
krt19<-WhichCells(object=combined,
                  expression = `ENSG00000171345` > 4.5,
                  slot="data")

FeaturePlot(object = combined, features = c("ENSG00000254647"), #INS
            raster=F,
            reduction = "umap.unintegrated")
ins<-WhichCells(object=combined,
                expression = `ENSG00000254647` > 7.5,
                slot="data")

FeaturePlot(object = combined, features = c("ENSG00000115263"), #GCG
            raster=F,
            reduction = "umap.unintegrated")
gcg<-WhichCells(object=combined,
                expression = `ENSG00000115263` > 7.5,
                slot="data")

#Check if highest expressors of canonical markers have expected cell type annotations
human_islets$celltype<-human_islets$predicted.celltype
sst<-sst[sst%in%colnames(human_islets)]
human_islets$predicted.celltype[sst]
human_islets$predicted.celltype.score[sst]
sst<-sst[human_islets$predicted.celltype.score[sst]<0.7]
human_islets$celltype[sst]<-"delta"

cpa1<-cpa1[cpa1%in%colnames(human_islets)]
human_islets$predicted.celltype[cpa1]
human_islets$predicted.celltype.score[cpa1]
cpa1<-cpa1[human_islets$predicted.celltype.score[cpa1]<0.7]
human_islets$celltype[cpa1]<-"acinar"

krt19<-krt19[krt19%in%colnames(human_islets)]
human_islets$predicted.celltype[krt19]
human_islets$predicted.celltype.score[krt19]
krt19<-krt19[human_islets$predicted.celltype.score[krt19]<0.7]
DimPlot(combined,
        raster=F,
        reduction = "umap.unintegrated",
        cells.highlight=krt19,
        cols.highlight = "#053059",
        order=T)
human_islets$celltype[krt19]<-"ductal"

pp<-pp[pp%in%colnames(human_islets)]
human_islets$predicted.celltype[pp]
human_islets$predicted.celltype.score[pp]
pp<-pp[human_islets$predicted.celltype.score[pp]<0.7]
DimPlot(combined,
        raster=F,
        reduction = "umap.unintegrated",
        cells.highlight=pp[4],
        cols.highlight = "#053059",
        order=T)
human_islets$celltype[pp[4]]<-"PP"

ghrl<-ghrl[ghrl%in%colnames(human_islets)]
human_islets$predicted.celltype[ghrl]
human_islets$predicted.celltype.score[ghrl]
ghrl<-ghrl[human_islets$predicted.celltype.score[ghrl]<0.7]
human_islets$celltype[ghrl]<-"epsilon"

ins<-ins[ins%in%colnames(human_islets)]
human_islets$predicted.celltype[ins]
human_islets$predicted.celltype.score[ins]
ins<-ins[human_islets$predicted.celltype.score[ins]<0.7]
DimPlot(combined,
        raster=F,
        reduction = "umap.unintegrated",
        cells.highlight=ins,
        cols.highlight = "#053059",
        order=T)
human_islets$celltype[ins]<-"beta"

gcg<-gcg[gcg%in%colnames(human_islets)]
human_islets$predicted.celltype[gcg]
human_islets$predicted.celltype.score[gcg]
gcg<-gcg[human_islets$predicted.celltype.score[gcg]<0.7]
DimPlot(combined,
        raster=F,
        reduction = "umap.unintegrated",
        cells.highlight=gcg,
        cols.highlight = "#053059",
        order=T)
human_islets$celltype[gcg]<-"alpha"

intersect(ins,gcg)

DimPlot(combined, reduction = "umap.unintegrated", cols=cluster_cols2,
        group.by = "celltype", label = F, label.size = 3,
        raster=F,
        repel = TRUE)

AG_seurobj_ensg<-subset(human_seurobj,User=="AG")
AG_meta<-human_meta%>%dplyr::filter(User=="AG")
AG_meta_clean <- AG_meta %>% 
  dplyr::select(where(~ any(!is.na(.))))
AG_meta_clean <- AG_meta_clean%>%dplyr::filter(CellID!="BLANK")
AG_meta_clean <- AG_meta_clean%>%dplyr::filter(NAME_FROM_MBSU%in%colnames(AG_seurobj_ensg))
AG_meta_clean$CellID<-str_replace_all(AG_meta_clean$CellID,"00","0") 
AG_meta_clean$Donor<-case_when(AG_meta_clean$Donor=="504"~"R504",
                               AG_meta_clean$Donor=="506"~"R506",
                               AG_meta_clean$Donor=="H252"~"H2522",
                               AG_meta_clean$Donor=="HPAP-150"~"HPAP150",
                               TRUE~AG_meta_clean$Donor)
AG_meta_clean$Donor<-str_replace_all(AG_meta_clean$Donor,"HPAP","HPAP-")

AG_meta_full<-full_join(AG_meta_clean,dat,keep = FALSE)
setdiff(AG_meta_full$NAME_FROM_MBSU,colnames(AG_seurobj_ensg))

AG_meta_final<-AG_meta_full%>%dplyr::select(
  -c(Species,Well,Group,Note,Comments,Sample,`Date (6 CHARACTER DATE)`,
     celltype_old, celltype,IGly,User,Date)
)
head(AG_meta_final)
AG_meta_final<-AG_meta_final%>%rename("Glycinecurrent+strychnine(pA/pF)"=I.Gly.Stryc,
                                      "Glycinecurrentwashout(pA/pF)"=IGly.w.)

write_csv(AG_meta_final,file="20241017_AG_metadata.csv")

AG_mat<-JM_mat[,AG_meta_final$NAME_FROM_MBSU]
write_csv(AG_mat,"20241017_AG_full_raw_counts.csv")

counts<-AG_seurobj_ensg[["RNA"]]$counts

rownames(AG_meta_full)<-AG_meta_full$NAME_FROM_MBSU

gene_ids<-rownames(counts)
tr2g<-tibble(gene=row.names(counts),
             gene_name=mapIds(EnsDb.Hsapiens.v86,
                              keys=row.names(counts),
                              column="SYMBOL",
                              keytype="GENEID",
                              multiVals="first"))

keeps<-tr2g[!is.na(tr2g$gene_name),]
counts<-counts[keeps$gene,]

rownames(counts)<-keeps$gene_name

counts2<- rowsum(counts, rownames(counts))%>%as("dgCMatrix")
counts2[1:5,1:5]

AG_seurobj <- Seurat::CreateSeuratObject(counts = counts2,
                                         assay = "RNA", 
                                         meta.data = AG_meta_full,
                                         min.cells=1,
                                         min.features=200)

AG_seurobj[["percent.mt"]] <- PercentageFeatureSet(object = AG_seurobj,
                                                   pattern = "^MT-")

AG_seurobj<-subset(AG_seurobj,percent.mt<=30)

AG_seurobj <- NormalizeData(object = AG_seurobj,
                            normalization.method = "RC",
                            scale.factor = 1e6)

AG_dat<-tibble(AG_ID=AG_seurobj$CellID,
               Genes=AG_seurobj$nFeature_RNA,
               Counts=AG_seurobj$nCount_RNA,
               percent.mt=round(AG_seurobj$percent.mt,2),
               celltype=AG_seurobj$celltype)

AG_dat %>% 
  DT::datatable()

AG_seurobj$Diabetes.Status[AG_seurobj$Diabetes.Status=="T2DM"]<-"T2D"

table(AG_seurobj$celltype)
AG_seurobj_beta<-subset(AG_seurobj,celltype=="beta")

counts<-AG_seurobj[["RNA"]]$counts
genes.percent.expression <- rowMeans(counts>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>0]) 

metadata<-AG_seurobj@meta.data

counts<-counts[genes.filter,]
AG_sce<-SingleCellExperiment(assays=list(counts=counts),
                             colData=metadata)
sceasy::convertFormat(AG_sce, from="sce", to="anndata",
                      main_layer="counts", 
                      drop_single_values=FALSE,
                      outFile='20241020_AG_symbs.h5ad')

###### create counts table for Spearman Rank Correlations ######
counts<-AG_seurobj_beta[["RNA"]]$counts
genes.percent.expression <- rowMeans(counts>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>50]) 

data<-AG_seurobj_beta[["RNA"]]$data
data<-data[genes.filter,]

metadata<-AG_seurobj_beta@meta.data

counts<-counts[genes.filter,]
AG_beta_sce<-SingleCellExperiment(assays=list(counts=counts),
                                  colData=metadata)
sceasy::convertFormat(AG_beta_sce, from="sce", to="anndata",
                      main_layer="counts", 
                      drop_single_values=FALSE,
                      outFile='20241020_AG_beta_symbs.h5ad')
data<-data.frame(data)

data<-data%>%mutate(Index=rownames(data),.keep="all",.before=everything())
data[1:5,1:5]

metadata<-t(metadata)%>%data.frame()
metadata[1:5,1:5]
metadata<-metadata%>%mutate(Index=rownames(metadata),.keep="all",
                            .before=everything())

bad_dat<-rbind(metadata,data)

write_csv(bad_dat,file="20241020_AG_cpm_beta_cells.csv")