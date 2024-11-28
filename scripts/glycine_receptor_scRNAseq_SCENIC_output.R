library("reticulate")
library("sceasy")
library("SCopeLoomR")
library("AUCell")
library("SCENIC")

sc <- import('scanpy', convert = FALSE)
adata = sc$read_h5ad("20241018_adata_AG_SCENIC.h5ad")
load("20241020_AG_seurobj_symbs.RData")

scenicLoomPath <- "./20241018_AG_pyscenic_final.loom"
motifEnrichmentFile <- "./AG_reg.csv"


loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellClusters <- get_clusterings(loom)
close_loom(loom)

auc_mtx<-adata$obsm[["RegulonsAUC"]]%>%py_to_r%>%as.data.frame%>%t
AUC_assay <- CreateAssayObject(counts = auc_mtx)
AG_seurobj[["AUC"]] <- AUC_assay

auc_umap<-adata$obsm[["X_SCENIC_AUC_UMAP"]]%>%py_to_r
rownames(auc_umap)<-colnames(AG_seurobj)

AG_seurobj[["SCENIC"]] <- CreateDimReducObject(embeddings = auc_umap, 
                                               key = "scenic_", assay = "AUC")

motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

AG_seurobj$Diabetes.Status["C12-Sept-23-Plate2-PEM_S132"]<-"T2D"

AG_beta<-subset(AG_seurobj,celltype=="beta")

DiabetesStatus<-AG_beta$Diabetes.Status

# Split the cells by cluster:
cellsPerCluster <- split(names(AG_seurobj$Diabetes.Status),as.character(DiabetesStatus))
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[complete.cases(regulonActivity_byCellType_Scaled), ]
# plot:
#options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- ComplexHeatmap::draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, 
                                                   name="Regulon activity",
                                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[ComplexHeatmap::row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "DiabetesStatus", "RelativeActivity")
topRegulators$DiabetesStatus <- factor(as.character(topRegulators$DiabetesStatus))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

topRegulators%>%dplyr::filter(DiabetesStatus=="ND")%>%slice_max(n=15, order_by = RelativeActivity)
topRegulators%>%dplyr::filter(DiabetesStatus=="T2D")%>%slice_max(n=15, order_by = RelativeActivity)

viewTable(topRegulators, options = list(pageLength = 10))
topRegulators<-topRegulators%>%dplyr::group_by(DiabetesStatus)%>%
  dplyr::arrange(desc(RelativeActivity),.by_group=T)

write_csv(topRegulators,file = "AG_SCENIC_top_regulons_beta_diabetes.csv")

cellsSelected<-colnames(AG_beta)

regulonAUC_subset <- regulonAUC[regulonOrder, which(colnames(regulonAUC) %in% cellsSelected)]

rss <- calcRSS(AUC=getAUC(regulonAUC_subset), 
               cellAnnotation=DiabetesStatus[colnames(regulonAUC_subset)])
rss <- rss[complete.cases(rss), ]

rss_df<-rss%>%data.frame%>%tibble%>%
  dplyr::mutate(regulon=rownames(rss),
                .keep="all",
                .before="T2D")
rssNorm_df<-rss%>%scale%>%data.frame%>%tibble%>%
  dplyr::mutate(regulon=rownames(rss),
                .keep="all",
                .before="T2D")
rss_tibble <- full_join(rss_df,
                        by="regulon",
                        rssNorm_df,
                        suffix = c(".RSS",".Z"))

write_csv(rss_tibble,file="20241106_rss_AG.csv")

rss_cols<-scico(3,palette = "roma")

rssPlot_filt <- plotRSS(rss, thr=0.2,
                        zThreshold = 1.105,
                        #revCol = T,
                        col.low = rss_cols[1],
                        col.mid = rss_cols[2],
                        col.high = rss_cols[3],)
p_filt<-plotly::ggplotly(rssPlot_filt$plot)
df_filt<-rssPlot_filt$df

htmlwidgets::saveWidget(p_filt, "20241106_ASG_RSS_filt.html")
orca(p_filt, "./Figures/20241106_ASG_RSS_filt_39.pdf",
     width = 400,
     height = 500)

png(filename = "./Figures/20241018_SCENIC_rss_ND_AG.png",width=12,height=12,units="cm",res=600)
plotRSS_oneSet(rss, setName = "ND")
dev.off()

png(filename = "./Figures/20241018_SCENIC_rss_T2D_AG.png",width=12,height=12,units="cm",res=600)
plotRSS_oneSet(rss, setName = "T2D")
dev.off()

regulons_df <- data.frame(
  Regulon = rep(names(regulons), sapply(regulons, length)), 
  Genes = unlist(regulons)
)

# Write the data frame to a CSV file
write.csv(regulons_df, file = "AG_regulons_output.csv", row.names = FALSE)
