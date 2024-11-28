##### AG mixOmics #####
library("mixOmics")
library("scico")
library("tidyverse")
library("Seurat")
library("reticulate")

load("20241020_AG_seurobj_symbs.RData")

sc <- import('scanpy', convert = FALSE)
adata = sc$read_h5ad("20241018_adata_AG_SCENIC.h5ad")

auc_mtx<-adata$obsm[["RegulonsAUC"]]%>%py_to_r%>%as.data.frame%>%t
AUC_assay <- CreateAssayObject(counts = auc_mtx)
AG_seurobj[["AUC"]] <- AUC_assay

AG_seurobj$circadian_score1<-human_islets$circadian_features1[colnames(AG_seurobj)]

AG_seurobj[["percent.ribo"]] <- PercentageFeatureSet(object = AG_seurobj, 
                                                     pattern = "^RP[LS]")
summary(AG_seurobj$circadian_score1)
summary(AG_seurobj$percent.ribo)

AG_beta<-subset(AG_seurobj,celltype=="beta")
AG_beta$Diabetes.Status[AG_beta$Diabetes.Status=="T2DM"]<-"T2D"

auc_mtx<-AG_beta[["AUC"]]$counts%>%t%>%as.matrix()
nzv_scenic<-nearZeroVar(auc_mtx, 
                        freqCut = 95/5, uniqueCut = 10)
auc_mtx<-auc_mtx[,-nzv_scenic$Position]

auc_mtx_scaled<-scale(auc_mtx)

AG_beta<-NormalizeData(AG_beta)
cor_trans<-read_csv("20241023_Updated_Amanda_Correlations.csv") #spearman rank correlations

features<-cor_trans%>%dplyr::filter(abs(Correlation)>0.3)%>%dplyr::select(Variable_2_Transcript)%>%unique

AG_beta<-ScaleData(AG_beta, verbose = T,
                   vars.to.regress = c("nFeature_RNA",
                                       "percent.mt",
                                       "percent.ribo",
                                       "circadian_score1"),
                   features = features$Variable_2_Transcript)

mrna_mtx<-t(AG_beta[["RNA"]]$scale.data)
mrna_mtx[1:5,1:5]
dim(mrna_mtx)

sum("INS"%in%colnames(mrna_mtx))

#remove ribosomal transcripts and pseudogenes to make interpretation easier with small sample size
sum(str_detect(colnames(mrna_mtx),"^RP[LS]"))
mrna_mtx<-mrna_mtx[,!str_detect(colnames(mrna_mtx),"^RP[LS]")]
mrna_mtx<-mrna_mtx[,!str_detect(colnames(mrna_mtx),"^RP11")]
mrna_mtx<-mrna_mtx[,!str_detect(colnames(mrna_mtx),"^RP3")]
mrna_mtx<-mrna_mtx[,!str_detect(colnames(mrna_mtx),"^RP5")]
sum(str_detect(colnames(mrna_mtx),"^MT"))
mrna_mtx<-mrna_mtx[,!str_detect(colnames(mrna_mtx),"^MT")]
mrna_mtx[1:5,1:5]
mrna_mtx<-mrna_mtx[,-(1)]
dim(mrna_mtx)

ephys<-AG_beta@meta.data[c(29:32)]
ephys_scaled <- scale(ephys)  # cell size normalized data

X <- list(mRNA = mrna_mtx, 
          ephys = ephys_scaled, 
          scenic = auc_mtx_scaled)

Y <-AG_beta$Diabetes.Status
summary(Y)

res1.pls.glyR <- pls(X$mRNA, X$ephys, ncomp = 1)
cor(res1.pls.glyR$variates$X, res1.pls.glyR$variates$Y)

res2.pls.glyR <- pls(X$mRNA, X$scenic, ncomp = 1)
cor(res2.pls.glyR$variates$X, res2.pls.glyR$variates$Y)

res3.pls.glyR <- pls(X$scenic, X$ephys, ncomp = 1)
cor(res3.pls.glyR$variates$X, res3.pls.glyR$variates$Y)

design <- matrix(1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design

diablo.glyR <- block.plsda(X, Y, ncomp = 4, design = design, scale = F)

perf.diablo.glyR = perf(diablo.glyR, validation = 'Mfold', folds = 10, nrepeat = 25)

perf.diablo.glyR$error.rate  # Lists the different types of error rates

# Plot of the error rates based on weighted vote
plot(perf.diablo.glyR)

perf.diablo.glyR$choice.ncomp$WeightedVote
ncomp <- perf.diablo.glyR$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]

set.seed(42) # for reproducibility
test.keepX <- list(mRNA = c(25, 50, 100, 200),
                   ephys = c(1, 2, 3),
                   scenic = c(25, 50, 100, 150))

tune.diablo.glyR <- tune.block.splsda(X, Y, ncomp = ncomp,
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = 10, nrepeat = 100, 
                                      BPPARAM = BiocParallel::SnowParam(workers = 8),
                                      dist = "centroids.dist",
                                      scale = F)

(list.keepX <- tune.diablo.glyR$choice.keepX)

diablo.glyR <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design,
                            scale = F)
save(diablo.glyR,file="20241029_ASG_mixOmics_filtered_norm.RData")


pop_colours <- c("black",scico::scico(2,palette = "roma")[2])
names(pop_colours)<-names(pop_colours)<-c("ND","T2D")
block_colours<-scico::scico(3,palette = "roma",end=0.8) 
ccp_cols<-scico::scico(20,palette = "roma")[c(3,20)]

corMat.P <- circosPlot(diablo.glyR, cutoff = 1, line = TRUE, #comp=3,
                       size.variables = 0.5,
                       blocks.link="ephys",
                       color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
                       color.cor = c("chocolate3","grey20"), size.labels = 1.25)
dim(corMat.P)
corMat.P[corMat.P < -1] <- -1
corMat.P[corMat.P > 1] <- 1

features_to_keep <- apply(abs(corMat.P) > 0.5, 1, any)
scenic_to_keep<-rownames(corMat.P)%>%str_detect("\\(\\+\\)")
feats<-rownames(corMat.P)[features_to_keep&scenic_to_keep]

pdf(file = "./Figures/20241106_ASG_circos_norm_no_mRNA_roma_roma_text.pdf",
    width=7.75,height=7.75)
circosPlot(diablo.glyR, cutoff = 0.7, line = F, #comp=1:1,
           size.variables = 0.475,
           var.adj=-0.7,
           blocks.link="ephys",
           blocks=c(2:3),
           color.Y = pop_colours,
           color.blocks = block_colours[2:3],
           color.cor = ccp_cols,
           size.labels = 1)

dev.off()

cor_df<-corMat.P%>%data.frame%>%tibble%>%
  dplyr::mutate(RowIndex=colnames(corMat.P),
                .keep="all",
                .before=everything())

cor_df

write_csv(cor_df,file="20241106_ASG_cormat.csv")

ephys_params <- diablo.glyR$loadings$ephys

rows_to_keep <- apply(abs(corMat.P) > 0.6, 1, any)
cols_to_keep <- apply(abs(corMat.P) > 0.6, 2, any)

ephys_rows <- rownames(corMat.P) %in% ephys_params
ephys_cols <- colnames(corMat.P) %in% ephys_params

rows_to_keep <- rows_to_keep | ephys_rows
cols_to_keep <- cols_to_keep | ephys_cols

corMat_filtered <- corMat.P[rows_to_keep, cols_to_keep]
dim(corMat_filtered)

cor_df_filt<-corMat_filtered%>%data.frame%>%tibble%>%
  dplyr::mutate(RowIndex=colnames(corMat_filtered),
                .keep="all",
                .before=everything())

write_csv(cor_df_filt%>%tibble(),file="20241128_ASG_cormat_filtered.csv")

hclust<-corrplot::corrplot(corMat_filtered, order="hclust",tl.col="black",
                           col=scico::scico(n=10, palette = "roma"),
                           addrect = 4,tl.cex = 0.55,
                           diag=F)
order<-rownames(hclust$corr)
order <- order[order != "I.Gly.Stryc"]
position <- which(order == "IGly.w.")
order <- append(order, "I.Gly.Stryc", after = position)
corMat_filtered<-corMat_filtered[order,order]

pdf(file = "./Figures/20241128_ASG_corrplot_norm_roma_modified.pdf",
    width=7,height=7,pointsize = 16)
corrplot::corrplot(corMat_filtered, order="original",tl.col="black",
                   col=scico::scico(n=10, palette = "roma"),
                   addrect = 4,tl.cex = 0.55,
                   diag=F)
dev.off()