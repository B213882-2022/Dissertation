library(SingleCellExperiment)
library(scater)
library(scran)
library(biomaRt)
library(pheatmap)
library(data.table)

# read scRNA seq data (raw count) and cell annotation
raw_count <- fread('counts.txt', skip=1, header=TRUE, data.table = FALSE)
rownames(raw_count) <- raw_count$Geneid
raw_count <- raw_count[,7:ncol(raw_count)]
raw_count <- as.matrix(raw_count)
anno <- read.csv('annotation.csv')
idx <- match(colnames(raw_count),anno$Dir)
colnames(raw_count) <- anno$Run[idx]
anno <- anno[idx,]
coldata <- anno[,c('Cell_type','sex','Strain')]
rownames(coldata) <- anno$Run
colnames(coldata) <- c('Cell_Type','Sex','Strain')
dim(raw_count)
sce <- SingleCellExperiment(assay=list(counts = raw_count),colData=coldata)
sce

# gene annotation
ah <- AnnotationHub()
#display(ah)  # search 'GRCm39' and 'EnsDb' to find v109 ID
ens.mm.v109 <- AnnotationHub()[["AH109655"]]
columns(ens.mm.v109)
gene_anno <- AnnotationDbi::select(ens.mm.v109, keys=rownames(sce), 
                                   keytype='GENEID', column='GENENAME')
rowData(sce)$GENENAME <- make.names(gene_anno$GENENAME[match(rownames(sce),gene_anno$GENEID)],
                                  unique = TRUE)
rowData(sce)$ENSEMBL <- rownames(sce)
rownames(sce) <- rowData(sce)$GENENAME
count(grepl("^ERCC", rowData(sce)$ENSEMBL))
is.spike <- grepl("^ERCC", rowData(sce)$ENSEMBL)
spike_names <- rowData(sce)$ENSEMBL[grepl("^ERCC", rowData(sce)$ENSEMBL)]
rownames(sce)[is.spike] <- spike_names
sce
head(rowData(sce))
length(grep('^NA',rowData(sce)$GENENAME))  # the number of failures in gene name annotation
length(grep('^NA',rownames(sce)))

# Get Mitochondrial genes' ENSEMBL ID
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
head(listDatasets(ensembl))
searchDatasets(mart = ensembl, pattern = "GRCm39")
ensembl <- useDataset(dataset='mmusculus_gene_ensembl', mart=ensembl)
head(listFilters(ensembl))
head(listAttributes(ensembl))
MT_genes <- getBM(filters = 'chromosome_name', values='MT',
                  attributes = c('ensembl_gene_id','external_gene_name'), 
                  mart = ensembl)
head(MT_genes)
MT_genes_id <- paste(MT_genes$ensembl_gene_id, collapse='|')

# cell QC
count(grepl(MT_genes_id, rowData(sce)$ENSEMBL))
is.mito <- grepl(MT_genes_id, rowData(sce)$ENSEMBL)
sce$stats <- perCellQCMetrics(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
# or we can store the stats directly to colData 
# by addPerCellQC(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
head(sce$stats)
#Plot the results
png('QC_raw.png')
par(mfrow=c(2,2))
hist(sce$stats$sum/1e6, xlab="Library sizes (millions)", main="", 
     breaks=50, col="grey80", ylab="Number of cells")
hist(sce$stats$detected, xlab="Number of detected genes", main="", 
     breaks=50, col="grey80", ylab="Number of cells")
hist(sce$stats$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=50, main="", col="grey80")
hist(sce$stats$subsets_ERCC_percent, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=50, main="", col="grey80")
dev.off()
# remove QC failed cells
reasons <- perCellQCFilters(sce$stats,sub.fields=c("subsets_Mt_percent", "subsets_ERCC_percent"))
colSums(as.matrix(reasons))  # check how many cells were dropped 
summary(reasons$discard)
sce_origin <- sce
#sce <- sce[, !reasons$discard]  # drop 40 cells in total
dim(sce)
# perCellQCFilters is equal to the following command
#libsize.drop <- isOutlier(sce$stats$sum, nmads=3, type='lower', log=TRUE)
#sum(libsize.drop)
#feature.drop <- isOutlier(sce$stats$detected, nmads=3, type='lower', log=TRUE)
#sum(feature.drop)
#mito.drop <- isOutlier(sce$stats$subsets_Mt_percent, nmads=3, type='higher')
#sum(mito.drop)
#spike.drop <- isOutlier(sce$stats$subsets_ERCC_percent, nmads=3, type='higher')
#sum(spike.drop)
#discard.all <- ncol(sce[,(libsize.drop | feature.drop | mito.drop | spike.drop)])
#discard.all

# check QC after filtering
sce$stats_selected_cells <- perCellQCMetrics(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
png('QC_filtered.png')
par(mfrow=c(2,2))
hist(sce$stats_selected_cells$sum/1e6, xlab="Library sizes (millions)", main="", 
     breaks=50, col="grey80", ylab="Number of cells")
hist(sce$stats_selected_cells$detected, xlab="Number of detected genes", main="", 
     breaks=50, col="grey80", ylab="Number of cells")
hist(sce$stats_selected_cells$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=50, main="", col="grey80")
hist(sce$stats_selected_cells$subsets_ERCC_percent, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=50, main="", col="grey80")
dev.off()
# QC summary
sce_origin$discard <- reasons$discard
sce_origin$sum <- sce_origin$stats$sum/1e6
sce_origin$detected <- sce_origin$stats$detected
sce_origin$mito <- sce_origin$stats$subsets_Mt_percent
sce_origin$ERCC <- sce_origin$stats$subsets_ERCC_percent
QC_sum <- gridExtra::grid.arrange(
  plotColData(sce_origin, x="Sex", y="sum", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
              ggtitle("Library sizes") + ylab('Total Reads Count (millions)') +
              xlab('Sex'),
  plotColData(sce_origin, x="Sex", y="detected", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
              ggtitle("Detected genes") + ylab('Number of detected genes') +
              xlab('Sex'),
  plotColData(sce_origin, x="Sex", y="mito", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
              ggtitle("Mito percent") + ylab('Mitochondrial proportion (%)') +
              xlab('Sex'),
  plotColData(sce_origin, x="Sex", y="ERCC", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
              ggtitle("ERCC percent") + ylab('ERCC proportion (%)') +
              xlab('Sex'),
  ncol=1
)
ggsave('QC_summary.jpg',QC_sum, device='jpg', width = 10, height = 30)

QC_sum2 <- gridExtra::grid.arrange(
  plotColData(sce_origin, x="sum", y="detected", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
              ggtitle("Total counts VS Detected genes per cell") +
              ylab('Number of detected genes') + xlab('Total Reads Count (millions)'),
  plotColData(sce_origin, x="sum", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
              ggtitle("Total counts VS Mitochondrial proportion") +
              ylab('Mitochondrial proportion (%)') + xlab('Total Reads Count (millions)'),
  plotColData(sce_origin, x="sum", y="ERCC", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
              ggtitle("Total counts VS ERCC proportion") +
              ylab('ERCC proportion (%)') + xlab('Total Reads Count (millions)'),
  plotColData(sce_origin, x="detected", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
              ggtitle("Detected genes per cell VS Mitochondrial proportion") +
              ylab('Mitochondrial proportion (%)') + xlab('Number of detected genes'),
  plotColData(sce_origin, x="detected", y="ERCC", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
              ggtitle("Detected genes per cell VS ERCC proportion") +
              ylab('ERCC proportion (%)') + xlab('Number of detected genes'),
  plotColData(sce_origin, x="ERCC", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
              ggtitle("ERCC proportion VS Mitochondrial proportion") +
              ylab('Mitochondrial proportion (%)') + xlab('ERCC proportion (%)'),
  ncol=1
)
ggsave('QC_summary2.jpg',QC_sum2, device='jpg', width = 10, height = 40, limitsize = FALSE)

# show top 100 highly expressed genes
plotHighestExprs(sce_origin, n=100, exprs_values = "counts",
                 feature_names_to_plot = "GENENAME", colour_cells_by="detected") +
                  xlab('Percentage of counts') + ylab('Gene names')
ggsave('Highly_expressed_genes.jpg', device='jpg',width=8, height = 20)

# seperate E3.5 and E4.5
table(sce$Cell_Type)
sce_34 <- sce[,sce$Cell_Type %in% c('E3.5 ICM','E4.5 Epiblast')]
sce_34$stats <- perCellQCMetrics(sce_34, subsets=list(ERCC=is.spike, Mt=is.mito))
sce_34_reasons <- perCellQCFilters(sce_34$stats,sub.fields=c("subsets_Mt_percent", "subsets_ERCC_percent"))
sce_34$discard <- sce_34_reasons$discard
sce_34$sum <- sce_34$stats$sum/1e6
sce_34$detected <- sce_34$stats$detected
sce_34$mito <- sce_34$stats$subsets_Mt_percent
sce_34$ERCC <- sce_34$stats$subsets_ERCC_percent
QC_34_sum <- gridExtra::grid.arrange(
  plotColData(sce_34, x="Sex", y="sum", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("Library sizes") + ylab('Total Reads Count (millions)') +
    xlab('Sex'),
  plotColData(sce_34, x="Sex", y="detected", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("Detected genes") + ylab('Number of detected genes') +
    xlab('Sex'),
  plotColData(sce_34, x="Sex", y="mito", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("Mito percent") + ylab('Mitochondrial proportion (%)') +
    xlab('Sex'),
  plotColData(sce_34, x="Sex", y="ERCC", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("ERCC percent") + ylab('ERCC proportion (%)') +
    xlab('Sex'),
  ncol=1
)
ggsave('QC_34_summary.jpg',QC_34_sum, device='jpg', width = 10, height = 20)

QC_34_sum2 <- gridExtra::grid.arrange(
  plotColData(sce_34, x="sum", y="detected", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Total counts VS Detected genes per cell") +
    ylab('Number of detected genes') + xlab('Total Reads Count (millions)'),
  plotColData(sce_34, x="sum", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Total counts VS Mitochondrial proportion") +
    ylab('Mitochondrial proportion (%)') + xlab('Total Reads Count (millions)'),
  plotColData(sce_34, x="sum", y="ERCC", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Total counts VS ERCC proportion") +
    ylab('ERCC proportion (%)') + xlab('Total Reads Count (millions)'),
  plotColData(sce_34, x="detected", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Detected genes per cell VS Mitochondrial proportion") +
    ylab('Mitochondrial proportion (%)') + xlab('Number of detected genes'),
  plotColData(sce_34, x="detected", y="ERCC", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Detected genes per cell VS ERCC proportion") +
    ylab('ERCC proportion (%)') + xlab('Number of detected genes'),
  plotColData(sce_34, x="ERCC", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("ERCC proportion VS Mitochondrial proportion") +
    ylab('Mitochondrial proportion (%)') + xlab('ERCC proportion (%)'),
  ncol=1
)
ggsave('QC_34_summary2.jpg',QC_34_sum2, device='jpg', width = 10, height = 40, limitsize = FALSE)

# E3.5 ICM and E4.5
# drop QC failed cells
sce_34_qc <- sce_34[, !sce_34_reasons$discard]
sce_34_qc
# move spikes to AltExp
sce_34_qc <- splitAltExps(sce_34_qc, is.spike)
altExpNames(sce_34_qc) <- 'spikes'
sce_34_qc
# Normalisation
E3.5_num <- sum(sce_34_qc$Cell_Type %in% 'E3.5 ICM')
E4.5_num <- sum(sce_34_qc$Cell_Type %in% 'E4.5 Epiblast')
sce_34_qc <- computeSumFactors(sce_34_qc, size=c(E3.5_num, E4.5_num))
#clust_sce_34_qc <- quickCluster(sce_34_qc, min.size=10)
#table(clust_sce_34_qc)
#sce_34_qc <- computeSumFactors(sce_34_qc, clusters = clust_sce_34_qc)
sizeFactors(sce_34_qc)
plot(sizeFactors(sce_34_qc),sce_34_qc$stats$sum/1e6, log='xy')
sce_34_qc <- logNormCounts(sce_34_qc)
assayNames(sce_34_qc)
#set.seed(213882)
#set.seed(NULL)
# by default use top 500 HVG, but it would be better to calculate HVG ourselves and use 'subset_row'
sce_34_qc <- runPCA(sce_34_qc, ntop=500)  
#plot(attr(reducedDim(sce_34_qc, 'PCA'), "percentVar"), 
#     log="y", xlab="PC", ylab="Variance explained (%)")
sce_34_qc <- runTSNE(sce_34_qc,perplexity=10)
sce_34_qc <- runUMAP(sce_34_qc)
reducedDimNames(sce_34_qc)
cluster_34_plot <- gridExtra::grid.arrange(
  plotReducedDim(sce_34_qc, "PCA", colour_by = 'Cell_Type')
                +ggtitle('PCA'),
  plotReducedDim(sce_34_qc, "TSNE", colour_by = 'Cell_Type')
                +ggtitle('TSNE'),
  plotReducedDim(sce_34_qc, "UMAP", colour_by = 'Cell_Type')
                +ggtitle('UMAP'),
ncol=1)
ggsave('cluster_34.jpg',cluster_34_plot, device='jpg', width = 6, height = 12)

# Rest cell
sce_rest <- sce[,!sce$Cell_Type %in% c('E3.5 ICM','E4.5 Epiblast')]
sce_rest
table(sce_rest$Cell_Type)
sce_rest$stats <- perCellQCMetrics(sce_rest, subsets=list(ERCC=is.spike, Mt=is.mito))
sce_rest_reasons <- perCellQCFilters(sce_rest$stats,sub.fields=c("subsets_Mt_percent", "subsets_ERCC_percent"))
sce_rest$discard <- sce_rest_reasons$discard
sce_rest$sum <- sce_rest$stats$sum/1e6
sce_rest$detected <- sce_rest$stats$detected
sce_rest$mito <- sce_rest$stats$subsets_Mt_percent
sce_rest$ERCC <- sce_rest$stats$subsets_ERCC_percent
QC_rest_sum <- gridExtra::grid.arrange(
  plotColData(sce_rest, x="Sex", y="sum", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("Library sizes") + ylab('Total Reads Count (millions)') +
    xlab('Sex'),
  plotColData(sce_rest, x="Sex", y="detected", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("Detected genes") + ylab('Number of detected genes') +
    xlab('Sex'),
  plotColData(sce_rest, x="Sex", y="mito", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("Mito percent") + ylab('Mitochondrial proportion (%)') +
    xlab('Sex'),
  plotColData(sce_rest, x="Sex", y="ERCC", colour_by="discard",
              other_fields="Cell_Type") + facet_wrap(~Cell_Type) +
    ggtitle("ERCC percent") + ylab('ERCC proportion (%)') +
    xlab('Sex'),
  ncol=1
)
ggsave('QC_rest_summary.jpg',QC_rest_sum, device='jpg', width = 10, height = 20)

QC_rest_sum2 <- gridExtra::grid.arrange(
  plotColData(sce_rest, x="sum", y="detected", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Total counts VS Detected genes per cell") +
    ylab('Number of detected genes') + xlab('Total Reads Count (millions)'),
  plotColData(sce_rest, x="sum", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Total counts VS Mitochondrial proportion") +
    ylab('Mitochondrial proportion (%)') + xlab('Total Reads Count (millions)'),
  plotColData(sce_rest, x="sum", y="ERCC", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Total counts VS ERCC proportion") +
    ylab('ERCC proportion (%)') + xlab('Total Reads Count (millions)'),
  plotColData(sce_rest, x="detected", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Detected genes per cell VS Mitochondrial proportion") +
    ylab('Mitochondrial proportion (%)') + xlab('Number of detected genes'),
  plotColData(sce_rest, x="detected", y="ERCC", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("Detected genes per cell VS ERCC proportion") +
    ylab('ERCC proportion (%)') + xlab('Number of detected genes'),
  plotColData(sce_rest, x="ERCC", y="mito", colour_by="discard", 
              other_fields = "Cell_Type") + facet_wrap(~Cell_Type) + 
    ggtitle("ERCC proportion VS Mitochondrial proportion") +
    ylab('Mitochondrial proportion (%)') + xlab('ERCC proportion (%)'),
  ncol=1
)
ggsave('QC_rest_summary2.jpg',QC_rest_sum2, device='jpg', width = 10, height = 40, limitsize = FALSE)
# cluster
sce_rest_qc <- sce_rest[, !sce_rest_reasons$discard]
sce_rest_qc <- splitAltExps(sce_rest_qc, is.spike)
altExpNames(sce_rest_qc) <- 'spikes'
clust_sce_rest_qc <- quickCluster(sce_rest_qc, min.size=10)
table(clust_sce_rest_qc)
sce_rest_qc <- computeSumFactors(sce_rest_qc, clusters = clust_sce_rest_qc)
sizeFactors(sce_rest_qc)
plot(sizeFactors(sce_rest_qc),sce_rest_qc$stats$sum/1e6, log='xy')
sce_rest_qc <- logNormCounts(sce_rest_qc)
assayNames(sce_rest_qc)
sce_rest_qc <- runPCA(sce_rest_qc)
sce_rest_qc <- runTSNE(sce_rest_qc,perplexity=10)
sce_rest_qc <- runUMAP(sce_rest_qc)
reducedDimNames(sce_rest_qc)
cluster_rest_plot <- gridExtra::grid.arrange(
  plotReducedDim(sce_rest_qc, "PCA", colour_by = 'Cell_Type')
  +ggtitle('PCA'),
  plotReducedDim(sce_rest_qc, "TSNE", colour_by = 'Cell_Type')
  +ggtitle('TSNE'),
  plotReducedDim(sce_rest_qc, "UMAP", colour_by = 'Cell_Type')
  +ggtitle('UMAP'),
  ncol=1)
ggsave('cluster_rest.jpg',cluster_rest_plot, device='jpg', width = 6, height = 12)

# all cells
length(c(colnames(sce_rest_qc),colnames(sce_34_qc)))
sce <- sce[,c(colnames(sce_rest_qc),colnames(sce_34_qc))]
sce
# cluster
sce <- splitAltExps(sce, is.spike)
altExpNames(sce) <- 'spikes'
clust_sce <- quickCluster(sce, min.size=10)
table(clust_sce)
sce <- computeSumFactors(sce, clusters = clust_sce)
sizeFactors(sce)
plot(sizeFactors(sce),sce$stats$sum/1e6, log='xy')
sce <- logNormCounts(sce)
assayNames(sce)
sce <- runPCA(sce)
sce <- runTSNE(sce,perplexity=10)
sce <- runUMAP(sce)
reducedDimNames(sce)
cluster_all_plot <- gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by = 'Cell_Type')
  +ggtitle('PCA'),
  plotReducedDim(sce, "TSNE", colour_by = 'Cell_Type')
  +ggtitle('TSNE'),
  plotReducedDim(sce, "UMAP", colour_by = 'Cell_Type')
  +ggtitle('UMAP'),
  ncol=1)
ggsave('cluster_all.jpg',cluster_all_plot, device='jpg', width = 6, height = 12)

# top 50 HVG
var.out <- modelGeneVarWithSpikes(sce_34_qc,'spikes')
head(var.out)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab='mean log-expression',
     ylab="Total Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
var.fit <- metadata(var.out)
points(var.fit$mean, var.fit$var, col="red", pch=16)
var.out[order(var.out$bio, decreasing=TRUE),] 
hvg <- getTopHVGs(var.out,n=50, var.threshold = 0.5, fdr.threshold = 0.2)
plotExpression(sce_34_qc,hvg)

# Spearman Corr
corr <- correlatePairs(sce_34_qc,subset.row=hvg)
#corr <- corr[corr$FDR < 0.01,]
corr
# Corr heatmap
m.hvg <- matrix(rep(1,50*50), 50, 50)
colnames(m.hvg) <- hvg
rownames(m.hvg) <- hvg
for(i in 1:50){
  for(j in 1:50){
    if(j==i){
      next
    }
    else{
      rho <- corr$rho[(corr$gene1==hvg[i] | corr$gene2==hvg[i]) & 
               (corr$gene1==hvg[j] | corr$gene2==hvg[j])]
      m.hvg[i,j]=rho
    }
  }
}
heatmap_corr <- pheatmap(m.hvg,
                         main = 'Spearman correlation of top 50 HVG')
ggsave('heatmap_corr.jpg',heatmap_corr, device='jpg', width = 10, height = 10)
# dot plot
x_value <- assays(sce_34_qc)$logcounts['Bex3',]
y_value <- assays(sce_34_qc)$logcounts['Tceal9',]
fdr <- corr$FDR[(corr$gene1=='Bex3' | corr$gene2=='Bex3') & 
                  (corr$gene1=='Tceal9' | corr$gene2=='Tceal9')]
png('dot_plot.png')
plot(x_value,y_value, pch=16, col='blue',
     xlab=paste('Bex3','log-normalized expression'),
     ylab=paste('Tceal9','log-normalized expression'))
mtext(side=3, line=0, adj=0, cex=1.2,  
      paste('rho=',round(m.hvg['Bex3','Tceal9'],3),' fdr=', format(fdr, scientific=TRUE)))
dev.off()
# logcounts heatmap
hm_logcount <- pheatmap(assays(sce_34_qc)$logcounts[hvg,],
         annotation_col = data.frame(colData(sce_34_qc)['Cell_Type']),
         main = 'Log count of top 50 HVG')
ggsave('heatmap_logcount.jpg',hm_logcount,device='jpg', width = 10, height = 10)

# Oct4 corr
get_corr <- function(name1,name2,sce){
  if(name1==name2){
    same <- t(as.matrix(c(1,0,0)))
    colnames(same) <- c('rho','p.value','FDR')
    rownames(same) <- name1
    return(same)
  }
  else{
    r <- as.matrix(correlatePairs(sce,subset.row=c(name1,name2))[c('rho','p.value','FDR')])
    rownames(r) <- name1
    return(r)
  }
}
none_zero_genes <- rownames(sce_34_qc)[nexprs(sce_34_qc, byrow=TRUE)>10]
oct4_corr <- lapply(none_zero_genes, FUN = get_corr, name2='Pou5f1', sce=sce_34_qc)
oct4_corr <- do.call(rbind, oct4_corr)
head(oct4_corr)
sum(is.na(oct4_corr))
#oct4_corr <- na.omit(oct4_corr)
oct4_corr <- oct4_corr[rownames(oct4_corr)!='Pou5f1',]
dim(oct4_corr)
oct4_corr_posi <- oct4_corr[oct4_corr[,1]>0,]
oct4_corr_posi <- oct4_corr_posi[order(oct4_corr_posi[,3], decreasing = FALSE),]
head(oct4_corr_posi)
dim(oct4_corr_posi)
oct4_corr_nega <- oct4_corr[oct4_corr[,1]<0,]
oct4_corr_nega <- oct4_corr_nega[order(oct4_corr_nega[,3], decreasing = FALSE),]
head(oct4_corr_nega)
dim(oct4_corr_nega)
oct4_corr_posi_genes <- rownames(oct4_corr_posi)[1:50]
oct4_corr_nega_genes <- rownames(oct4_corr_nega)[1:50]
selected_genes_posi <- c('Pou5f1', oct4_corr_posi_genes)
selected_genes_nega <- c('Pou5f1', oct4_corr_nega_genes)
sorted_cells <- names(sort(assays(sce_34_qc)$logcounts['Pou5f1',]))
heatmap_posi <- pheatmap(assays(sce_34_qc)$logcounts[selected_genes_posi,sorted_cells],
         cluster_row=FALSE,cluster_cols=FALSE,scale='row',
         annotation_col = data.frame(colData(sce_34_qc)['Cell_Type']),
         main = 'Log counts of top 50 positive correlated genes')
ggsave('heatmap_posi.jpg',heatmap_posi,device='jpg', width = 10, height = 10)
heatmap_nega <- pheatmap(assays(sce_34_qc)$logcounts[selected_genes_nega,sorted_cells],
         cluster_row=FALSE,cluster_cols=FALSE,scale='row',
         annotation_col = data.frame(colData(sce_34_qc)['Cell_Type']),
         main = 'Log counts of top 50 negative correlated genes')
ggsave('heatmap_nega.jpg',heatmap_nega,device='jpg', width = 10, height = 10)
