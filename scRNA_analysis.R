library(SingleCellExperiment)
library(scater)
library(scran)
library(AnnotationHub)
library(biomaRt)
library(pheatmap)
library(data.table)
library(gplots)
library(gridExtra)
library(latex2exp)
library(tidyr)
library(AnnotationDbi)
library(scry)
library(scDblFinder)
library(ggVennDiagram)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(ggpubr)
library(PCAtools)
library(intrinsicDimension)
library(bluster)
library(Seurat)
library(SC3)
library(cluster)
library(dplyr)
library(viridisLite)
library(viridis)
library(sctransform)
library(glmGamPoi)
library(parallel)
library(MAST)
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# preparetion ##################################################################
# set working directory
setwd('/home/s2321661/Dissertation')

# set random seed
seed <- 1000


# load data ####################################################################
# read scRNA seq data (raw count) and cell annotation
raw_count <- fread('counts.txt', skip=1, header=TRUE, data.table = FALSE)
rownames(raw_count) <- raw_count$Geneid
raw_count <- raw_count[,7:ncol(raw_count)]
raw_count <- as.matrix(raw_count)
anno <- read.csv('annotation.csv')
idx <- match(colnames(raw_count),anno$Dir)
colnames(raw_count) <- anno$Source.Name[idx]
anno <- anno[idx,]
coldata <- anno[,c('Cell.Type','Source.Name')]
rownames(coldata) <- anno$Source.Name
dim(raw_count)
sce <- SingleCellExperiment(assay=list(counts = raw_count),colData=coldata)
sce
sce_origin <- sce # backup
rm(idx)
rm(anno)
rm(coldata)
rm(raw_count)

# gene name annotation
ah <- AnnotationHub()
#display(ah)  # search 'GRCm39' and 'EnsDb' to find v109 ID
ens.mm.v109 <- AnnotationHub()[["AH109655"]]
columns(ens.mm.v109)
gene_anno <- AnnotationDbi::select(ens.mm.v109, keys=rownames(sce), 
                                   keytype='GENEID', column='GENENAME')
head(gene_anno)
sum(gene_anno$GENENAME == '')  # the number of failures in gene name annotation
rowData(sce)$GENENAME <- make.names(gene_anno$GENENAME[match(rownames(sce),gene_anno$GENEID)],
                                    unique = TRUE)
sum(grepl('^NA',rowData(sce)$GENENAME))  
sum(grepl('^X\\.|^X$',rowData(sce)$GENENAME) == 1)  # the number of failures in gene name annotation
rowData(sce)$ENSEMBL <- rownames(sce)

# ERCC name annotation
is.spike <- grepl("^ERCC", rowData(sce)$ENSEMBL)
sum(is.spike)
idx <-  which(is.spike)
rowData(sce)$GENENAME[idx] <- rowData(sce)$ENSEMBL[idx]
rowData(sce)$GENENAME[idx]

# change sce rownames from ENSEMBL IDs to gene names and ERCC names
rownames(sce) <- rowData(sce)$GENENAME
sce
head(rowData(sce))
length(grep('^NA',rownames(sce)))
sum(colSums(assays(sce[is.spike, ])$counts) == 0)  # check how many cells do not have spike-ins
sce_origin <- sce  #backup
rm(idx)
rm(ah)
rm(gene_anno)
save(list = ls(), file = 'sce_preperation.RData')
#load("sce_preperation.RData")

# Get Mitochondrial genes' ENSEMBL ID from Biomart
#listEnsembl()
#ensembl <- useEnsembl(biomart = "genes")
#head(listDatasets(ensembl))
#searchDatasets(mart = ensembl, pattern = "GRCm39")
#ensembl <- useDataset(dataset='mmusculus_gene_ensembl', mart=ensembl)
#head(listFilters(ensembl))
#head(listAttributes(ensembl))
#MT_genes <- getBM(filters = 'chromosome_name', values='MT',
#                  attributes = c('ensembl_gene_id','external_gene_name'), 
#                  mart = ensembl)
#head(MT_genes)
#MT_id <- paste(MT_genes$ensembl_gene_id, collapse='|')
#head(MT_id)
#is.mito <- grepl(MT_id, rowData(sce)$ENSEMBL)
#count(is.mito)
#rm(ensembl)
#rm(MT_id)

# get mitochondrial genes from annotation files got from Biomart Website
MT_genes <- read.csv("MT_genes.csv")
head(MT_genes)
MT_id <- paste(MT_genes$Gene.stable.ID, collapse='|')
head(MT_id)
is.mito <- grepl(MT_id, rowData(sce)$ENSEMBL)
sum(is.mito)
rm(MT_id)

# cell QC ######################################################################
QC_stats.cell <- perCellQCMetrics(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
# or we can store the stats directly to colData 
# by addPerCellQC(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
QC_stats.cell
sce$total_counts <- QC_stats.cell$sum/1e6
sce$detected_genes <- QC_stats.cell$detected
sce$mito_percent <- QC_stats.cell$subsets_Mt_percent
sce$ERCC_percent <- QC_stats.cell$subsets_ERCC_percent

# check QC failed cells (filter cells)
cell_filter <- perCellQCFilters(QC_stats.cell,sub.fields=c("subsets_Mt_percent",
                                                           "subsets_ERCC_percent"),
                                nmad=5)
# or use quickPerCellQC(QC_stats.cell,sub.fields=c("subsets_Mt_percent", "subsets_ERCC_percent"))
cell_filter  
colSums(as.matrix(cell_filter))  # check how many cells were dropped 
summary(cell_filter$discard)  # check total discarded cells count
sce$discard <- cell_filter$discard
sce$low_lib_size <- cell_filter$low_lib_size
sce$low_detected_genes <- cell_filter$low_n_features
sce$high_mito_percent <- cell_filter$high_subsets_Mt_percent
sce$high_ERCC_percent <- cell_filter$high_subsets_ERCC_percent
sce_origin <- sce  # backup
sce

# Scatter plots that summarize cell QC
thres_total_counts <- attributes(cell_filter$low_lib_size)$thresholds[1]/1e6
thres_detected_genes <- attributes(cell_filter$low_n_features)$thresholds[1]
thres_mito_percent <- attributes(cell_filter$high_subsets_Mt_percent)$thresholds[2]
thres_ERCC_percent <- attributes(cell_filter$high_subsets_ERCC_percent)$thresholds[2]
QC_scatter <- gridExtra::grid.arrange(
  plotColData(sce, x="Cell.Type", y="total_counts", colour_by="discard") +
    ggtitle("Library sizes") + ylab('Total Reads Count') + xlab('Cell Type') + 
    geom_hline(yintercept = thres_total_counts,
               linetype = "dashed", color = "red") + 
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(labels = c(0,
                                  sapply(sort(seq(2,max(sce$total_counts),2)), 
                                         function(x){latex2exp::TeX(paste0('$',x,'\\times 10^6$'))})),
                       breaks = sort(c(0,seq(2,max(sce$total_counts),2)))) +
    annotate(geom='text', x = 3.6, y = thres_total_counts,
             label = latex2exp::TeX(paste0('$',round(thres_total_counts,2),'\\times 10^6$'),
                                    output = 'character'), 
             vjust = 1.2, size=3.9, color='red',parse = TRUE) + 
    coord_cartesian(clip = 'off'),
  plotColData(sce, x="Cell.Type", y="detected_genes", colour_by="discard")+
    ggtitle("Detected genes") + ylab('Number of detected genes') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_detected_genes,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce$detected_genes),2000)) +
    annotate('text',x = 3.5, y = thres_detected_genes,
             label = round(thres_detected_genes,2), 
             vjust = 1.7, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce, x="Cell.Type", y="mito_percent", colour_by="discard")+
    ggtitle("Mito percent") + ylab('Mitochondrial proportion (%)') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_mito_percent,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce$mito_percent),10)) +
    annotate('text',x = 3.5, y = thres_mito_percent,
             label = round(thres_mito_percent,2), 
             vjust = -0.8, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce, x="Cell.Type", y="ERCC_percent", colour_by="discard")+
    ggtitle("ERCC percent") + ylab('ERCC proportion (%)') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_ERCC_percent,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce$ERCC_percent),10)) +
    annotate('text',x = 3.5, y = thres_ERCC_percent,
             label = round(thres_ERCC_percent,2), 
             vjust = -0.8, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  ncol=2,
  nrow=2
)
ggsave('figures/QC_summary.pdf',QC_scatter, device='pdf', width = 15, height = 15)

# Table summary of cell QC
total.cells <- table(sce$Cell.Type)
total.cells <- c(total.cells,Total=sum(total.cells))
# lib size
low_lib_size.cells <- as.data.frame(colData(sce)) %>% group_by(Cell.Type) %>% 
  summarise(count=sum(low_lib_size==TRUE)) %>% tibble::deframe()
low_lib_size.cells <- c(low_lib_size.cells, Total=sum(low_lib_size.cells))
low_lib_size.prop <- round(low_lib_size.cells/total.cells,4)*100
# detected genes
low_detected_genes.cells <- as.data.frame(colData(sce)) %>% group_by(Cell.Type) %>% 
  summarise(count=sum(low_detected_genes==TRUE)) %>% tibble::deframe()
low_detected_genes.cells <- c(low_detected_genes.cells, Total=sum(low_detected_genes.cells))
low_detected_genes.prop <- round(low_detected_genes.cells/total.cells,4)*100
# mito percent
high_mito_percent.cells <- as.data.frame(colData(sce)) %>% group_by(Cell.Type) %>% 
  summarise(count=sum(high_mito_percent==TRUE)) %>% tibble::deframe()
high_mito_percent.cells <- c(high_mito_percent.cells, Total=sum(high_mito_percent.cells))
high_mito_percent.prop <- round(high_mito_percent.cells/total.cells,4)*100
# ERCC percent
high_ERCC_percent.cells <- as.data.frame(colData(sce)) %>% group_by(Cell.Type) %>% 
  summarise(count=sum(high_ERCC_percent==TRUE)) %>% tibble::deframe()
high_ERCC_percent.cells <- c(high_ERCC_percent.cells, Total=sum(high_ERCC_percent.cells))
high_ERCC_percent.prop <- round(high_ERCC_percent.cells/total.cells,4)*100
# discard
discard.cells <- as.data.frame(colData(sce)) %>% group_by(Cell.Type) %>% 
  summarise(count=sum(discard==TRUE)) %>% tibble::deframe()
discard.cells <- c(discard.cells, Total=sum(discard.cells))
discard.prop <- round(discard.cells/total.cells,4)*100
# keep
keep.cells <- as.data.frame(colData(sce)) %>% group_by(Cell.Type) %>% 
  summarise(count=sum(discard==FALSE)) %>% tibble::deframe()
keep.cells <- c(keep.cells, Total=sum(keep.cells))
keep.prop <- round(keep.cells/total.cells,4)*100
# create table
QC_table <- cbind(total.cells, 
                  low_lib_size.cells,
                  low_lib_size.prop,
                  low_detected_genes.cells,
                  low_detected_genes.prop,
                  high_mito_percent.cells,
                  high_mito_percent.prop,
                  high_ERCC_percent.cells,
                  high_ERCC_percent.prop,
                  discard.cells, 
                  discard.prop,
                  keep.cells,
                  keep.prop)
QC_table <- as.data.frame(QC_table)
QC_table
# lib size
QC_table['low_lib_size.prop'] <- sapply(QC_table['low_lib_size.prop'],
                                        function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, low_lib_size, 
                         c(low_lib_size.cells, low_lib_size.prop), sep=' ')
# detected genes
QC_table['low_detected_genes.prop'] <- sapply(QC_table['low_detected_genes.prop'],
                                              function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, low_detected_genes, 
                         c(low_detected_genes.cells, low_detected_genes.prop), sep=' ')
# mito percent
QC_table['high_mito_percent.prop'] <- sapply(QC_table['high_mito_percent.prop'],
                                             function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, high_mito_percent, 
                         c(high_mito_percent.cells, high_mito_percent.prop), sep=' ')
# ERCC percent
QC_table['high_ERCC_percent.prop'] <- sapply(QC_table['high_ERCC_percent.prop'],
                                             function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, high_ERCC_percent, 
                         c(high_ERCC_percent.cells, high_ERCC_percent.prop), sep=' ')
# discard
QC_table['discard.prop'] <- sapply(QC_table['discard.prop'],
                                   function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, discard, 
                         c(discard.cells, discard.prop), sep=' ')
# keep 
QC_table['keep.prop'] <- sapply(QC_table['keep.prop'],
                                function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, keep, 
                         c(keep.cells, keep.prop), sep=' ')
# add index column
QC_table <- cbind(rownames(QC_table),QC_table)
# change column names
colnames(QC_table) <- c('Cell Types', 
                        'Total', 
                        'Low Library Size',
                        'Low Detected Genes',
                        'High Mito Percent',
                        'High ERCC Percent',
                        'Total Discard', 
                        'Keep')
QC_table
write.csv(QC_table ,'figures/QC_table.csv', row.names = FALSE)
rm(total.cells)
rm(low_lib_size.cells)
rm(low_lib_size.prop)
rm(low_detected_genes.cells)
rm(low_detected_genes.prop)
rm(high_mito_percent.cells)
rm(high_mito_percent.prop)
rm(high_ERCC_percent.cells)
rm(high_ERCC_percent.prop)
rm(discard.cells)
rm(discard.prop)
rm(keep.cells)
rm(keep.prop)

# Venn Diagram summary of failures in cell QC
venn <- as.data.frame(cell_filter)
QC_venn <- ggVennDiagram(apply(venn[1:4], 2, function(x) which(x == TRUE)),
                         label_alpha=0, 
                         set_color = c("deepskyblue2","darkolivegreen","darkorange","darkorchid"), 
                         category.names = c(paste("Low Library Size (",sum(venn[1]),")",sep = ''),
                                            paste("Low Detected Genes (",sum(venn[2]),")",sep = ''),
                                            paste("High Mito Percent (",sum(venn[3]),")",sep = ''),
                                            paste("High ERCC percent (",sum(venn[4]),")",sep = ''))) + 
  scale_fill_gradient(low="white",high = "coral1") +
  ggtitle('Summary of failures in cell QC') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_color_manual(values = c("deepskyblue2","darkolivegreen","darkorange","darkorchid")) +
  annotate('table',label = QC_table[,c(1,2,7,8)],x = 0, y = 0,vjust = 0.35, hjust = -0.4) +
  theme(legend.position = c(0.95, 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
QC_venn
ggsave('figures/QC_venn.jpg',QC_venn, device='jpg', width = 8, height = 8)
rm(venn)

# Histogram summary of cell QC
QC_hist <- gridExtra::grid.arrange(
  ggplot(as.data.frame(QC_stats.cell), aes(x=sum/1e6))+
    geom_histogram(color='grey80',bins=50) +
    geom_density(alpha=.2, fill="blue") +
    xlab('Library sizes (millions)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_total_counts), color = 'red', linetype="dashed") +
    annotate('text',x=thres_total_counts+0.8,y=100, label=round(thres_total_counts,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(QC_stats.cell), aes(x=detected))+
    geom_histogram(color='grey80',bins=50) +
    xlab('Number of detected genes') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_detected_genes), color = 'red', linetype="dashed") +
    annotate('text',x=thres_detected_genes+1400,y=60, label=round(thres_detected_genes,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(QC_stats.cell), aes(x=subsets_Mt_percent))+
    geom_histogram(color='grey80',bins=50) +
    xlab('Mitochondrial proportion (%)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_mito_percent), color = 'red', linetype="dashed") +
    annotate('text',x=thres_mito_percent+5,y=200, label=round(thres_mito_percent,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(QC_stats.cell), aes(x=subsets_ERCC_percent))+
    geom_histogram(color='grey80',bins=50) +
    xlab('ERCC proportion (%)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_ERCC_percent), color = 'red', linetype="dashed") +
    annotate('text',x=thres_ERCC_percent+1.5,y=300, label=round(thres_ERCC_percent,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ncol=2
)
ggsave('figures/QC_hist.jpg',QC_hist, device='jpg', width = 10, height = 8)

# drop cells 
# cell QC fail
sce <- sce[, !cell_filter$discard]

# remove cells with no expression of Oct4 and Sox2 (cell QC 2)
# check how many cells are not expressing Oct4 in each cell type
table(sce$Cell.Type)
Oct4_fail <- c(
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'E3.5 ICM'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'E4.5 Epiblast'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'E5.5 Epiblast'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'Epi'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'ES'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'ES2i'] == 0)
)
Oct4_fail
# check how many cells are not expressing Sox2 in each cell type
Sox2_fail <- c(
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'E3.5 ICM'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'E4.5 Epiblast'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'E5.5 Epiblast'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'Epi'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'ES'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'ES2i'] == 0)
)
Sox2_fail
# summary
Cell_Types <- c('E3.5 ICM','E4.5 Epiblast','E5.5 Epiblast', 'Epi', 'ES', 'ES2i')
Total <- c(
  sum(colData(sce)$Cell.Type == 'E3.5 ICM'),
  sum(colData(sce)$Cell.Type == 'E4.5 Epiblast'),
  sum(colData(sce)$Cell.Type == 'E5.5 Epiblast'),
  sum(colData(sce)$Cell.Type == 'Epi'),
  sum(colData(sce)$Cell.Type == 'ES'),
  sum(colData(sce)$Cell.Type == 'ES2i')
)
QC_table_Oct4_Sox2 <- data.frame(Cell_Types, Total, Oct4_fail, Sox2_fail)
QC_table_Oct4_Sox2
write.csv(QC_table_Oct4_Sox2 ,'figures/QC_table_Oct4_Sox2.csv', row.names = FALSE)
rm(Oct4_fail)
rm(Sox2_fail)
rm(Cell_Types)
rm(Total)


# gene filtering ###############################################################
# remove spikes
sce <- splitAltExps(sce, is.spike)
altExpNames(sce) <- 'spikes'
sce

# average drop-out rate
mean(rowSums(counts(sce)==0)/ncol(sce))
ggplot(data.frame(numbers=1:nrow(sce),
                  drop_out_rate=rowSums(counts(sce)==0)/ncol(sce))) +
  geom_histogram(aes(x=drop_out_rate),color='grey80',bins=100)

# remove low expression genes (express in less than 3 cells)
QC_stats.gene <- perFeatureQCMetrics(sce)
detected_cell_prop.hist <- ggplot(as.data.frame(QC_stats.gene), aes(x=detected)) +
  geom_histogram(color='grey80',bins=100) +
  xlab('Detected %') +
  ylab('Counts') + 
  ggtitle('Genes with expression >0')
detected_cell_prop.hist
ggsave('figures/detected_cell_prop.jpg',detected_cell_prop.hist, 
       device='jpg', width = 8, height = 6)
rowData(sce)$mean <- QC_stats.gene$mean
rowData(sce)$detected_prop <- QC_stats.gene$detected
sum(rowSums(assays(sce)$counts > 0) >= 3)  # check how many genes are left 
# remove genes
sce <- sce[rowSums(assays(sce)$counts > 0) >= 3,]
sce


# Normalisation ################################################################
sce <- computeSumFactors(sce, cluster = quickCluster(sce))
sce <- logNormCounts(sce)
sce  # assays(sce) has 'logcounts'


# check oct4 variation #########################################################
# log counts
oct4_variation_logcounts <- plotExpression(sce,'Pou5f1',x='Cell.Type',exprs_values = "logcounts") +
  xlab('Cell Type') + ylab('Log Counts') + 
  ggtitle('Oct4 expression level after Normalisation (log counts)')
oct4_variation_logcounts
ggsave('figures/oct4_variation_logcounts.jpg',oct4_variation_logcounts, device='jpg', width = 8, height = 10)


# famous correlated & non-correlated genes(scatter plot) #######################
# positive control
targets <- c('Nanog', 'Sox2','Klf4','Zfp42','Utf1','Esrrb')
#rownames(sce)[grepl('Esrrb',rownames(sce))]

# log counts
oct4_corr_scatter_logcounts <- plotExpression(sce,targets,x='Pou5f1',exprs_values = "logcounts",color_by = 'Cell.Type') + 
  xlab('Oct4 expression level') + 
  ylab('Expression level (log counts after normalisation)') +
  ggtitle('Expression level of Oct4 and some of its target genes (after normalisation)') +
  guides(color=guide_legend("Cell Type"))
oct4_corr_scatter_logcounts
ggsave('figures/oct4_corr_scatter_logcounts.jpg',oct4_corr_scatter_logcounts, device='jpg', width = 8, height = 10)
rm(targets)

# negative control
hk_genes <- c('Actb','Tbp','Pgk1','Ppia','Rpl38','Hmbs')

# log counts
oct4_corr_scatter_logcounts.hk <- plotExpression(sce,hk_genes,x='Pou5f1',exprs_values = "logcounts",color_by = 'Cell.Type') + 
  xlab('Oct4 expression level') + 
  ylab('Expression level (log counts after normalisation)') +
  ggtitle('Expression level of Oct4 and some housekeeping genes (after normalisation)') +
  guides(color=guide_legend("Cell Type"))
oct4_corr_scatter_logcounts.hk
ggsave('figures/oct4_corr_scatter_logcounts_hk.jpg',oct4_corr_scatter_logcounts.hk, device='jpg', width = 8, height = 10)

# Deviance (select genes that are both highly expressed and highly variable)
dev <- rowData(devianceFeatureSelection(sce,assay = "counts", fam='binomial'))$binomial_deviance

# log counts of deviance
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:100]
plotExpression(sce,features_dev,exprs_values = "logcounts")+ ylim(0,20)
ggsave('figures/top100_dev_log.jpg',device='jpg', width = 20, height = 6)


# cluster cells (deviance feature) #############################################
# cluster by seurat package (using Leiden/Louvain algorithm)
#features_hvg <- getTopHVGs(var.out,n=5000)
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:5000]
sce_seurat <- sce_origin
sce_seurat <- sce_seurat[, !cell_filter$discard]
sce_seurat <- sce_seurat[!is.spike,]
sce_seurat <- sce_seurat[rowSums(assays(sce_seurat)$counts > 0) >= 3,]
sce_seurat
sce_seurat <- computeSumFactors(sce_seurat, cluster = quickCluster(sce_seurat))
sce_seurat <- logNormCounts(sce_seurat)
seurat <- as.Seurat(sce_seurat, counts = "counts", data = "logcounts")
#seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
#VariableFeatures(seurat)[1:100]
#VlnPlot(object = seurat, features = c('Pou5f1','Nanog'))
seurat <- ScaleData(seurat, features = rownames(seurat))

# select the best PC number/selection method
feat_num_sele <- 4000
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:feat_num_sele]
seurat <- RunPCA(seurat, features = features_dev, seed.use = seed)  # here the features are chosen by deviance
#DimHeatmap(seurat, reduction = "pca",dims = 1:3)
PC_num.elbow <- findElbowPoint(seurat[['pca']]@stdev)
PC_num.elbow
PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
PC_num.gml
PC_num.gml <- round(PC_num.gml$dim.est,0)
PC_num.gml
# check if these 2 methods provide convincing results
ElbowPlot(seurat,ndims = 30) +
  geom_vline(xintercept = PC_num.elbow, color = 'blue', linetype="dashed") +
  geom_vline(xintercept = PC_num.gml, color = 'coral', linetype="dashed") +
  annotate('text',x=PC_num.elbow, y=10,label='elbow point',color = 'red') +
  annotate('text',x=PC_num.gml, y=12,label='global maximum likelihood',color = 'red') +
  ggtitle('Variation represented by each PC (feature number: 4000)')
ggsave('figures/PC_num_sele.jpg',device='jpg', width = 10, height = 6)
# identify significant PC (JackStraw Plot)
feat_num_sele <- 4000  # can change from 1000-5000 and check results
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:feat_num_sele]
seurat <- RunPCA(seurat, features = features_dev, seed.use = seed)  # here the features are chosen by deviance
seurat <- JackStraw(seurat, num.replicate = 100, dims = 30)
seurat <- ScoreJackStraw(seurat, dims = 1:30)
JackStrawPlot(seurat, dims = 1:30) +
  ggtitle('JackStraw Plot (fearture number: 4000)')  # find sharp drop-off in p-value
ggsave('figures/PC_JSP_4000.jpg',device='jpg', width = 10, height = 6)
# plot p-value changes
PC_pvalue <- -log10(seurat@reductions[["pca"]]@jackstraw@overall.p.values[,'Score'])
PC_pvalue[is.infinite(PC_pvalue)] <- 200
ggplot(data.frame(PC_num=1:length(PC_pvalue), log10_pvalue=PC_pvalue)) +
  geom_point(aes(y=log10_pvalue,x=PC_num)) +
  geom_line(aes(y=log10_pvalue,x=PC_num)) +
  ylab('-log10 P-value') + xlab('PC number') +
  ggtitle('P-value changes in JackStraw analysis (feature number: 4000)')
ggsave('figures/PC_JSP_4000_pvalue.jpg',device='jpg', width = 10, height = 6)
# PC numbers selected by GML/elbow are not satisfying here, so set to 14 based on JackStraw Plot
PC_num.gml <- 14
rm(PC_pvalue)

# test for resolution of leiden and features number 
sweep_para <- function(seurat, seed, features, res, feature_num, PC_method=1, set_PC_num=15){
  # PC number selection can be changed!!!! (Here I change to 14 based on JackStraw plot)
  # PC_method: 1='GML', 2='elbow', 3='JS'
  para_res <- c()
  sil_score <- c()
  cluster_num <- c()
  ft_num <- c()
  PC_num <- c()
  for(j in feature_num){
    print(paste0('feature number = ', j))
    seurat <- RunPCA(seurat, features = features[1:j], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.elbow <- findElbowPoint(seurat[['pca']]@stdev)
    PC_sele <- switch(PC_method,round(PC_num.gml$dim.est,0),PC_num.elbow,set_PC_num)
    seurat <- FindNeighbors(seurat, dims = 1:PC_sele)
    for(i in res){
      print(paste0('res = ', i))
      para_res <- c(para_res,i)
      seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
      cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
      sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
      sil.data <- as.data.frame(sil)
      score <- mean(sil.data$sil_width)
      sil_score <- c(sil_score,score)
      cluster_num <- c(cluster_num, length(table(cluster)))
      ft_num <- c(ft_num, j)
      PC_num <- c(PC_num, PC_sele)
    }
    cat('\n')
  }
  r <- data.frame(res=para_res, sil_score=sil_score, clust_num=cluster_num,
                  feature_num=ft_num, PC_num=PC_num)
  return(r)
}
out <- sweep_para(seurat, seed, features_dev, 
                  res=seq(0.1,2.0,0.1),
                  feature_num = c(seq(1000,5000,1000)),
                  3,PC_num.gml)
out
ggplot(out, aes(x=res,y=sil_score,label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_line(aes(color = factor(feature_num))) +
  ylab('Silhouette Score') +
  xlab('Leiden Resolution') +
  guides(color=guide_legend(title="Feature Number")) +
  scale_color_brewer(palette="Paired")
ggsave('figures/leiden_para_1_dev.jpg',device='jpg', width = 8, height = 6)
ggplot(out, aes(x = factor(feature_num), y = sil_score, 
                label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_bar(stat = "summary", fun = "var")+
  ylab('Variance of Silhouette') +
  xlab('Feature Number')
ggsave('figures/leiden_para_2_dev.jpg',device='jpg', width = 8, height = 5)
rm(out)

# comparison of different cluster results (different resolutions) but with same cluster number
compare_para <- function(seurat, seed, resolution, feature_num){
  clusters <- list()
  for(i in c(1,2)){
    res <- resolution[i]
    ft_num <- feature_num[i]
    seurat <- RunPCA(seurat, features = features_dev[1:ft_num], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
    seurat <- FindClusters(seurat, resolution = res, algorithm = 4)
    name <- paste0('originalexp_snn_res.',res)
    cluster <- seurat@meta.data[[name]]
    clusters[[i]] <- cluster
  }
  return(clusters)
}
com_para_r <- compare_para(seurat, seed, c(1.1,1.1),c(2000, 4000))
table(com_para_r[[1]],com_para_r[[2]])
rm(com_para_r)

# select the best feature number and PC number
seurat <- RunTSNE(seurat, seed.use = seed, perplexity = 20, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "tsne" , group.by = 'Cell.Type')
#ggsave('figures/leiden_tsne.jpg', device='jpg', width = 8, height = 7)
seurat <- RunUMAP(seurat,seed.use = seed, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "umap" , group.by = 'Cell.Type')

#  convey those dim-reduction plots to a SCE object (to use some functions in bioconductor)
sce_dev <- as.SingleCellExperiment(seurat)
plotReducedDim(sce_dev, "TSNE", colour_by="Cell.Type")  # make sure it's the same as seurat's tSNE
counts(sce_dev) <- as.matrix(counts(sce_dev))
assays(sce_dev)$norm_counts <- 2^(logcounts(sce_dev))-1
sce_dev

# select the best resolution and cluster cells using Leiden algorithm
res_sele <-1.1  # 0.1, 0.4, 0.7, 1.1, 1.2
seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
seurat <- FindClusters(seurat, resolution = res_sele, algorithm = 4)  # algorithm 4 is "Leiden"; 1 is "Louvain"
DimPlot(seurat, reduction = "tsne", label = TRUE, shape.by = 'Cell.Type')

# # use cluster results by SC3 (7 clusters)
# seurat$sc3_7_clusters <- sce_dev_SC3$sc3_7_clusters
# Idents(seurat) <- sce_dev_SC3$sc3_7_clusters
# colLabels(sce_dev) <- sce_dev_SC3$sc3_7_clusters

# pass the clustering result back to SCE object
colLabels(sce_dev) <- seurat@meta.data[[paste0('originalexp_snn_res.',res_sele)]]

# find marker genes of each cluster (heatmap)
group_marker_wilcox <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25, test.use = 'wilcox')
group_marker_wilcox <- group_marker_wilcox %>% dplyr::filter(p_val_adj < 0.01)
group_marker_mast <- FindAllMarkers(seurat, only.pos = TRUE, test.use = 'MAST')
group_marker_mast <- group_marker_mast %>% dplyr::filter(p_val_adj < 0.01)

# # find markers intersection in 2i groups (2 clusters) (intersection of wilcox and mast)
# ggVennDiagram(list(wil_clust4 = group_marker_wilcox[group_marker_wilcox$cluster == 4,'gene'], 
#                    wil_clust6 = group_marker_wilcox[group_marker_wilcox$cluster == 6,'gene'], 
#                    mast_clust4 = group_marker_mast[group_marker_mast$cluster == 4,'gene'], 
#                    mast_clust6 = group_marker_mast[group_marker_mast$cluster == 6,'gene']),
#               label = "count", label_size = 5) +
#   theme(legend.position = "none") +
#   ggtitle("intersections of markers (FDR < 0.01) in Kim's and Chen's dataset (2i)") +
#   theme(plot.title = element_text(hjust = 0.5))
# find_key_markers <- function(a1,a2,b1,b2){
#   a1a2 <- intersect(a1,a2)
#   a1 <- a1[! a1 %in% a1a2]
#   a2 <- a2[! a2 %in% a1a2]
#   b1b2 <- intersect(b1,b2)
#   b1 <- b1[! b1 %in% b1b2]
#   b2 <- b2[! b2 %in% b1b2]
#   return(list(a1b1 = intersect(a1,b1), a2b2 = intersect(a2,b2)))
# }
# key_markers <- find_key_markers(group_marker_wilcox[group_marker_wilcox$cluster == 4,'gene'],
#                                 group_marker_wilcox[group_marker_wilcox$cluster == 6,'gene'],
#                                 group_marker_mast[group_marker_mast$cluster == 4,'gene'],
#                                 group_marker_mast[group_marker_mast$cluster == 6,'gene'])
# head(key_markers$a1b1)

write.csv(group_marker_mast, 'leiden_markers_mast.csv')
write.csv(group_marker_wilcox, 'leiden_markers_wilcox.csv')
group_marker <- group_marker_mast
group_marker %>% group_by(cluster) %>% count()

# check markers intersection (Kim' dataset and Chen' dataset: 2i groups and serum groups)
# wilcox test
test_markers <- read.csv('/home/s2321661/Test/figures/leiden_markers_wilcox.csv')
rownames(test_markers) <- test_markers$X
test_markers <- test_markers[,2:ncol(test_markers)]
test_markers %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 100, order_by = p_val_adj) -> test_top100
group_marker_wilcox %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 100, order_by = p_val_adj) -> group_top100
ggVennDiagram(list(Kim_clust1 = test_top100[test_top100$cluster == 1,] %>% .$gene, 
                   Kim_clust4 = test_top100[test_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust4 = group_top100[group_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust6 = group_top100[group_top100$cluster == 6,] %>% .$gene),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle("intersections of top100 markers (wilcox) in Kim's and Chen's dataset (2i)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/markers_2i_inter_Kim&Chen_wilcox.pdf",device='pdf', width = 12, height = 12)
ggVennDiagram(list(Kim_clust3 = test_top100[test_top100$cluster == 1,] %>% .$gene, 
                   Kim_clust5 = test_top100[test_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust5 = group_top100[group_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust8 = group_top100[group_top100$cluster == 6,] %>% .$gene,
                   Chen_clust10 = group_top100[group_top100$cluster == 10,] %>% .$gene),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle("intersections of top 100 markers (wilcox) in Kim's and Chen's dataset (serum)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/markers_serum_inter_Kim&Chen_wilcox.pdf",device='pdf', width = 12, height = 12)

# MAST
test_markers <- read.csv('/home/s2321661/Test/figures/leiden_markers_mast.csv')
rownames(test_markers) <- test_markers$X
test_markers <- test_markers[,2:ncol(test_markers)]
test_markers %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 100, order_by = p_val_adj) -> test_top100
group_marker_mast %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 100, order_by = p_val_adj) -> group_top100
ggVennDiagram(list(Kim_clust1 = test_top100[test_top100$cluster == 1,] %>% .$gene, 
                   Kim_clust4 = test_top100[test_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust4 = group_top100[group_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust6 = group_top100[group_top100$cluster == 6,] %>% .$gene),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle("intersections of top100 markers (MAST) in Kim's and Chen's dataset (2i)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/markers_2i_inter_Kim&Chen_mast.pdf",device='pdf', width = 12, height = 12)
ggVennDiagram(list(Kim_clust3 = test_top100[test_top100$cluster == 1,] %>% .$gene, 
                   Kim_clust5 = test_top100[test_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust5 = group_top100[group_top100$cluster == 4,] %>% .$gene, 
                   Chen_clust8 = group_top100[group_top100$cluster == 6,] %>% .$gene,
                   Chen_clust10 = group_top100[group_top100$cluster == 10,] %>% .$gene),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle("intersections of top 100 markers (MAST) in Kim's and Chen's dataset (serum)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/markers_serum_inter_Kim&Chen_mast.pdf",device='pdf', width = 12, height = 12)
rm(test_markers)
rm(test_top100)
rm(group_top100)

# plot marker genes of each cluster (heatmap)
# top 15 is created in case duplications appear
group_marker %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 15, order_by = p_val_adj) -> top15
top15 <- top15[!duplicated(top15$gene),]
top15 %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 10, order_by = p_val_adj) -> top10
top10
DoHeatmap(seurat,features=top10$gene, slot = "scale.data") #+ scale_fill_virdis()
#FeaturePlot(seurat, reduction = 'TSNE', features=top10[top10['cluster']==1,]$gene[1:4])
leiden_markers <- plotHeatmap(sce_dev, exprs_values = "logcounts", 
                              order_columns_by=c("label", "Cell.Type"), 
                              features=top10$gene, cluster_rows = FALSE, center = TRUE,
                              gaps_col = cumsum(as.numeric(table(colLabels(sce_dev)))),
                              gaps_row = cumsum(top10 %>% group_by(cluster) %>% count() %>% .$n), 
                              color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40))
ggsave('figures/leiden_markers_10clust_mast.jpg', leiden_markers, device='jpg', width = 14, height = 14)

# details of Silhouette width (and tSNE plots of clustering results)
#sil <- cluster::silhouette(as.integer(cluster.ld), dist(reducedDim(sce_dev, "PCA")))
detail_sil <- function(res, seurat, PC_sele){
  for(i in res){
    print(paste0('res = ', i))
    slot_name <- paste0('originalexp_snn_res.',i)
    dir1 <- paste0('figures/leiden_res_',i,'_1.jpg')
    dir2 <- paste0('figures/leiden_res_',i,'_2.jpg')
    dir3 <- paste0('figures/leiden_res_',i,'_tsne.jpg')
    seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
    cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
    sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
    sil.data <- as.data.frame(sil)
    sil.data$closest <- factor(ifelse(sil.data$sil_width > 0, cluster, sil.data$neighbor))
    #sil.data$cluster <- cluster.ld
    ggplot(sil.data, aes(x=cluster, y=sil_width, colour=closest)) +
      ggbeeswarm::geom_quasirandom(method="smiley") + 
      ylab('Silhouette width') +
      geom_hline(aes(yintercept=0), color = 'red', linetype="dashed")
    ggsave(dir1, device='jpg', width = 6, height = 5)
    ggplot(as.data.frame(sil.data), aes(x=sil_width))+
      geom_histogram(color='grey80',bins=50) +
      geom_density(alpha=.2, fill="blue") + xlab('Silhouette Width')+
      ggtitle(paste0('avarage Silhouette Width: ',round(mean(sil.data$sil_width),4))) +
      geom_vline(aes(xintercept=0), color = 'red', linetype="dashed")
    ggsave(dir2, device='jpg', width = 5, height = 4)
    DimPlot(seurat, reduction = "tsne" ,label = TRUE, shape.by = 'Cell.Type')
    ggsave(dir3, device='jpg', width = 8, height = 6)
  }
}
detail_sil(res=c(0.1,0.4,1.1,1.2), seurat = seurat, PC_sele=PC_num.gml) # need to generate cluster data with specific resolution in advance to Seurat object

# clean up variable
rm(group_diff)
rm(group_marker)
rm(top10)
rm(top15)
rm(group_diff_serum_2)
rm(leiden_markers)
rm(up_reg)
rm(down_reg)
rm(sce_seurat)


# cluster by SC3 package (k-mean based)
sce_dev_SC3 <- sce_dev
rowData(sce_dev_SC3)$feature_symbol <- rownames(sce_dev_SC3)
sce_dev_SC3
sce_dev_SC3 <- sc3(sce_dev_SC3, gene_filter=TRUE, ks = 5:12, biology = TRUE, rand_seed=seed,
                   d_region_min=0.01)
#sc3_plot_consensus(sce_dev_SC3, k = 6, show_pdata = c("Cell.Type", "sc3_6_clusters"))

# evaluate the best clustering number for SC3 (check the Silhouette score in picture!)
eva_SC3_cluster <- function(sce,cluster_num){
  for(i in cluster_num){
    filename <- paste0('figures/SC3_sil_',i,'.png')
    png(filename,1000,1000)
    sc3_plot_silhouette(sce, k = i)
    dev.off()
  }
}
eva_SC3_cluster(sce_dev_SC3,5:12)
#sc3_plot_cluster_stability(sce_dev_SC3, k = 6)

# plot cluster result
plotTSNE(sce_dev_SC3,shape_by = 'Cell.Type', colour_by = 'sc3_12_clusters')
ggsave('figures/SC3_tsne_12_clust.jpg', device='jpg', width = 6, height = 4.5)
#sc3_plot_expression(sce_dev_SC3, k = 6,show_pdata = c("sc3_6_clusters","Cell.Type"))

# markers of each cluster
#sc3_plot_de_genes(sce_dev_SC3, k = 6,show_pdata = c("sc3_6_clusters", "Cell.Type"))
SC3_markerplot <- sc3_plot_markers(sce_dev_SC3, k = 7, auroc = 0.8, p.val = 0.01,
                                   show_pdata = c("sc3_7_clusters", "Cell.Type"))
ggsave('figures/SC3_markerplot_7_clust.jpg', SC3_markerplot, device='jpg', width = 10, height = 12)
cluster.sc3 <- colData(sce_dev_SC3)$sc3_12_clusters
colLabels(sce_dev_SC3) <- cluster.sc3
rm(SC3_markerplot)

# # Find Markers by SC3 (not as good as seurat)
# find_marker_sc3 <- function(sce, k, padj, auroc, clust){
#   clust_name <- paste0('sc3_',k,'_markers_clusts')
#   padj_name <- paste0('sc3_',k,'_markers_padj')
#   auroc_name <- paste0('sc3_',k,'_markers_auroc')
#   df <- rowData(sce)[,c(clust_name, padj_name, auroc_name)]
#   df[is.na(df[[clust_name]]) |
#        is.na(df[[padj_name]]) |
#        is.na(df[[auroc_name]]),] <- 0
#   df <- df[df[[clust_name]] == clust &
#              df[[padj_name]] < padj &
#              df[[auroc_name]] > auroc,]
#   df <- df[order(df[[padj_name]], decreasing = FALSE),]
#   markers <- rownames(df)
#   return(markers)
# }
# sc3_markers <- find_marker_sc3(sce_dev_SC3, 6, padj=0.01, auroc=0.7, clust=1)
# sc3_markers
# plotHeatmap(sce_dev_SC3, exprs_values = "logcounts",
#             order_columns_by=c("label", "Cell.Type"), features=c(sc3_markers),
#             cluster_rows = FALSE, center=TRUE)


# Spearman Correlation of all Oct4 potential targets############################
# load potential Oct4 targets
oct4_tar.potent <- read.csv('oct4_targets.csv', na.strings = c("", "NA"))
oct4_tar.potent <- na.omit(oct4_tar.potent[,1:2])
oct4_tar.potent <- oct4_tar.potent$Gene.name
table(oct4_tar.potent %in% rownames(sce_dev))  
oct4_tar.potent[!oct4_tar.potent %in% rownames(sce_dev)]  # check which genes are not available
oct4_tar.potent <- c('Pou5f1',oct4_tar.potent[oct4_tar.potent %in% rownames(sce_dev)])
tar_length <- length(oct4_tar.potent)
tar_length

# calculate Spearman correlation
oct4_tar.corr <- correlatePairs(sce_dev[,sce_dev$Cell.Type == 'ES'], subset.row=oct4_tar.potent)
oct4_tar.corr <- oct4_tar.corr[apply(is.na(oct4_tar.corr),1,sum)==0,]
oct4_tar.corr

# build correlation matrix
build_corr <- function(corr){
  df <- as.data.frame(corr[c('gene1','gene2','rho')])
  df2 <- data.frame(gene1=df$gene2, gene2=df$gene1, rho=df$rho)
  df <- rbind(df,df2)
  corr.matrix <- spread(df, key = gene2, value = rho)
  rownames(corr.matrix) <- corr.matrix$gene1
  corr.matrix$gene1 <- NULL  # remove 1st column
  diag(corr.matrix) <- 1
  return(corr.matrix)
}
oct4_tar_corr.matrix <- build_corr(oct4_tar.corr)

# plot correlation
oct4_potent_tar <- pheatmap(oct4_tar_corr.matrix,fontsize_row=3,fontsize_col=3, 
                            main = 'Spearman correlation of potential Oct4 targets in serum',
                            breaks = seq(-1,1,0.05), 
                            color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(40))
# change Oct4 text color to red
names = oct4_potent_tar$gtable$grobs[[6]]$label
color = rep('black', length(oct4_tar.potent))
color[grep('Pou5f1',names)] <- 'red'
color[grep('Sox2',names)] <- 'red'
color[grep('Nanog',names)] <- 'red'
oct4_potent_tar$gtable$grobs[[5]]$gp=grid::gpar(col=color, fontsize=rep(3,length(oct4_tar.potent)))
oct4_potent_tar$gtable$grobs[[6]]$gp=grid::gpar(col=color, fontsize=rep(3,length(oct4_tar.potent)))
oct4_potent_tar
ggsave('figures/oct4_potent_tar_serum.jpg',oct4_potent_tar, device='jpg', width = 15, height = 15, dpi=400)
rm(names)
rm(color)
rm(tar_length)


# permutation test -> most correlated genes of Oct4 ############################
get_corr <- function(name1,name2,sce){
  if(name1==name2){
    same <- t(as.matrix(c(1,0)))
    colnames(same) <- c('rho','p.value')
    rownames(same) <- name1
    return(same)
  }
  else{
    r <- as.matrix(correlatePairs(sce,subset.row=c(name1,name2), 
                                  assay.type = "logcounts")[c('rho','p.value')])
    rownames(r) <- name1
    return(r)
  }
}

spearman_corr <- function(sce_obj, clust_num, gene_name, ctrl_genes, clust_name, top_n){
  # gene_name is the name of the core gene used for correlation calculation
  # con_genes means control gene names
  sce <- sce_obj[,sce_obj$label %in% clust_num]
  
  # filter genes that are expressed only in less than 20% cells
  cell_num <- round(ncol(sce)/5,0)
  filtered_genes <- rownames(sce)[nexprs(sce, byrow=TRUE) > cell_num]
  print(paste0('Total number of genes: ',length(filtered_genes)))
  
  # calculate corr of Oct4 and all genes
  corr <- mclapply(filtered_genes, FUN = get_corr, name2=gene_name, sce=sce, mc.cores=48)
  corr <- do.call(rbind, corr)
  corr <- as.data.frame(corr)
  corr$FDR <- p.adjust(corr$p.value, method='fdr')
  original_corr <- corr
  
  # filter genes with FDR > 0.1
  corr <- corr[rownames(corr)!=gene_name,]
  corr <- corr[corr$FDR < 0.1,]
  print(paste0('totol number of genes with FDR < 0.1: ', nrow(corr)))
  
  # select and sort positive Rho genes
  corr_posi <- corr[corr$rho >0,]
  corr_posi <- corr_posi[order(corr_posi$p.value, decreasing = FALSE),]
  print(paste0('positive genes number: ', nrow(corr_posi)))
  
  # select and sort negative Rho genes
  corr_nega <- corr[corr$rho <0,]
  corr_nega <- corr_nega[order(corr_nega$p.value, decreasing = FALSE),]
  print(paste0('negative genes number: ', nrow(corr_nega)))
  
  # check control genes' correlations
  ctrl_genes <- ctrl_genes[ctrl_genes != gene_name]
  filtered_ctrl_genes <- ctrl_genes[!ctrl_genes %in% filtered_genes]
  print(paste0('Express < 20% cells control genes : ',paste(filtered_ctrl_genes, collapse = ', ')))
  top_ctrl_genes <- ctrl_genes[ctrl_genes %in% c(rownames(corr_posi[1:top_n,]),rownames(corr_nega[1:top_n,]))]
  print(paste0('Top 30 correlated control genes: ',paste(top_ctrl_genes, collapse = ', ')))
  ctrl_genes <- ctrl_genes[!ctrl_genes %in% c(rownames(corr_posi[1:top_n,]),rownames(corr_nega[1:top_n,]))]
  
  # sort positive control genes
  posi_ctrl <- ctrl_genes[ctrl_genes %in% rownames(corr_posi)]
  posi_order <- order(corr_posi[rownames(corr_posi) %in% posi_ctrl,]$p.value, decreasing = FALSE)
  posi_ctrl <- posi_ctrl[posi_order]
  print(paste0('Positively correlated control genes with FDR < 0.1: ',paste(posi_ctrl, collapse = ', ')))
  
  # sort negative control genes
  nega_ctrl <- ctrl_genes[ctrl_genes %in% rownames(corr_nega)]
  nega_order <- order(corr_nega[rownames(corr_nega) %in% nega_ctrl,]$p.value, decreasing = FALSE)
  nega_ctrl <- nega_ctrl[nega_order]
  print(paste0('Negatively correlated control genes with FDR < 0.1: ',paste(nega_ctrl, collapse = ', ')))
  
  # check FDR < 0.1 number (corr_df is genes with FDR < 0.1)
  if(nrow(corr) == 0){
    return(list(original_df=original_corr, corr_df=corr, heatmap=NULL))
  }
  
  # sum up and plot gene logcounts
  top_corr_genes <- c(rownames(na.omit(corr_posi[1:top_n,])),
                      gene_name,posi_ctrl,nega_ctrl,
                      rownames(na.omit(corr_nega[1:top_n,])))
  gap1 <- nrow(na.omit(corr_posi[1:top_n,]))
  gap2 <- gap1+length(c(gene_name,posi_ctrl))
  gap3 <- gap2+length(nega_ctrl)
  ifelse(gap3==gap2, gap <- c(gap1,gap2), gap <- c(gap1,gap2,gap3))
  top_corr_plot <- plotHeatmap(sce, exprs_values = 'logcounts', 
                               features = top_corr_genes,
                               columns = names(sort(assays(sce)$logcounts[gene_name,])),
                               cluster_rows = FALSE, cluster_cols=FALSE, 
                               center = TRUE, zlim = c(-5,5),
                               color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40),
                               main = paste0('Top ',top_n,' positively & negatively correlated genes in ', clust_name),
                               gaps_row = gap,
                               color_columns_by=c("label", "Cell.Type"))
  
  # change the color of target gene
  fontsize = top_corr_plot$gtable$grobs[[3]]$gp$fontsize
  names = top_corr_plot$gtable$grobs[[3]]$label
  color = rep('black', length(top_corr_genes))
  color[grep(paste0('^',gene_name,'$'),names)] <- 'red'
  for(n in c(top_ctrl_genes,posi_ctrl,nega_ctrl)){
    color[grep(paste0('^',n,'$'),names)] <- 'blue'
  }
  top_corr_plot$gtable$grobs[[3]]$gp <- grid::gpar(col=color, fontsize=fontsize)
  
  # (corr_df is genes with FDR < 0.1)
  return(list(original_df=original_corr, corr_df=corr, heatmap=top_corr_plot))
}

spear_corr_heatmap <- function(sce,clust,target,ctrl,
                               title_suffix='',top_num,file_prefix){
  # 'clust' is a list
  results <- list()
  for(i in clust){
    name <- ifelse(length(i)==1, as.character(i), paste(i,collapse = '&'))
    print(name)
    title_name <- paste0('cluster ',name)
    file_name <- paste0(file_prefix,'clust',name,'.jpg')
    spear_corr$out <- capture.output(spear_corr <- spearman_corr(sce,i,target,
                                                                 ctrl,title_name,
                                                                 top_num))
    print(spear_corr$out)
    ggsave(file_name,spear_corr$heatmap,device='jpg', width = 12, height = 12)
    results[[name]] <- spear_corr
  }
  return(results)
}

positive_control_genes <- c("Pou5f1",'Sox2','Nanog','Zfp42','Klf4','Klf2','Klf5',
                            'Esrrb','Eras','Nacc1','Utf1','Lefty1', 'Sall4',
                            'Smad3','Smad1','Gdf3', 'Lin28a','Tfcp2l1','Id3',
                            'Fgf4','Tcl1','Spp1','Upp1','Fbxo15','Dppa3','Trp53',
                            'Dppa5a','Nr0b1','Stat3','Cdyl', 'Mycbp', 'Tbx3', 
                            'Zfx', "Prdm14",'Foxd3', 'Gbx2','Zfp143','Otx2',
                            'Cdx2', "Gata3", "Eomes", 'Tcf15',"Tcf3",'Dnmt3a','Dnmt3b')


# oct4 spearman correlation ####################################################
oct4_spear_corr <- spear_corr_heatmap(sce_dev,list(1,2,3,4,5,6,7,8,9,10,
                                                   c(3,7),c(4,6),c(5,8,10),
                                                   c(2,5,8,10)),
                                      'Pou5f1',
                                      positive_control_genes,'',
                                      30,'figures/spear_corr_')
oct4_spear_corr[['1']]$out
oct4_spear_corr[['1']]$heatmap

# Nanog spearman correlation ###################################################
nanog_spear_corr <- spear_corr_heatmap(sce_dev,list(1,2,3,4,5,6,7,8,9,10,
                                                    c(3,7),c(4,6),c(5,8,10),
                                                    c(2,5,8,10)),
                                       'Nanog',
                                       positive_control_genes,'',
                                       30,'figures/spear_corr_nanog_')
nanog_spear_corr[['1']]$out
nanog_spear_corr[['1']]$heatmap

# Sox2 spearman correlation ###################################################
sox2_spear_corr <- spear_corr_heatmap(sce_dev,list(1,2,3,4,5,6,7,8,9,10,
                                                   c(3,7),c(4,6),c(5,8,10),
                                                   c(2,5,8,10)),
                                      'Sox2',
                                      positive_control_genes,'',
                                      30,'figures/spear_corr_sox2_')
sox2_spear_corr[['1']]$out
sox2_spear_corr[['1']]$heatmap

# venn diagram comparison ######################################################
# check in vivo intersection
find_inter <- function(li){
  r <- list()
  for(i in 2:length(li)){
    arr <- combn(1:length(li),i)
    names <- c()
    for(j in 1:ncol(arr)){
      n <- Reduce(intersect, li[arr[,j]])
      names <- c(names, n)
    }
    r[[i]] <- unique(names)
  }
  return(r)
}

# # in vivo intersections (oct4)
# ggVennDiagram(list(E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0,]), 
#                    E4_nega = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho < 0,]),
#                    E5.5_posi = rownames(oct4_spear_corr[['1']]$corr_df[oct4_spear_corr[['1']]$corr_df$rho > 0,]), 
#                    E5.5_nega = rownames(oct4_spear_corr[['1']]$corr_df[oct4_spear_corr[['1']]$corr_df$rho < 0,])),
#               label = "count", label_size = 5) +
#   theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
#   ggtitle('intersections of correlated genes in E3.5&E4.5 and E5.5') +
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave("figures/oct4_invivo_venn.pdf",device='pdf', width = 12, height = 12)


# # positive to negative (in vivo)
# posi2nega_vivo <- find_inter(list(cluster6_posi = rownames(oct4_spear_corr[['6']]$corr_df[oct4_spear_corr[['6']]$corr_df$rho > 0,]),
#                                   cluster2_nega = rownames(oct4_spear_corr[['2']]$corr_df[oct4_spear_corr[['2']]$corr_df$rho < 0,])))[[2]]
# posi2nega_vivo
# posi2nega_vivo[posi2nega_vivo %in% positive_control_genes]
# 
# # positive to positive (in vivo)
# posi2posi_vivo <- find_inter(list(cluster6_posi = rownames(oct4_spear_corr[['6']]$corr_df[oct4_spear_corr[['6']]$corr_df$rho > 0,]),
#                                   cluster2_posi = rownames(oct4_spear_corr[['2']]$corr_df[oct4_spear_corr[['2']]$corr_df$rho > 0,])))[[2]]
# posi2posi_vivo
# posi2posi_vivo[posi2posi_vivo %in% positive_control_genes]
# 
# # negative to positive (in vivo)
# nega2posi_vivo <- find_inter(list(cluster6_nega = rownames(oct4_spear_corr[['6']]$corr_df[oct4_spear_corr[['6']]$corr_df$rho < 0,]),
#                                   cluster2_posi = rownames(oct4_spear_corr[['2']]$corr_df[oct4_spear_corr[['2']]$corr_df$rho > 0,])))[[2]]
# nega2posi_vivo
# nega2posi_vivo[nega2posi_vivo %in% positive_control_genes]
# 
# # negative to negative (in vivo)
# nega2nega_vivo <- find_inter(list(cluster6_nega = rownames(oct4_spear_corr[['6']]$corr_df[oct4_spear_corr[['6']]$corr_df$rho < 0,]),
#                                   cluster2_posi = rownames(oct4_spear_corr[['2']]$corr_df[oct4_spear_corr[['2']]$corr_df$rho < 0,])))[[2]]
# nega2nega_vivo
# nega2nega_vivo[nega2nega_vivo %in% positive_control_genes]
# 
# # in vitro intersections
# ggVennDiagram(list(Serum_posi = rownames(oct4_spear_corr[['4']]$corr_df[oct4_spear_corr[['4']]$corr_df$rho > 0,]), 
#                    Serum_nega = rownames(oct4_spear_corr[['4']]$corr_df[oct4_spear_corr[['4']]$corr_df$rho < 0,]),
#                    Epi_posi = rownames(oct4_spear_corr[['5&7']]$corr_df[oct4_spear_corr[['5&7']]$corr_df$rho > 0,]), 
#                    Epi_nega = rownames(oct4_spear_corr[['5&7']]$corr_df[oct4_spear_corr[['5&7']]$corr_df$rho < 0,])),
#               label = "count", label_size = 5) +
#   theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
#   ggtitle('intersections of correlated genes in Serum and EpiStemCell') +
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave("figures/oct4_invitro_venn.pdf",device='pdf', width = 10, height = 10)
# 
# # check in vitro intersection
# posi2nega_vitro <- find_inter(list(cluster6_posi = rownames(oct4_spear_corr[['4']]$corr_df[oct4_spear_corr[['4']]$corr_df$rho > 0,]),
#                                    cluster2_nega = rownames(oct4_spear_corr[['5&7']]$corr_df[oct4_spear_corr[['5&7']]$corr_df$rho < 0,])))[[2]]
# posi2nega_vitro
# posi2nega_vitro[posi2nega_vitro %in% positive_control_genes]
# nega2posi_vitro <- find_inter(list(cluster6_nega = rownames(oct4_spear_corr[['4']]$corr_df[oct4_spear_corr[['4']]$corr_df$rho < 0,]),
#                                    cluster2_posi = rownames(oct4_spear_corr[['5&7']]$corr_df[oct4_spear_corr[['5&7']]$corr_df$rho > 0,])))[[2]]
# nega2posi_vitro
# nega2posi_vitro[nega2posi_vitro %in% positive_control_genes]
# 
# # check intersections of in vivo intersections and in vitro intersection
# find_inter(list(vitro=posi2nega_vitro,vivo=posi2nega_vivo))[[2]]
# find_inter(list(vitro=nega2posi_vitro,vivo=nega2posi_vivo))[[2]]

# check in vivo and in vitro intersections #####################################
# oct4
ggVennDiagram(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                               oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]), 
                   Serum_nega = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho < 0 & 
                                                                               oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                   E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                       oct4_spear_corr[['9']]$corr_df$FDR <0.1,]), 
                   E4_nega = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho < 0 & 
                                                                       oct4_spear_corr[['9']]$corr_df$FDR <0.1,]),
                   predicted = oct4_tar.potent),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of correlated genes in Serum and E4 (Oct4)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/oct4_serum_E4_venn.pdf",device='pdf', width = 10, height = 10)

# oct4_posi_serum_E4 <- find_inter(list(Serum_posi = rownames(oct4_spear_corr[['4']]$corr_df[oct4_spear_corr[['4']]$corr_df$rho > 0,]),
#                 E4_posi = rownames(oct4_spear_corr[['6']]$corr_df[oct4_spear_corr[['6']]$corr_df$rho > 0,])))[[2]]
# oct4_posi_serum_E4 <- oct4_spear_corr[['4']]$corr_df[rownames(oct4_spear_corr[['4']]$corr_df) %in% oct4_posi_serum_E4,]
# oct4_posi_serum_E4 <- oct4_posi_serum_E4[order(oct4_posi_serum_E4$FDR, decreasing = FALSE),]
# find_inter(list(rownames(oct4_posi_serum_E4),oct4_tar.potent))[[2]]
# 
# oct4_nega_serum_E4 <- find_inter(list(Serum_nega = rownames(oct4_spear_corr[['4']]$corr_df[oct4_spear_corr[['4']]$corr_df$rho < 0,]),
#                                       E4_nega = rownames(oct4_spear_corr[['6']]$corr_df[oct4_spear_corr[['6']]$corr_df$rho < 0,])))[[2]]
# oct4_nega_serum_E4 <- oct4_spear_corr[['4']]$corr_df[rownames(oct4_spear_corr[['4']]$corr_df) %in% oct4_nega_serum_E4,]
# oct4_nega_serum_E4 <- oct4_nega_serum_E4[order(oct4_nega_serum_E4$FDR, decreasing = FALSE),]
# find_inter(list(rownames(oct4_nega_serum_E4),oct4_tar.potent))[[2]]
# 
find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                            oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                    oct4_spear_corr[['9']]$corr_df$FDR <0.1,]),
                predicted = oct4_tar.potent))[[3]]

# nanog 
ggVennDiagram(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]), 
                   Serum_nega = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho < 0 & 
                                                                                nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                   E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                        nanog_spear_corr[['9']]$corr_df$FDR <0.1,]), 
                   E4_nega = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho < 0 & 
                                                                        nanog_spear_corr[['9']]$corr_df$FDR <0.1,]),
                   predicted = oct4_tar.potent),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of correlated genes in Serum and E4 (Nanog)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/nanog_serum_E4_venn.pdf",device='pdf', width = 10, height = 10)

find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                             nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                     nanog_spear_corr[['9']]$corr_df$FDR <0.01,]),
                predicted = oct4_tar.potent))[[3]]

# sox2
ggVennDiagram(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                               sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]), 
                   Serum_nega = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho < 0 & 
                                                                               sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                   E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                       sox2_spear_corr[['9']]$corr_df$FDR <0.1,]), 
                   E4_nega = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho < 0 & 
                                                                       sox2_spear_corr[['9']]$corr_df$FDR <0.1,]),
                   predicted = oct4_tar.potent),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of correlated genes in Serum and E4 (Sox2)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/sox2_serum_E4_venn.pdf",device='pdf', width = 10, height = 10)

find_inter(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                            sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                    sox2_spear_corr[['9']]$corr_df$FDR <0.01,]),
                predicted = oct4_tar.potent))[[3]]

# oct4 & sox2 & nanog
ggVennDiagram(list(Oct4_targets = find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                                              oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                                                  E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                                      oct4_spear_corr[['9']]$corr_df$FDR <0.1,])))[[2]], 
                   Nanog_targets = find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                                                nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                                                   E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                                        nanog_spear_corr[['9']]$corr_df$FDR <0.1,])))[[2]],
                   Sox2_targets = find_inter(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                                              sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                                                  E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                                      sox2_spear_corr[['9']]$corr_df$FDR <0.1,])))[[2]], 
                   predicted = oct4_tar.potent),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of correlated genes in Serum and E4 (Oct4, Sox2, Nanog)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/all_serum_E4_venn.pdf",device='pdf', width = 10, height = 10)

# # oct4 sox2
# find_inter(list(find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                                 E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0,])))[[2]],
#                 find_inter(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                                 E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0,])))[[2]]))[[2]]
# 
# # oct4 nanog
# find_inter(list(find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                                 E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0,])))[[2]],
#                 find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                                 E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0,])))[[2]]))[[2]]
# 
# # sox2 nanog
# find_inter(list(find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                                 E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0,])))[[2]],
#                 find_inter(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                                 E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0,])))[[2]]))[[2]]

# oct4 sox2 nanog
find_inter(list(find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                            oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                    oct4_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                             nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                     nanog_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                find_inter(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                            sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                    sox2_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]]))[[3]]

# oct4 sox2 nanog candidates
find_inter(list(find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                            oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                    oct4_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                             nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                     nanog_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                find_inter(list(Serum_posi = rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                            sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                    sox2_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                predicted = oct4_tar.potent))[[4]]

# oct4 nanog candidates
find_inter(list(find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                            oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                    oct4_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                find_inter(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                             nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.01,]),
                                E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                     nanog_spear_corr[['9']]$corr_df$FDR <0.01,])))[[2]],
                predicted = oct4_tar.potent))[[3]]

# oct4 depletion bulk RNA-seq targets
oct4_RNA_seq_targets <- read.csv('oct4 depletion RNA seq.csv',header=FALSE)
oct4_RNA_seq_targets <- oct4_RNA_seq_targets[,1]
head(oct4_RNA_seq_targets)
ggVennDiagram(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                               oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]), 
                   E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                       oct4_spear_corr[['9']]$corr_df$FDR <0.1,]), 
                   RNAseq_targets = oct4_RNA_seq_targets,
                   predicted = oct4_tar.potent),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of correlated genes in Serum, E4 and RNA-seq (Oct4)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/oct4_serum_E4_RNAseq_venn.pdf",device='pdf', width = 10, height = 10)

# nanog depletion bulk RNA-seq targets
ggVennDiagram(list(Serum_posi = rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]), 
                   E4_posi = rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                        nanog_spear_corr[['9']]$corr_df$FDR <0.1,]), 
                   RNAseq_targets = nanog_RNA_seq_targets,
                   predicted = oct4_tar.potent),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of correlated genes in Serum, E4 and RNAseq + ChIPseq (Nanog)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/nanog_serum_E4_RNAseq_venn.pdf",device='pdf', width = 10, height = 10)

# find genes with oct4/nanog/sox2 in their top30 correlated genes list###################
rank_by_pvalue <- function(sce_obj,clust_num, examiners,gene_name, target){
  # calculate the fdr rank of target in all genes correlated to 'gene_name'
  sce <- sce_obj[,sce_obj$label %in% clust_num]
  corr <- mclapply(examiners, FUN = get_corr, name2=gene_name, sce=sce, mc.cores=16)
  corr <- do.call(rbind, corr)
  corr <- as.data.frame(corr)
  corr$FDR <- p.adjust(corr$p.value, method='fdr')
  corr <- corr[rownames(corr)!=gene_name,]
  corr <- corr[rowSums(is.na(corr)) == 0, ]
  corr <- corr[corr$FDR < 0.1,]
  corr_posi <- corr[corr$rho >0,]
  if(! target %in% rownames(corr_posi)){
    result <- t(as.matrix(c(NA,NA,NA,NA,NA,NA,NA,NA)))
    colnames(result) <- c('rho','p.value','FDR','rank','posi_num','rho_diff','pvalue_diff','expr_ratio')
    rownames(result) <- gene_name
    return(result)
  }
  corr_posi <- corr_posi[order(corr_posi$p.value, decreasing = FALSE),]
  corr_posi$rank <- rank(corr_posi$p.value)
  corr_posi$posi_num <- nrow(corr_posi)
  result <- corr_posi[rownames(corr_posi) == target,]
  if(result$rank == result$posi_num){
    rho_diff <- 0
    pvalue_diff <- 0
  }
  else{
    next_rank_index <- match(target,rownames(corr_posi))+1
    next_rho <- corr_posi[next_rank_index, 'rho']
    rho_diff <- result$rho - next_rho
    next_pvalue <- corr_posi[next_rank_index, 'p.value']
    pvalue_diff <- next_pvalue - result$p.value
  }
  result$rho_diff <- rho_diff
  result$pvalue_diff <- pvalue_diff
  result$expr_ratio <- nexprs(sce[rownames(sce)==gene_name,], byrow=TRUE)/ncol(sce)
  rownames(result) <- gene_name
  return(result)
}

filter_genes <- function(sce_obj,clust_num,ratio){
  sce <- sce_obj[,sce_obj$label %in% clust_num]
  cell_num <- round(ncol(sce)*ratio,0)
  filtered_genes <- rownames(sce)[nexprs(sce, byrow=TRUE) > cell_num]
  return(filtered_genes)
}

rank_by_pvalue(sce_dev,9,oct4_tar.potent,'Sf3b4','Pou5f1')
# candidates <- find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0,]),
#                               E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0,])))[[2]]

# Oct4 E4
candidates <- rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                        oct4_spear_corr[['9']]$corr_df$FDR <0.01,])
print(system.time(oct4_rank_E4 <- mclapply(candidates, 
                                        FUN = rank_by_pvalue, 
                                        sce_obj=sce_dev, 
                                        clust_num=c(9),
                                        examiners = filter_genes(sce_dev,
                                                                 clust_num = c(9),
                                                                 ratio = 0.5), 
                                        target = 'Pou5f1',
                                        mc.cores=3)))
oct4_rank_E4 <- do.call(rbind,oct4_rank_E4)
head(oct4_rank_E4)
write.csv(oct4_rank_E4, 'figures/oct4_rank_E4.csv')
rm(candidates)

oct4_rank_E4 %>% 
  arrange(rank) %>% 
  dplyr::filter(rank <= 100 & expr_ratio > 0.5 & rho > 0.3) %>% 
  rownames(.)

# Oct4 serum
candidates <- rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                             oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.01,])
print(system.time(oct4_rank_serum <- mclapply(candidates, 
                                           FUN = rank_by_pvalue, 
                                           sce_obj=sce_dev, 
                                           clust_num=c(5,8,10),
                                           examiners = filter_genes(sce_dev,
                                                                    clust_num = c(5,8,10),
                                                                    ratio = 0.5), 
                                           target = 'Pou5f1',
                                           mc.cores=3)))
oct4_rank_serum <- do.call(rbind,oct4_rank_serum)
head(oct4_rank_serum)
write.csv(oct4_rank_serum, 'figures/oct4_rank_serum.csv')
rm(candidates)

oct4_rank_serum %>% 
  arrange(rank) %>% 
  dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
  rownames(.)

# Nanog E4
candidates <- rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                         nanog_spear_corr[['9']]$corr_df$FDR <0.01,])
print(system.time(nanog_rank_E4 <- mclapply(candidates, 
                                           FUN = rank_by_pvalue, 
                                           sce_obj=sce_dev, 
                                           clust_num=c(9),
                                           examiners = filter_genes(sce_dev,
                                                                    clust_num = c(9),
                                                                    ratio = 0.5), 
                                           target = 'Pou5f1',
                                           mc.cores=3)))
nanog_rank_E4 <- do.call(rbind,nanog_rank_E4)
head(nanog_rank_E4)
write.csv(nanog_rank_E4, 'figures/nanog_rank_E4.csv')
rm(candidates)

# Nanog serum
candidates <- rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                              nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.01,])
print(system.time(nanog_rank_serum <- mclapply(candidates, 
                                              FUN = rank_by_pvalue, 
                                              sce_obj=sce_dev, 
                                              clust_num=c(5,8,10),
                                              examiners = filter_genes(sce_dev,
                                                                       clust_num = c(5,8,10),
                                                                       ratio = 0.5), 
                                              target = 'Pou5f1',
                                              mc.cores=3)))
nanog_rank_serum <- do.call(rbind,nanog_rank_serum)
head(nanog_rank_serum)
write.csv(nanog_rank_serum, 'figures/nanog_rank_serum.csv')

# find overlaps
ggVennDiagram(list(Oct4_E4 = oct4_rank_E4 %>% 
                               arrange(rank) %>% 
                               dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
                               rownames(.), 
                   Oct4_serum = oct4_rank_serum %>% 
                                  arrange(rank) %>% 
                                  dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
                                  rownames(.),
                   Nanog_E4 = nanog_rank_E4 %>% 
                                arrange(rank) %>% 
                                dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
                                rownames(.),
                   Nanog_serum = nanog_rank_serum %>% 
                                   arrange(rank) %>% 
                                   dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
                                   rownames(.)),
              label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersections of tightly correlated genes') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2))
ggsave("figures/tightly_correlated_genes_venn.pdf",device='pdf', width = 10, height = 10)

find_inter(list(oct4_rank_E4 %>% 
                  arrange(rank) %>% 
                  dplyr::filter(rank <= 10 & expr_ratio > 0.5 & rho > 0.3) %>% 
                  rownames(.),
                nanog_rank_E4 %>% 
                  arrange(rank) %>% 
                  dplyr::filter(rank <= 10 & expr_ratio > 0.5 & rho > 0.3) %>% 
                  rownames(.)
))[[2]]

# # Sox2 E4
# candidates <- rownames(sox2_spear_corr[['9']]$corr_df[sox2_spear_corr[['9']]$corr_df$rho > 0 & 
#                                                         sox2_spear_corr[['9']]$corr_df$FDR <0.01,])
# print(system.time(sox2_rank_E4 <- mclapply(candidates, 
#                                             FUN = rank_by_pvalue, 
#                                             sce_obj=sce_dev, 
#                                             clust_num=c(9),
#                                             examiners = filter_genes(sce_dev,
#                                                                      clust_num = c(9),
#                                                                      ratio = 0.5), 
#                                             target = 'Pou5f1',
#                                             mc.cores=3)))
# sox2_rank_E4 <- do.call(rbind,sox2_rank_E4)
# head(sox2_rank_E4)
# write.csv(sox2_rank_E4, 'figures/sox2_rank_E4.csv')
# rm(candidates)
# 
# # Sox2 serum
# candidates <- rownames(sox2_spear_corr[['5&8&10']]$corr_df[sox2_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
#                                                              sox2_spear_corr[['5&8&10']]$corr_df$FDR <0.01,])
# print(system.time(sox2_rank_serum <- mclapply(candidates, 
#                                                FUN = rank_by_pvalue, 
#                                                sce_obj=sce_dev, 
#                                                clust_num=c(5,8,10),
#                                                examiners = filter_genes(sce_dev,
#                                                                         clust_num = c(5,8,10),
#                                                                         ratio = 0.5), 
#                                                target = 'Pou5f1',
#                                                mc.cores=3)))
# sox2_rank_serum <- do.call(rbind,sox2_rank_serum)
# head(sox2_rank_serum)
# write.csv(sox2_rank_serum, 'figures/sox2_rank_serum.csv')
# rm(candidates)

# check correlated genes number (FDR < 0.1)##################################################
fdr_lt_0.1_num <- function(sce_obj,clust_num, examiners,gene_name){
  sce <- sce_obj[,sce_obj$label %in% clust_num]
  corr <- mclapply(examiners, FUN = get_corr, name2=gene_name, sce=sce, mc.cores=16)
  corr <- do.call(rbind, corr)
  corr <- as.data.frame(corr)
  corr$FDR <- p.adjust(corr$p.value, method='fdr')
  corr <- corr[rownames(corr)!=gene_name,]
  corr <- corr[rowSums(is.na(corr)) == 0, ]
  corr <- corr[corr$FDR < 0.1,]
  corr_posi <- corr[corr$rho >0,]
  corr_nega <- corr[corr$rho <0,]
  count <- t(as.matrix(c(nrow(corr_posi),nrow(corr_nega))))
  colnames(count) <- c('posi_counts','nega_counts')
  rownames(count) <- gene_name
  return(count)
}
# oct4 potential target genes (E4)
print(system.time(fdr_count_oct4_E4 <- mclapply(unique(c(positive_control_genes,oct4_tar.potent)), 
                                           FUN = fdr_lt_0.1_num, 
                                           sce_obj=sce_dev, 
                                           clust_num=c(9),
                                           examiners = filter_genes(sce_dev,
                                                                    clust_num = c(9),
                                                                    ratio = 0.2), 
                                           mc.cores=3)))
fdr_count_oct4_E4 <- do.call(rbind,fdr_count_oct4_E4)
fdr_count_oct4_E4 <- as.data.frame(fdr_count_oct4_E4)
fdr_count_oct4_E4$total <- fdr_count_oct4_E4$posi_counts + fdr_count_oct4_E4$nega_counts
head(fdr_count_oct4_E4)

# oct4 potential target genes (serum)
print(system.time(fdr_count_oct4_serum <- mclapply(unique(c(positive_control_genes,oct4_tar.potent)), 
                                                FUN = fdr_lt_0.1_num, 
                                                sce_obj=sce_dev, 
                                                clust_num=c(5,8,10),
                                                examiners = filter_genes(sce_dev,
                                                                         clust_num = c(5,8,10),
                                                                         ratio = 0.2), 
                                                mc.cores=3)))
fdr_count_oct4_serum <- do.call(rbind,fdr_count_oct4_serum)
fdr_count_oct4_serum <- as.data.frame(fdr_count_oct4_serum)
fdr_count_oct4_serum$total <- fdr_count_oct4_serum$posi_counts+fdr_count_oct4_serum$nega_counts
head(fdr_count_oct4_serum)

# check overlaps
find_inter(list(E4 = fdr_count_oct4_E4 %>% dplyr::slice_max(n = 30, order_by = total) %>% rownames(),
                serum = fdr_count_oct4_serum %>% dplyr::slice_max(n = 30, order_by = total) %>% rownames()))[[2]]


# enrichment analysis ##########################################################
enrich_analysis <- function(candidates,ensdb,orgdb,GO_term,row_number,organism,file_name){
  # candidates are gene names
  # GO term is a list containing: 'CC'/'MF'/'BP'
  # row_number is the row number for each GO term subplot and KEGG plot
  target_df <- AnnotationDbi::select(ensdb, keys=candidates, 
                                     keytype='GENENAME', column='ENTREZID')
  dupli_num <- sum(duplicated(target_df$GENENAME))  # check duplication number
  print(paste0('duplicates number: ',dupli_num))  
  target_df <- target_df[!duplicated(target_df$GENENAME),]
  
  # GO
  csv_filename <- paste0('figures/enrich_go_',file_name,'.csv')
  jpg_suffix <- paste(GO_term,collapse='_')
  jpg_filename <- paste0('figures/enrich_go_',file_name,'_',jpg_suffix,'.jpg')
  enrich_go = enrichGO(gene = target_df$ENTREZID, 
                               OrgDb = orgdb, 
                               keyType = "ENTREZID", 
                               ont = "ALL", # can be selected in "CC\MF\BP\ALL"
                               pAdjustMethod = "fdr", 
                               pvalueCutoff = 0.05,  
                               qvalueCutoff = 0.1, # default is 0.2
                               readable = T) 
  enrich_go <- data.frame(enrich_go)
  write.csv(enrich_go,csv_filename,row.names = TRUE)
  enrich_go$ratio <- sapply(enrich_go$GeneRatio, FUN = function(x){eval(parse(text=x))})
  ggplot(enrich_go %>% 
           dplyr::filter(ONTOLOGY %in% GO_term) %>%
           group_by(ONTOLOGY) %>% 
           dplyr::slice_min(n = row_number, order_by = p.adjust),
         aes(y=reorder(Description,-p.adjust),x=ratio))+
    geom_point(aes(size=Count,color=p.adjust))+
    facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
    scale_color_gradient(low = "red",high ="blue")+
    labs(color=expression(p.adjust,size="Count"), 
         x="Gene Ratio",y="GO term",title="GO Enrichment (Biological Process)")+
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12))
  ggsave(jpg_filename,device='jpg', width = 14, height = 8)
  
  # KEGG
  csv_filename <- paste0('figures/enrich_kegg_',file_name,'.csv')
  jpg_filename <- paste0('figures/enrich_kegg_',file_name,'.jpg')
  enrich_kegg = enrichKEGG(gene = target_df$ENTREZID, 
                                   keyType = "kegg", 
                                   organism = organism,
                                   pAdjustMethod = "fdr", 
                                   pvalueCutoff = 0.05,  
                                   qvalueCutoff = 0.1 # default is 0.2
  ) 
  enrich_kegg <- data.frame(enrich_kegg)
  write.csv(enrich_kegg,csv_filename,row.names = TRUE)
  enrich_kegg$Description <- gsub('- Mus musculus \\(house mouse\\)','',enrich_kegg$Description)
  enrich_kegg$ratio <- sapply(enrich_kegg$GeneRatio, FUN = function(x){eval(parse(text=x))})
  ggplot(enrich_kegg %>% 
           dplyr::slice_min(n = 30, order_by = p.adjust),
         aes(y=reorder(Description,-p.adjust),x=ratio))+
    geom_point(aes(size=Count,color=p.adjust))+
    scale_color_gradient(low = "red",high ="blue")+
    labs(color=expression(p.adjust,size="Count"), 
         x="Gene Ratio",y="KEGG term",title="KEGG Enrichment")+
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12))
  ggsave(jpg_filename,device='jpg', width = 10, height = 8)
}

# get orgdb
orgdb.mm <- AnnotationHub()[["AH107060"]]

# oct4 RNA seq
enrich_analysis(oct4_RNA_seq_targets,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_RNAseq')

# oct4 direct targets (ChIP-seq)
enrich_analysis(oct4_tar.potent,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_ChIPseq')

# oct4 E4 
candidates <- rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                        oct4_spear_corr[['9']]$corr_df$FDR <0.1,])
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_E4')

# oct4 serum 
candidates <- rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                        oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,])
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_serum')

# oct4 serum&E4
candidates <- find_inter(list(Serum_posi = rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                                                          oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
                              E4_posi = rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                                                  oct4_spear_corr[['9']]$corr_df$FDR <0.1,])))[[2]]
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_serum&E4')

# oct4 tightly correlated genes in E4
candidates <- oct4_rank_E4 %>% 
  arrange(rank) %>% 
  dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
  rownames(.)
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_E4_tight')

# oct4 tightly correlated genes in serum
candidates <- oct4_rank_serum %>% 
  arrange(rank) %>% 
  dplyr::filter(rank <= 30 & expr_ratio > 0.5 & rho > 0.3) %>% 
  rownames(.)
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','oct4_serum_tight')

# nanog E4 
candidates <- rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                         nanog_spear_corr[['9']]$corr_df$FDR <0.1,])
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','nanog_E4')

# nanog serum 
candidates <- rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                             nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.1,])
enrich_analysis(candidates,ens.mm.v109,orgdb.mm,c('BP'),30,'mmu','nanog_serum')

# clean up
rm(candidates)


# # ChIP-seq peak annotation #####################################################
# peak_anno_nanog <- annotatePeak("./ChIP_nanog_mm10.bed", 
#                                 TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
#                                 tssRegion = c(-1000, 1000),
#                                 annoDb="org.Mm.eg.db",
#                                 addFlankGeneInfo = TRUE,
#                                 flankDistance=5000)
# 
# # nanog ChIPseq targets
# as.data.frame(peak_anno_nanog) %>% head()
# # candidates <- as.data.frame(peak_anno_nanog) %>% 
# #   dplyr::filter(distanceToTSS < 1000 & distanceToTSS > -1000) %>% .$geneId
# candidates <- as.data.frame(peak_anno_nanog) %>% .$geneId
# length(unique(candidates))
# chip_targets_nanog <- AnnotationDbi::select(ens.mm.v109, 
#                                             keys=candidates, 
#                                             keytype='ENTREZID', column='GENENAME') %>% .$GENENAME
# head(chip_targets_nanog)
# length(chip_targets_nanog)
# 
# # # check genes within the flankDistance 
# # id <- as.data.frame(peak_anno_nanog)$flank_geneIds
# # sum(is.na(id))
# # id <- id[!is.na(id)]
# # length(id)
# # sum(stringr::str_detect(id,';'))  # check how many peaks contain multiple genes in their flank distance
# # id <- id[!stringr::str_detect(id,';')]
# # length(id)
# # rm(id)
# 
# # nanog RNAseq targets
# nanog_RNA_seq_targets <- read.csv('RNAseq_nanog.csv',header=FALSE)
# nanog_RNA_seq_targets <- nanog_RNA_seq_targets[,2]
# head(nanog_RNA_seq_targets)
# 
# # clean up
# rm(candidates)

# other enrichment analysis tools (web-based method) ###########################
# recommended web-based tools: ChEA3, Enrichr-KG
# we need to record the gene number successfully recognized by web-tools (e.g. 1309 below)

# oct4 serum ------------------------------------------------------------
write.table(rownames(oct4_spear_corr[['5&8&10']]$corr_df[oct4_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                           oct4_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
            'figures/targets_oct4_serum.txt',
            row.names = FALSE,col.names = FALSE, quote=FALSE)
# literature enrichment (ChIP-X by ChEA3)
liter_oct4_serum_df <- read.csv('figures/Literature_ChIP-seq_oct4_serum_1309.tsv',sep='\t')
liter_oct4_serum_df$ratio <- liter_oct4_serum_df$Intersect / 1309  
ggplot(liter_oct4_serum_df %>% 
         dplyr::slice_min(n = 30, order_by = FDR),
       aes(y=reorder(TF,-FDR),x=ratio))+
  geom_point(aes(size=Intersect,color=FDR))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(FDR,size="Intersect"), 
       x="Gene Ratio",y="Transcription Factor",title="literature ChIP-X Enrichment (Oct4 serum)")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  guides(size=guide_legend("Count"))
ggsave('figures/liter_oct4_serum.jpg',device='jpg', width = 8, height = 8)
# MGI enrichment (by Enrichr-KG)
mgi_oct4_serum_df <- read.csv('figures/Enrichr-KG_MGI_oct4_serum_1316.csv')
ggplot(mgi_oct4_serum_df %>%
         dplyr::slice_min(n = 30, order_by = q.value),
       aes(x=reorder(Term,-q.value),y=-log10(q.value))) +
  geom_bar(stat="identity") +
  labs(x="MGI term",y="-log10(qvalue)",title="MGI Enrichment (Oct4 serum)") + 
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))
ggsave('figures/mgi_oct4_serum.jpg',device='jpg', width = 12, height = 8)

# oct4 E4 -----------------------------------------------
write.table(rownames(oct4_spear_corr[['9']]$corr_df[oct4_spear_corr[['9']]$corr_df$rho > 0 & 
                                                           oct4_spear_corr[['9']]$corr_df$FDR <0.1,]),
            'figures/targets_oct4_E4.txt',
            row.names = FALSE,col.names = FALSE, quote=FALSE)
# literature enrichment (ChIP-X by ChEA3)
liter_oct4_E4_df <- read.csv('figures/Literature_ChIP-seq_oct4_E4_1380.tsv',sep='\t')
liter_oct4_E4_df$ratio <- liter_oct4_E4_df$Intersect / 1380  
ggplot(liter_oct4_E4_df %>% 
         dplyr::slice_min(n = 30, order_by = FDR),
       aes(y=reorder(TF,-FDR),x=ratio))+
  geom_point(aes(size=Intersect,color=FDR))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(FDR,size="Intersect"), 
       x="Gene Ratio",y="Transcription Factor",title="literature ChIP-X Enrichment (Oct4 E4)")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  guides(size=guide_legend("Count"))
ggsave('figures/liter_oct4_E4.jpg',device='jpg', width = 8, height = 8)
# MGI enrichment (by Enrichr-KG)
mgi_oct4_E4_df <- read.csv('figures/Enrichr-KG_MGI_oct4_E4_1396.csv')
ggplot(mgi_oct4_E4_df %>%
         dplyr::slice_min(n = 30, order_by = q.value),
       aes(x=reorder(Term,-q.value),y=-log10(q.value))) +
  geom_bar(stat="identity") +
  labs(x="MGI term",y="-log10(qvalue)",title="MGI Enrichment (Oct4 E4)") + 
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))
ggsave('figures/mgi_oct4_E4.jpg',device='jpg', width = 12, height = 8)

# nanog serum --------------------------------------------
write.table(rownames(nanog_spear_corr[['5&8&10']]$corr_df[nanog_spear_corr[['5&8&10']]$corr_df$rho > 0 & 
                                                            nanog_spear_corr[['5&8&10']]$corr_df$FDR <0.1,]),
            'figures/targets_nanog_serum.txt',
            row.names = FALSE,col.names = FALSE, quote=FALSE)
# literature enrichment (ChIP-X by ChEA3)
liter_nanog_serum_df <- read.csv('figures/Literature_ChIP-seq_nanog_serum_1694.tsv',sep='\t')
liter_nanog_serum_df$ratio <- liter_nanog_serum_df$Intersect / 1694  
ggplot(liter_nanog_serum_df %>% 
         dplyr::slice_min(n = 30, order_by = FDR),
       aes(y=reorder(TF,-FDR),x=ratio))+
  geom_point(aes(size=Intersect,color=FDR))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(FDR,size="Intersect"), 
       x="Gene Ratio",y="Transcription Factor",title="literature ChIP-X Enrichment (Nanog serum)")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  guides(size=guide_legend("Count"))
ggsave('figures/liter_nanog_serum.jpg',device='jpg', width = 8, height = 8)
# MGI enrichment (by Enrichr-KG)
mgi_nanog_serum_df <- read.csv('figures/Enrichr-KG_MGI_nanog_serum_1696.csv')
ggplot(mgi_nanog_serum_df %>%
         dplyr::slice_min(n = 30, order_by = q.value),
       aes(x=reorder(Term,-q.value),y=-log10(q.value))) +
  geom_bar(stat="identity") +
  labs(x="MGI term",y="-log10(qvalue)",title="MGI Enrichment (Nanog serum)") + 
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))
ggsave('figures/mgi_nanog_serum.jpg',device='jpg', width = 12, height = 8)

# nanog E4 ---------------------------------------------
write.table(rownames(nanog_spear_corr[['9']]$corr_df[nanog_spear_corr[['9']]$corr_df$rho > 0 & 
                                                       nanog_spear_corr[['9']]$corr_df$FDR <0.1,]),
            'figures/targets_nanog_E4.txt',
            row.names = FALSE,col.names = FALSE, quote=FALSE)
# literature enrichment (ChIP-X by ChEA3)
liter_nanog_E4_df <- read.csv('figures/Literature_ChIP-seq_nanog_E4_1010.tsv',sep='\t')
liter_nanog_E4_df$ratio <- liter_nanog_E4_df$Intersect / 1010  
ggplot(liter_nanog_E4_df %>% 
         dplyr::slice_min(n = 30, order_by = FDR),
       aes(y=reorder(TF,-FDR),x=ratio))+
  geom_point(aes(size=Intersect,color=FDR))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(FDR,size="Intersect"), 
       x="Gene Ratio",y="Transcription Factor",title="literature ChIP-X Enrichment (Nanog E4)")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  guides(size=guide_legend("Count"))
ggsave('figures/liter_nanog_E4.jpg',device='jpg', width = 8, height = 8)
# MGI enrichment (by Enrichr-KG)
mgi_nanog_E4_df <- read.csv('figures/Enrichr-KG_MGI_nanog_E4_1020.csv')
ggplot(mgi_nanog_E4_df %>%
         dplyr::slice_min(n = 30, order_by = q.value),
       aes(x=reorder(Term,-q.value),y=-log10(q.value))) +
  geom_bar(stat="identity") +
  labs(x="MGI term",y="-log10(qvalue)",title="MGI Enrichment (Oct4 E4)") + 
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))
ggsave('figures/mgi_nanog_E4.jpg',device='jpg', width = 12, height = 8)


###############################################################################
save(list = ls(), file = 'temp.RData')
load('temp.RData')
