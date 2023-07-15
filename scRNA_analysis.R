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
rm(ens.mm.v109)
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
                                  sapply(sort(seq(5,max(sce$total_counts),5)), 
                                         function(x){latex2exp::TeX(paste0('$',x,'\\times 10^6$'))})),
                       breaks = sort(c(0,seq(5,max(sce$total_counts),5)))) +
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
  theme(legend.position = c(0.95, 0.5))
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

# test for resolution of leiden and features number (can also test other parameters!)
sweep_para <- function(seurat, seed, features, res, feature_num){
  para_res <- c()
  sil_score <- c()
  cluster_num <- c()
  ft_num <- c()
  PC_num <- c()
  for(j in feature_num){
    print(paste0('feature number = ', j))
    seurat <- RunPCA(seurat, features = features[1:j], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml 
    PC_sele <- 15 #round(PC_num.gml$dim.est,0)
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
                  feature_num = c(100, 250, 500, seq(1000,5000,500)))
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

# comparison of different cluster results (different resolutions)
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
com_para_r <- compare_para(seurat, seed, c(0.3,0.6),c(4000, 4000))
table(com_para_r[[1]],com_para_r[[2]])
rm(com_para_r)

# select the best feature number and best PC number
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
ElbowPlot(seurat) +
  geom_vline(xintercept = PC_num.elbow, color = 'blue', linetype="dashed") +
  geom_vline(xintercept = PC_num.gml, color = 'coral', linetype="dashed") +
  annotate('text',x=PC_num.elbow, y=10,label='elbow point',color = 'red') +
  annotate('text',x=PC_num.gml, y=12,label='global maximum likelihood',color = 'red')
ggsave('figures/PC_num_sele.jpg',device='jpg', width = 10, height = 6)

# PC numbers selected by GML/elboe are not satisfying here, so arbitrarily change to 15
#PC_num.gml <- 15

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
res_sele <- 0.3
seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
seurat <- FindClusters(seurat, resolution = res_sele, algorithm = 4)  # algorithm 4 is "Leiden"; 1 is "Louvain"
DimPlot(seurat, reduction = "tsne", label = TRUE, shape.by = 'Cell.Type')

# pass the clustering result back to SCE object
cluster.leiden <- seurat@meta.data[[paste0('originalexp_snn_res.',res_sele)]]
colLabels(sce_dev) <- cluster.leiden

# plot marker genes of each cluster (heatmap)
group_marker <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
group_marker %>% dplyr::filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 10, order_by = p_val_adj) -> top10
DoHeatmap(seurat,features=top10$gene, slot = "scale.data") #+ scale_fill_virdis()
#FeaturePlot(seurat, reduction = 'TSNE', features=top10[top10['cluster']==1,]$gene[1:4])
leiden_markers <- plotHeatmap(sce_dev, exprs_values = "logcounts", 
                              order_columns_by=c("label", "Cell.Type"), 
                              features=top10$gene, cluster_rows = FALSE, center = TRUE,
                              gaps_col = cumsum(as.numeric(table(colLabels(sce_dev)))),
                              gaps_row = seq(0,50,10), 
                              color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40))
ggsave('figures/leiden_markers.jpg', leiden_markers, device='jpg', width = 10, height = 8)

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
detail_sil(res=c(0.3,0.6), seurat = seurat, PC_sele=PC_num.gml)

# clean up variable
rm(group_diff)
rm(group_marker)
rm(top10)
rm(group_diff_serum_2)
rm(leiden_markers)
rm(up_reg)
rm(down_reg)
rm(sce_seurat)