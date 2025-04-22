library(Matrix); library(stringr); library(dplyr)
library(readr); library(here); 
library(DESeq2); library(reshape2); library("vsn")
library(ggplot2); library(ggpubr); library(ggrepel); library(RColorBrewer); library(patchwork)

setwd("/fs/cbsuvlaminck3/workdir/mm2937/ActivatedPBMC/")


counts2 <- read.csv("Unstranded_CountMatrix2_STAR_GRCh38.csv", sep = ",", header = TRUE, row.names = 1)
# counts <- read.csv("Unstranded_CountMatrix_removeDupl.csv", sep = ",", header = TRUE, row.names = 1)
colnames(counts2) = str_split_fixed(colnames(counts2), "_", 5)[,4]

counts = counts2
# counts = cbind(counts1, counts2)

# meta.data <- read.csv("metadata_corrected.csv", sep = ",", header = TRUE)
meta.data <- read.csv("metadata_corrected_updated.csv", sep = ",", header = TRUE)
meta.data

row.names(meta.data) <- meta.data$Sample.ID.to.core.lab

meta.data = meta.data[colnames(counts),]
counts = counts[,rownames(meta.data)]

# Calculate total genes detected and total counts
meta.data$TotalCounts = colSums(counts)
meta.data$GenesDetected = colSums(counts > 0)


# Rename the conditions
meta.data$Condition.of.Experiment <- as.character(meta.data$Condition.of.Expt)
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells alone at 0 hour"] <- "PBMCs (0 hours)"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells alone at 16 hours"] <- "PBMCs (16 hours)"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells+DMSO"] <- "PBMCs (Control)"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells+BTP2"] <- "PBMCs + BTP2"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells+PHA"] <- "PBMCs + PHA"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells+PHA+BTP2"] <- "PBMCs + PHA + BTP2"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells+CM4620"] <- "PBMCs + CM4620"
meta.data$Condition.of.Experiment[meta.data$Condition.of.Experiment == "Cells+PHA+CM4620"] <- "PBMCs + PHA + CM4620"
table(meta.data$Condition.of.Experiment)


meta.data$Condition.of.Experiment = factor(meta.data$Condition.of.Experiment, levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + CM4620", "PBMCs + PHA", "PBMCs + PHA + BTP2", "PBMCs + PHA + CM4620"))

# Remove extra rows from counts and add them to meta data
meta.data <- cbind(meta.data, t(counts[1:4,]))
counts = counts[5:nrow(counts),]

genenames = read.csv(file = "./geneid_genename_map.txt", sep = "\t", header = FALSE, col.names = c("geneid", "genename"), row.names = 1)
genenames$genename = make.names(genenames$genename, unique = T)
rownames(counts) <- genenames[rownames(counts),]

counts <- counts[,row.names(meta.data)]
table(meta.data$Condition.of.Experiment)

meta.data$S.no <- NULL
meta.data$Sample.ID <- meta.data$Sample.ID.to.core.lab

table(meta.data$Expt..)
meta.data$DonorID <- as.character(meta.data$DonorID)
meta.data$DonorID[meta.data$Expt.. == "BTP2#9"] <- "Donor1"
meta.data$DonorID[meta.data$Expt.. == "BTP2#10"] <- "Donor2"
meta.data$DonorID[meta.data$Expt.. == "BTP2#11"] <- "Donor3"
meta.data$DonorID[meta.data$Expt.. == "BTP2#12"] <- "Donor4"
meta.data$DonorID = factor(meta.data$DonorID, levels = c("Donor1", "Donor2", "Donor3", "Donor4"))

table(meta.data$Condition.of.Expt)
condition_levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + CM4620", "PBMCs + PHA", "PBMCs + PHA + BTP2", "PBMCs + PHA + CM4620")
meta.data$Condition.of.Experiment = factor(meta.data$Condition.of.Experiment, levels = condition_levels)

write.csv(meta.data, "./csvsDEseq2/metadata_full_export.csv")
write.csv(counts, "./csvsDEseq2/Unstranded_CountMatrix_full_export.csv")


pdf("plots2/data_qc.pdf", height = 2, width = 5)
(ggplot(data = meta.data) + geom_jitter(mapping = aes(x = DonorID, y = TotalCounts), width = 0.0) + geom_boxplot(mapping = aes(x = DonorID, y = TotalCounts)) | ggplot(data = meta.data) + geom_jitter(mapping = aes(x = DonorID, y = GenesDetected), width = 0.0) + geom_boxplot(mapping = aes(x = DonorID, y = GenesDetected))) & theme_classic(base_size = 6) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black")) & xlab("Donor ID")
dev.off()

pdf("plots2/data_qc2.pdf", height = 2.5, width = 5)
(ggplot(data = meta.data) + geom_jitter(mapping = aes(x = Condition.of.Experiment, y = TotalCounts), width = 0.0) + geom_boxplot(mapping = aes(x = Condition.of.Experiment, y = TotalCounts)) | ggplot(data = meta.data) + geom_jitter(mapping = aes(x = Condition.of.Experiment, y = GenesDetected), width = 0.0) + geom_boxplot(mapping = aes(x = Condition.of.Experiment, y = GenesDetected))) & theme_classic(base_size = 6)  & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0)) & xlab("Condition")
dev.off()



dataset <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = ~Condition.of.Experiment)

# Count normalization and Wald statistical tests
dataset = estimateSizeFactors(dataset) #Sequencing depth normalization (calculation of size factors)
dataset = estimateDispersions(dataset) #Estimation of dispersion parameters for the negative binomial distribution
keep <- rowSums(counts(dataset)) >= 10
dataset <- dataset[keep,]
levels(dataset$Condition.of.Experiment)
dataset$Condition.of.Experiment <- as.factor(dataset$Condition.of.Experiment)
dataset$Condition.of.Experiment <- droplevels(dataset$Condition.of.Experiment)

# Extracting transformed values
vsd <- vst(dataset, blind=FALSE)
# rld <- rlog(dataset, blind=FALSE)
head(assay(vsd), 3)

# Effects of transformations on the variance
# this gives log2(n + 1)
ntd <- normTransform(dataset)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

# Principal component plot of the samples
pcaData <- plotPCA(ntd, intgroup=c("DonorID", "Condition.of.Experiment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, mapping = aes(PC1, PC2, color=DonorID, shape=Condition.of.Experiment, label = name)) + geom_point(size = 3) + geom_text_repel()+ xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_classic(base_size = 16)

ggplot(pcaData, mapping = aes(PC1, PC2, color=Condition.of.Experiment, shape=DonorID)) + geom_point(size = 3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_classic(base_size = 16)

pdf("./plots2/pca_beforebc.pdf", width = 3.8, height = 2.0)
ggplot(pcaData, mapping = aes(PC1, PC2, color = Condition.of.Experiment, shape=DonorID)) + geom_point(size = 2) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))  & theme_classic(base_size = 8) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt")) & scale_color_manual(values =condition_colors)
dev.off()

brewer.pal(name = "Dark2", n = 6)
condition_colors = c('#E78AC3','#FFD92F','#377EB8','#7570B3','#E41A1C','#4DAF4A')


## Data quality assessment by sample clustering and visualization
library(limma)
## Removing batch effect using Limma
mat <- assay(vsd)
mm <- model.matrix(~Condition.of.Experiment, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$DonorID, design=mm)
assay(vsd) <- mat
pcaData <- plotPCA(vsd, intgroup=c("DonorID", "Condition.of.Experiment", "S.no"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, mapping = aes(PC1, PC2, color=DonorID, shape=Condition.of.Experiment)) + geom_point(size = 3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_classic(base_size = 16)

pdf("./plots2/pca_afterbc.pdf", width = 3.5, height = 2.0)
ggplot(pcaData, mapping = aes(PC1, PC2, color=Condition.of.Experiment, shape=DonorID)) + geom_point(size = 2)  + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) & theme_classic(base_size = 8) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt")) & scale_color_manual(values = condition_colors)
dev.off()



# Heatmap of the count matrix
# 
select <- order(rowMeans(counts(dataset,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(pcaData[,c("Condition.of.Experiment","DonorID")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


sampleDists <- cor(assay(vsd))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)
colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255)
colors

df <- as.data.frame(pcaData[,c("Condition.of.Experiment","DonorID")])
annotation_df = df
annotation_df$DonorID <- as.factor(annotation_df$DonorID)

length(unique(annotation_df$Condition.of.Experiment))

condition_colors = brewer.pal(name = "Dark2", n = length(unique(annotation_df$Condition.of.Experiment)))
condition_colors = c('#E78AC3','#FFD92F','#377EB8','#7570B3','#E41A1C','#4DAF4A')

names(condition_colors) <- levels(annotation_df$Condition.of.Experiment)
donor_colors = brewer.pal(name = "Set1", n = length(unique(annotation_df$DonorID)))
donor_colors = brewer.pal(name = "Set2", n = length(unique(annotation_df$DonorID)))

names(donor_colors) <- levels(annotation_df$DonorID)
my_colour = list(
  Condition.of.Experiment = condition_colors,
  DonorID = donor_colors
)

library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", annotation_names_row = F, show_rownames = F, 
         show_colnames = F, treeheight_row = 0, treeheight_col = 5,
         col=colors, annotation_col =  annotation_df, annotation_colors = my_colour,
         fontsize = 8)

pheatmap(sampleDistMatrix,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation", annotation_names_row = F, show_rownames = F, 
         show_colnames = F, treeheight_row = 0, treeheight_col = 5,
         col=colors, annotation_col =  annotation_df, annotation_colors = my_colour,
         fontsize = 6, height =  2.0, width = 3.5, filename = "./plots2/heatmap.pdf")
dev.off()

dataset$Condition.of.Experiment
dataset$Condition.of.Experiment <- relevel(dataset$Condition.of.Experiment, ref = "PBMCs (Control)")
table((dataset$Condition.of.Experiment))

# Running DESeq2
dataset <- DESeq(dataset)
resultsNames(dataset)

# Make a boxplot of the Cooks distances to see if one sample is consistently higher than others
par(mar=c(8,5,2,2))
boxplot(log10(assays(dataset)[["cooks"]]), range=0, las=2)
# cooks_data = melt(log10(assays(dataset)[["cooks"]]), varnames = "Sample",value.name = "CooksDistance")
# ggplot(data = cooks_data) + geom_boxplot(aes(Sample, CooksDistance))

save(dataset, file = "./robjs/deseq_object2.RObj")
load("./robjs/deseq_object2.RObj")

table(dataset$Condition.of.Experiment)
resultsNames(dataset) 

dataset_subset = dataset

table(dataset$Condition.of.Experiment)
resultsNames(dataset) 

# DESeq2 differential expression analysis for effects of BTP2
dataset_results = results(dataset_subset, contrast=c('Condition.of.Experiment',"PBMCs + BTP2","PBMCs (Control)"))
dataset_results = dataset_results[complete.cases(dataset_results),] #removes any rows with NA
summary(dataset_results)
dataset_results = as.data.frame(dataset_results[order(dataset_results$padj),]) #orders by adjusted p-value
dataset_results["STIM1",]

write.csv(dataset_results, file = "csvsDEseq2/BTP2-DMSO.csv")
resultsNames(dataset_subset)
resLFC <- lfcShrink(dataset_subset, coef="Condition.of.Experiment_PBMCs...BTP2_vs_PBMCs..Control.", type="apeglm")
summary(resLFC)
resLFC = as.data.frame(resLFC[order(resLFC$padj),]) #orders by adjusted p-value
write.csv(resLFC, file = "csvsDEseq2/BTP2-DMSO_shrunk.csv")

# DESeq2 differential expression analysis for effects of CM4620
dataset_results = results(dataset_subset, contrast=c('Condition.of.Experiment',"PBMCs + CM4620","PBMCs (Control)"))
dataset_results = dataset_results[complete.cases(dataset_results),] #removes any rows with NA
summary(dataset_results)
dataset_results = as.data.frame(dataset_results[order(dataset_results$padj),]) #orders by adjusted p-value
dataset_results["STIM1",]

write.csv(dataset_results, file = "csvsDEseq2/CM4620-DMSO.csv")
resultsNames(dataset_subset)
resLFC <- lfcShrink(dataset_subset, coef="Condition.of.Experiment_PBMCs...CM4620_vs_PBMCs..Control.", type="apeglm")
summary(resLFC)
resLFC = as.data.frame(resLFC[order(resLFC$padj),]) #orders by adjusted p-value
write.csv(resLFC, file = "csvsDEseq2/CM4620-DMSO_shrunk.csv")


# DESeq2 differential expression analysis for effects of PHA
dataset_results = results(dataset_subset, contrast=c('Condition.of.Experiment','PBMCs + PHA',"PBMCs (Control)"))
dataset_results = dataset_results[complete.cases(dataset_results),] #removes any rows with NA
summary(dataset_results)
dataset_results = as.data.frame(dataset_results[order(dataset_results$padj),]) #orders by adjusted p-value
dataset_results["STIM1",]

write.csv(dataset_results, file = "csvsDEseq2/PHA-DMSO.csv")
resultsNames(dataset_subset)
resLFC <- lfcShrink(dataset_subset, coef="Condition.of.Experiment_PBMCs...PHA_vs_PBMCs..Control.", type="apeglm")
summary(resLFC)
resLFC = as.data.frame(resLFC[order(resLFC$padj),]) #orders by adjusted p-value
write.csv(resLFC, file = "csvsDEseq2/PHA-DMSO_shrunk.csv")

###########################################################################
########################################################################### 

dataset_subset$Condition.of.Experiment <- droplevels(dataset_subset$Condition.of.Experiment)
dataset_subset$Condition.of.Experiment <- relevel(dataset_subset$Condition.of.Experiment, ref = 'PBMCs + PHA')
dataset_subset <- DESeq(dataset_subset)
resultsNames(dataset_subset)

dataset_results = results(dataset_subset, contrast=c('Condition.of.Experiment',"PBMCs + PHA + BTP2", 'PBMCs + PHA'))
dataset_results = dataset_results[complete.cases(dataset_results),] #removes any rows with NA
summary(dataset_results)
dataset_results = as.data.frame(dataset_results[order(dataset_results$padj),]) #orders by adjusted p-value
write.csv(dataset_results, file = "csvsDEseq2/BTP2-PHA.csv")
resultsNames(dataset_subset)
resLFC <- lfcShrink(dataset_subset, coef="Condition.of.Experiment_PBMCs...PHA...BTP2_vs_PBMCs...PHA", type="apeglm")
resLFC = resLFC[complete.cases(resLFC),] #removes any rows with NA
summary(resLFC)
resLFC = as.data.frame(resLFC[order(resLFC$padj),]) #orders by adjusted p-value
write.csv(resLFC, file = "csvsDEseq2/BTP2-PHA_shrunk.csv")

# DESeq2 differential expression analysis for effects of CM4620 on activation with PHA
dataset_results = results(dataset_subset, contrast=c('Condition.of.Experiment',"PBMCs + PHA + CM4620", 'PBMCs + PHA'))
dataset_results = dataset_results[complete.cases(dataset_results),] #removes any rows with NA
summary(dataset_results)
dataset_results = as.data.frame(dataset_results[order(dataset_results$padj),]) #orders by adjusted p-value
write.csv(dataset_results, file = "csvsDEseq2/CM4620-PHA.csv")

resultsNames(dataset_subset)
resLFC <- lfcShrink(dataset_subset, coef="Condition.of.Experiment_PBMCs...PHA...CM4620_vs_PBMCs...PHA", type="apeglm")
resLFC = resLFC[complete.cases(resLFC),] #removes any rows with NA
summary(resLFC)
resLFC = as.data.frame(resLFC[order(resLFC$padj),]) #orders by adjusted p-value
write.csv(resLFC, file = "csvsDEseq2/CM4620-PHA_shrunk.csv")

##########################################################################
########################################################################### 

dataset_subset$Condition.of.Experiment <- droplevels(dataset_subset$Condition.of.Experiment)
dataset_subset$Condition.of.Experiment <- relevel(dataset_subset$Condition.of.Experiment, ref = 'PBMCs + PHA + BTP2')
dataset_subset <- DESeq(dataset_subset)
resultsNames(dataset_subset)

dataset_results = results(dataset_subset, contrast=c('Condition.of.Experiment',"PBMCs + PHA + CM4620", 'PBMCs + PHA + BTP2'))
dataset_results = dataset_results[complete.cases(dataset_results),] #removes any rows with NA
summary(dataset_results)
dataset_results = as.data.frame(dataset_results[order(dataset_results$padj),]) #orders by adjusted p-value
write.csv(dataset_results, file = "csvsDEseq2/CM4620-BTP2.csv")
resultsNames(dataset_subset)
resLFC <- lfcShrink(dataset_subset, coef="Condition.of.Experiment_PBMCs...PHA...CM4620_vs_PBMCs...PHA...BTP2", type="apeglm")
resLFC = resLFC[complete.cases(resLFC),] #removes any rows with NA
summary(resLFC)
resLFC = as.data.frame(resLFC[order(resLFC$padj),]) #orders by adjusted p-value
write.csv(resLFC, file = "csvsDEseq2/CM4620-BTP2_shrunk.csv")

colData(dataset_subset)
gene_expression_df <- as.data.frame(t(counts(dataset_subset, normalized=TRUE)))
metadata_df <- as.data.frame(colData(dataset_subset))
metadata_df = metadata_df[c("Condition.of.Experiment")]
colnames(metadata_df) <- c("Condition")
gene_expression_df = cbind(gene_expression_df, metadata_df)
mean_expression <- aggregate(. ~ Condition, data = gene_expression_df, FUN = mean)
rownames(mean_expression)  = mean_expression$Condition
mean_expression$Condition = NULL
mean_expression = t(mean_expression)
mean_expression_condition = as.data.frame(mean_expression)
colnames(mean_expression_condition) <- paste0("baseMean: ", colnames(mean_expression_condition))
write.csv(mean_expression_condition, "./csvsDEseq2/baseMeanCounts_Conditions.csv")
dim(mean_expression_condition)


temp_dgea = read.csv("csvsDEseq2/CM4620-BTP2_shrunk.csv", row.names = 1)
dim(temp_dgea)
temp_dgea <- na.omit(temp_dgea)
dim(temp_dgea)
write.csv(cbind(mean_expression_condition[rownames(temp_dgea),], temp_dgea), "./csvsDEseq2/CM4620-BTP2_shrunk_withbasemeans.csv")
###########################################################################
########################################################################### 

#########################
library(EnhancedVolcano)
volcano_input <- read.csv(file = "./csvsDEseq2/CM4620-BTP2_shrunk.csv", row.names = 1)
colnames(volcano_input)
volcano_input <- na.omit(volcano_input)

pdf("plots2/CM4620_BTP2_volcano.pdf", height = 2.2, width = 3.0)
EnhancedVolcano(volcano_input,
                lab = rownames(volcano_input),
                x = 'log2FoldChange',
                y = 'padj',
                title = "",
                pCutoff = 10e-2, pointSize = 0.5,
                FCcutoff = 2.0, 
                subtitle = NULL, axisLabSize = 6, titleLabSize = 12, labSize = 2.0, # xlim = c(-10,13),
                legendPosition = "none", caption = NULL) & theme_classic(base_size = 6) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt"))
dev.off()


volcano_input <- read.csv(file = "./csvsDEseq2/PHA-DMSO_shrunk.csv", row.names = 1)
colnames(volcano_input)
volcano_input <- na.omit(volcano_input)

pdf("plots2/PHA_DMSO_volcano.pdf", height = 2.0, width = 3.4)
EnhancedVolcano(volcano_input,
                lab = rownames(volcano_input),
                x = 'log2FoldChange',
                y = 'padj',
                title = "",
                pCutoff = 10e-2, pointSize = 0.5,
                FCcutoff = 2.0, 
                subtitle = NULL, axisLabSize = 6, titleLabSize = 12, labSize = 2.0, xlim = c(-13.5,14.0),
                legendPosition = "none", caption = NULL) & theme_classic(base_size = 6) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt")) 
dev.off()


volcano_input <- read.csv(file = "./csvsDEseq2/BTP2-PHA_shrunk.csv", row.names = 1)
dim(volcano_input)
colnames(volcano_input)
volcano_input <- na.omit(volcano_input)

pdf("plots2/BTP2_PHA_volcano.pdf", height = 2.0, width = 3.4)
EnhancedVolcano(volcano_input,
                lab = rownames(volcano_input),
                x = 'log2FoldChange',
                y = 'padj',
                title = "",
                pCutoff = 10e-2, pointSize = 0.5,
                FCcutoff = 2.0, 
                subtitle = NULL, axisLabSize = 6, titleLabSize = 12, labSize = 2.0, xlim = c(-7.0,12.0),
                legendPosition = "none", caption = NULL) & theme_classic(base_size = 6) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt"))
dev.off()


volcano_input <- read.csv(file = "./csvsDEseq2/CM4620-PHA_shrunk.csv", row.names = 1)
dim(volcano_input)
colnames(volcano_input)
volcano_input <- na.omit(volcano_input)

pdf("plots2/CM4620_PHA_volcano.pdf", height = 2.0, width = 3.4)
EnhancedVolcano(volcano_input,
                lab = rownames(volcano_input),
                x = 'log2FoldChange',
                y = 'padj',
                title = "",
                pCutoff = 10e-2, pointSize = 0.5,
                FCcutoff = 2.0, 
                subtitle = NULL, axisLabSize = 6, titleLabSize = 12, labSize = 2.0, xlim = c(-8.0,12.0),
                legendPosition = "none", caption = NULL) & theme_classic(base_size = 6) & theme(plot.margin  = margin(0,0,0,0), axis.text = element_text(color = "black"), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt"))
dev.off()

volcano_input <- read.csv(file = "./csvsDEseq2/BTP2-PHA_shrunk.csv", row.names = 1)
combined_dgea_coor1 = data.frame(BTP2_lfc = volcano_input$log2FoldChange, BTP2_adj = volcano_input$padj, row.names = rownames(volcano_input))
volcano_input <- read.csv(file = "./csvsDEseq2/CM4620-PHA_shrunk.csv", row.names = 1)
combined_dgea_coor2 = data.frame(CM4620_lfc = volcano_input$log2FoldChange, CM4620_padj = volcano_input$padj, row.names = rownames(volcano_input))

combined_dgea_coor = merge(combined_dgea_coor1, combined_dgea_coor2, by = 'row.names', all = TRUE)
head(combined_dgea_coor)

combined_dgea_coor_sig = combined_dgea_coor[(combined_dgea_coor$BTP2_adj < 0.01) | (combined_dgea_coor$CM4620_padj < 0.01),]

(combined_dgea_coor$CM4620_lfc - combined_dgea_coor$BTP2_lfc) >= 1.0
combined_dge$label = NA

combined_dgea_coor_sig$label <- NA
combined_dgea_coor_sig$label <- ifelse(combined_dgea_coor_sig$CM4620_lfc - combined_dgea_coor_sig$BTP2_lfc > 2.0, combined_dgea_coor$Row.names, NA)
combined_dgea_coor_sig$label <- ifelse(combined_dgea_coor_sig$CM4620_lfc - combined_dgea_coor_sig$BTP2_lfc < -2.0, combined_dgea_coor$Row.names, combined_dgea_coor_sig$label)
head(combined_dgea_coor_sig)

pdf("plots2/BTP2_CM4620_PHA_volcano_scatter.pdf", height = 2.5, width = 3.0)
library(ggrepel)
ggplot(data = combined_dgea_coor_sig, mapping = aes(x = BTP2_lfc, y = CM4620_lfc)) + geom_point(size = 0.5, alpha = 0.5) & stat_cor(method="pearson") & theme_classic(base_size = 6) & theme(plot.margin  = margin(0,0,0,0), text = element_text(color = "black", size = 6), axis.text = element_text(color = "black"), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt"), legend.key.size = unit(2, "pt")) & xlab("Log2 fold change: PHA vs. PHA + BTP2") & ylab("Log2 fold change: PHA vs. PHA + CM4620") & geom_abline(intercept = -2, slope = 1, color="red", 
              linetype="dashed", size=0.5) & geom_abline(intercept = 2, slope = 1, color="red", 
              linetype="dashed", size=0.5) & geom_text_repel(data = combined_dgea_coor_sig, mapping = aes(x = BTP2_lfc, y = CM4620_lfc, label = label), size = 2, max.overlaps = 20) 
dev.off()

write.csv(combined_dgea_coor_sig, file = "./csvsDEseq2/BTP2_CM4620_combined_DGEA.csv")

# c(-5.8,6)
# c(-10,13)
# Plot gene expression data
d <- plotCounts(dataset_subset, gene=rownames(volcano_input)[which.min(volcano_input$padj)], intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID)) + geom_point(position=position_jitter(w=0.0,h=0)) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1.0))

# + stat_compare_means(data = d, method = "wilcox.test", comparisons = list(c("PBMCs (Control)", "PBMCs + PHA"), c("PBMCs + PHA", "PBMCs + PHA + BTP2")), size = 2, label = "p.signif", tip.length = 0.01, step.increase = 0.02, vjust = 0.5)

table(d$Condition.of.Experiment)


condition_colors = brewer.pal(name = "Dark2", n = length(unique(annotation_df$Condition.of.Experiment)))
condition_colors = c('#E78AC3','#FFD92F','#377EB8','#7570B3','#E41A1C','#4DAF4A')
condition_levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + CM4620", "PBMCs + PHA", "PBMCs + PHA + BTP2", "PBMCs + PHA + CM4620")
names(condition_colors) <- condition_levels

gene_name = c("ORAI1")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("ORAI2")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2


gene_name = c("ORAI3")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p3 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p3

gene_name = c("STIM1")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels =condition_levels)
table(d$Condition.of.Experiment)
p4 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p4

gene_name = c("STIM2")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p5 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p5


pdf("./plots2/Orai_Stim_plots_arr2.pdf", width = 5, height = 1.5)
((p4 | p5 | p1 | p2 | p3)) & theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()


gene_name = c("IL2")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("IL2RA")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2


pdf("./plots2/il2_il2ra_plots.pdf", width = 2, height = 1.5)
(p1 | p2)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()

gene_name = c("CD27")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("CD28")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2


pdf("./plots2/cd27_cd28_plots.pdf", width = 2, height = 1.5)
(p1 | p2)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()

gene_name = c("PRF1")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("GZMB")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2

pdf("./plots2/prf1_gzma_cd96_plots.pdf", width = 2.0, height = 1.5)
(p1 | p2)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()

gene_name = c("IFNG")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p0 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p0

gene_name = c("CXCL9")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels =condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("CXCL10")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2

gene_name = c("CXCR3")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p3 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p3

pdf("./plots2/cxcl9_10.pdf", width = 4, height = 1.5)
(p0 | p1 | p2 | p3)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()

gene_name = c("CD3E")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1


gene_name = c("CD4")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2

gene_name = c("CD8A")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p3 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p3



pdf("./plots2/co-receptors.pdf", width = 3, height = 1.5)
(p1 | p2 | p3)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()


gene_name = c("ITGAE")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p3 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p3

pdf("./plots2/contact.pdf", width = 2.0, height = 1.5)
( p3 | p3) & theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()


gene_name = c("TGFB1")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("IL10")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2

gene_name = c("FOXP3")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p3 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p3


gene_name = c("CTLA4")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = condition_levels)
table(d$Condition.of.Experiment)
p4 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name)  & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p4

pdf("./plots2/tolerance.pdf", width = 4, height = 1.5)
(p1 | p2| p3 | p4)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()


###################################################################################################
load(file = "robjs/deseq_object_corrected.Robj")
meta.data <- meta.data[!(meta.data$Expt.. %in% c("BTP2#3")),]

dim(dataset_subset@assays@data$counts)
dim(dataset_subset)


gene_name = c("CD3D")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + PHA", "PBMCs + PHA + BTP2"))
table(d$Condition.of.Experiment)
p1 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name) + stat_compare_means(data = d, mapping = aes(x=Condition.of.Experiment, y=count), inherit.aes = F, label.y = log10(max(d$count) * 1.20), size = 2) & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p1

gene_name = c("CD3E")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + PHA", "PBMCs + PHA + BTP2"))
table(d$Condition.of.Experiment)
p2 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name) + stat_compare_means(data = d, mapping = aes(x=Condition.of.Experiment, y=count), inherit.aes = F, label.y = log10(max(d$count) * 1.20), size = 2) & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p2


gene_name = c("CD3G")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + PHA", "PBMCs + PHA + BTP2"))
table(d$Condition.of.Experiment)
p3 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name) + stat_compare_means(data = d, mapping = aes(x=Condition.of.Experiment, y=count), inherit.aes = F, label.y = log10(max(d$count) * 1.20), size = 2) & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "none") & scale_color_manual(values = condition_colors)
p3


gene_name = c("CD247")
d <- plotCounts(dataset_subset, gene=gene_name, intgroup=c("Condition.of.Experiment", "DonorID"), returnData = T)
d$Condition.of.Experiment <- factor(d$Condition.of.Experiment, levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + PHA", "PBMCs + PHA + BTP2"))
table(d$Condition.of.Experiment)
p4 <- ggplot(d, aes(x=Condition.of.Experiment, y=count, shape = DonorID, color=Condition.of.Experiment)) + geom_point(position=position_jitter(w=0.0,h=0), size = 1) + scale_y_log10() + ylab(gene_name) + stat_compare_means(data = d, mapping = aes(x=Condition.of.Experiment, y=count), inherit.aes = F, label.y = log10(max(d$count) * 1.20), size = 2) & theme_classic(base_size = 6) & theme(plot.margin  = margin(3,3,3,3), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1.0), legend.position = "right") & scale_color_manual(values = condition_colors)
p4


pdf("./plots2/cd3.pdf", width = 5.0, height = 1.5)
(p1 | p2 | p3 | p4)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()

