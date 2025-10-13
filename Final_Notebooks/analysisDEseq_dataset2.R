library(Matrix); library(stringr); library(dplyr)
library(readr); library(here); 
library(DESeq2); library(reshape2); library("vsn")
library(ggplot2); library(ggpubr); library(ggrepel); library(RColorBrewer); library(patchwork)

setwd("/workdir/mm2937/ActivatedPBMC/")

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


condition_colors = brewer.pal(name = "Dark2", n = length(unique(d$Condition.of.Experiment)))
condition_colors = c('#E78AC3','#FFD92F','#377EB8','#7570B3','#E41A1C','#4DAF4A')
condition_levels = c("PBMCs (Control)", "PBMCs + BTP2", "PBMCs + CM4620", "PBMCs + PHA", "PBMCs + PHA + BTP2", "PBMCs + PHA + CM4620")
names(condition_colors) <- condition_levels

# --- Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)

# --- Inputs you already have / want ---
genes_keep <- c("CD3E", "CD4", "CD8A",
                "CD27","CD28",
                "ITGAE",
                "IL2","IL2RA",
                "IFNG","CXCL9","CXCL10","CXCR3",
                "PRF1","GZMB",
                "TGFB1","IL10","FOXP3","CTLA4")

condition_levels <- c(
  "PBMCs (Control)",
  "PBMCs + BTP2",
  "PBMCs + CM4620",
  "PBMCs + PHA",
  "PBMCs + PHA + BTP2",
  "PBMCs + PHA + CM4620"
)

condition_colors <- c(
  "PBMCs (Control)"       = "#E78AC3",
  "PBMCs + BTP2"          = "#FFD92F",
  "PBMCs + CM4620"        = "#377EB8",
  "PBMCs + PHA"           = "#7570B3",
  "PBMCs + PHA + BTP2"    = "#E41A1C",
  "PBMCs + PHA + CM4620"  = "#4DAF4A"
)

comp_use <- list(
  c("PBMCs (Control)", "PBMCs + BTP2"),
  c("PBMCs (Control)", "PBMCs + CM4620"),
  c("PBMCs (Control)", "PBMCs + PHA"),
  c("PBMCs + PHA", "PBMCs + PHA + BTP2"),
  c("PBMCs + PHA", "PBMCs + PHA + CM4620")
)
# --- 1) Get counts for all genes ---
# 1) Collect counts for all genes
get_counts_one_gene <- function(g) {
  d <- plotCounts(
    dataset_subset,
    gene = g,
    intgroup = c("Condition.of.Experiment", "DonorID"),
    returnData = TRUE
  )
  d$Gene <- g
  d
}

df_all <- bind_rows(lapply(genes_keep, get_counts_one_gene)) %>%
  mutate(
    Gene = factor(Gene, levels = genes_keep),
    Condition.of.Experiment = factor(Condition.of.Experiment, levels = condition_levels),
    count_plot = ifelse(count <= 0, NA_real_, count)  # log-safe
  )

# 2) Stats per gene (unpaired Wilcoxon on donor-level summaries), BH per gene
agg_method <- "median"  # or "mean"
summ_fun <- if (agg_method == "median") median else mean

get_stats_for_gene <- function(gdat) {
  # donor-level summary (median or mean as you set above)
  d_sum <- gdat %>%
    dplyr::group_by(DonorID, Condition.of.Experiment) %>%
    dplyr::summarise(value = summ_fun(count, na.rm = TRUE), .groups = "drop")
  
  out <- lapply(comp_use, function(pr) {
    # pivot to ensure rows are matched by DonorID and both conditions are present
    wide <- d_sum %>%
      dplyr::filter(Condition.of.Experiment %in% pr) %>%
      tidyr::pivot_wider(
        id_cols   = DonorID,
        names_from = Condition.of.Experiment,
        values_from = value
      ) %>%
      # keep only donors with both values
      dplyr::filter(is.finite(.data[[pr[1]]]), is.finite(.data[[pr[2]]]))
    
    if (nrow(wide) < 2L) {  # need at least 2 pairs for a sensible paired test
      return(tibble(group1 = pr[1], group2 = pr[2], p = NA_real_))
    }
    
    x <- wide[[pr[1]]]
    y <- wide[[pr[2]]]
    diffs <- x - y
    
    # exact Wilcoxon is not available with ties/zeros; detect and fall back
    has_ties_or_zeros <- any(diffs == 0) || (length(unique(abs(diffs))) < length(diffs))
    exact_flag <- if (has_ties_or_zeros) FALSE else TRUE
    
    pval <- tryCatch(
      stats::wilcox.test(x, y, paired = TRUE, exact = exact_flag,
                         conf.int = TRUE, conf.level = 0.95)$p.value,
      error = function(e) NA_real_
    )
    tibble(group1 = pr[1], group2 = pr[2], p = pval)
  })
  
  stats <- dplyr::bind_rows(out)
  stats$p_adj <- dplyr::if_else(is.finite(stats$p), p.adjust(stats$p, method = "BH"), NA_real_)
  stats$p_label <- dplyr::case_when(
    is.na(stats$p_adj)   ~ "n/a",
    stats$p_adj < 1e-4   ~ "****",
    stats$p_adj < 1e-3   ~ "***",
    stats$p_adj < 1e-2   ~ "**",
    stats$p_adj < 5e-2   ~ "*",
    TRUE                 ~ "ns"
  )
  stats
}

stats_all <- df_all %>%
  dplyr::group_by(Gene) %>%
  dplyr::group_modify(~ get_stats_for_gene(.x)) %>%
  dplyr::ungroup()

stats_all
# --- 3) Measure panel height in inches from a base plot ---
p0 <- ggplot(df_all, aes(x = Condition.of.Experiment, y = count_plot,
                         color = Condition.of.Experiment)) +
  geom_point(aes(shape = DonorID),
             position = position_jitter(width = 0.12, height = 0),
             size = 2, alpha = 1) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 6) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  coord_trans(y = "log10", clip = "off") +
  scale_y_continuous(
    breaks = c(10,50,100,500,1000,5000,10000,50000,100000,500000,1000000),
    labels = scales::label_number(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 8) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line.x   = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.margin = margin(4, 6, 4, 6)
  )

gt <- ggplotGrob(p0)
panel_rows <- unique(gt$layout$t[grepl("^panel", gt$layout$name)])
panel_height_in <- convertHeight(gt$heights[panel_rows[1]], "in", valueOnly = TRUE)
panel_height_in <- ifelse(is.finite(panel_height_in) && panel_height_in > 0, panel_height_in, 2) # fallback

# --- 4) Robust per-facet ranges + convert 1 inch to log-data units ---
# Guards: ensure y_min < y_max and log_span > 0; supply fallbacks if degenerate
y_range <- df_all %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    y_min_raw = suppressWarnings(min(count_plot[count_plot > 0], na.rm = TRUE)),
    y_max_raw = suppressWarnings(max(count_plot, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_max = dplyr::if_else(is.finite(y_max_raw) & y_max_raw > 0, y_max_raw, 1),
    y_min = dplyr::if_else(is.finite(y_min_raw) & y_min_raw > 0 & y_min_raw < y_max,
                           y_min_raw, pmax(y_max/3, 1e-6)),
    log_span = pmax(log10(y_max) - log10(y_min), 1e-3),  # never zero
    # 1 inch in data-log units; enforce a minimum band height (e.g., 0.25 decades)
    delta_log_1in = pmax((1 / panel_height_in) * log_span, 0.25),
    log_cap  = log10(y_max),
    log_top  = log_cap + delta_log_1in
  )

# --- 5) Bracket positions: evenly spaced across the 1-inch band with padding ---
valid_x <- levels(df_all$Condition.of.Experiment)

pad_top_in <- 0.05   # keep 0.05 in below the top of the band
pad_bot_in <- 0.05   # keep 0.05 in above the data max
stats_all_draw <- stats_all %>%
  dplyr::filter(is.finite(p_adj), group1 %in% valid_x, group2 %in% valid_x) %>%
  dplyr::left_join(y_range, by = "Gene") %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(
    n_br = dplyr::n(),
    idx  = dplyr::row_number(),
    # convert inch paddings to log units using the same inches?log mapping
    log_pad_top = (pad_top_in / panel_height_in) * log_span,
    log_pad_bot = (pad_bot_in / panel_height_in) * log_span,
    band_start  = log_cap + log_pad_bot,
    band_end    = log_top - log_pad_top,
    # ensure the band has width; if not, widen to a minimal 0.1 decade
    band_end    = dplyr::if_else(band_end <= band_start, band_start + 0.1, band_end),
    frac        = dplyr::if_else(n_br == 1, 0.5, (idx - 1) / pmax(n_br - 1, 1)),
    y.position  = 10^(band_start + frac * (band_end - band_start))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(y.position) & y.position > 0)

# --- 6) Per-facet upper limits: just above the top bracket (with tiny pad) ---
df_limits <- stats_all_draw %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(y_upper = max(y.position, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::right_join(y_range, by = "Gene") %>%
  dplyr::mutate(
    # if a facet had no valid brackets, use top of the 1-inch band
    y_upper = dplyr::if_else(is.finite(y_upper), y_upper, 10^log_top)
  ) %>%
  dplyr::select(Gene, y_upper)

# --- 7) Final plot (jitter-only) with per-facet y limits via geom_blank ---
p <- ggplot(df_all, aes(x = Condition.of.Experiment, y = count_plot,
                        color = Condition.of.Experiment)) +
  geom_blank(data = df_limits, aes(y = y_upper), inherit.aes = FALSE) +
  geom_point(aes(shape = DonorID),
             position = position_jitter(width = 0.12, height = 0),
             size = 1, alpha = 1) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 6) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  coord_trans(y = "log10", clip = "off") +
  scale_y_continuous(
    breaks = c(10,50,100,500,1000,5000,10000,50000,100000,500000,1000000),
    labels = scales::label_number(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  ggpubr::stat_pvalue_manual(
    stats_all_draw,
    label = "p_label",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0, size = 2
  ) +
  theme_classic(base_size = 6) +   # set default text size to 6
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line.x   = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6),
    plot.margin = margin(4, 6, 4, 6)
  ) +
  labs(y = "bulk RNA count (log10)")

p

### ONLY PLOT/ NO STATS

p <- ggplot(df_all, aes(x = Condition.of.Experiment, y = count_plot,
                        color = Condition.of.Experiment)) +
  geom_point(aes(shape = DonorID),
             position = position_jitter(width = 0.12, height = 0),
             size = 1, alpha = 1) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 6) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  coord_trans(y = "log10", clip = "off") +
  scale_y_continuous(
    breaks = c(10,50,100,500,1000,5000,10000,50000,100000,500000,1000000),
    labels = scales::label_number(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  theme_classic(base_size = 6) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line.x   = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6),
    plot.margin = margin(4, 6, 4, 6)
  ) +
  labs(y = "")

p

pdf("./plots2/Bulk_RNA_all_NoStats.pdf", width = 6, height = 3.5)
(p)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
dev.off()


###################################################################################################
# Load qPCR CNVs and create plots

qPCR_cnv = read_csv('qPCR/CNVs_qPCR.csv')

library(tidyr)
df_long <- qPCR_cnv %>%
  pivot_longer(
    cols = starts_with("EXPT"),
    names_to = "DonorID",
    values_to = "count"
  ) %>%
  mutate(
    DonorID = recode(DonorID,
                     "EXPT 9"  = "Donor1",
                     "EXPT 10" = "Donor2",
                     "EXPT 11" = "Donor3",
                     "EXPT 12" = "Donor4"),
    Condition.of.Experiment = recode(Condition,
                                     "PBMC"                = "PBMCs (Control)",
                                     "PBMC + BTP2"         = "PBMCs + BTP2",
                                     "PBMC + CM4620"       = "PBMCs + CM4620",
                                     "PBMC + PHA"          = "PBMCs + PHA",
                                     "PBMC + PHA + BTP2"   = "PBMCs + PHA + BTP2",
                                     "PBMC + PHA + CM4620" = "PBMCs + PHA + CM4620"
    )
  ) %>%
  select(Gene, count, Condition.of.Experiment, DonorID)

# Needed packages
library(dplyr)
library(ggplot2)
library(ggpubr)   # for stat_pvalue_manual
library(rstatix)  # for wilcox_test / add_xy_position

# Defaults you asked for
condition_levels <- c(
  "PBMCs (Control)",
  "PBMCs + BTP2",
  "PBMCs + CM4620",
  "PBMCs + PHA",
  "PBMCs + PHA + BTP2",
  "PBMCs + PHA + CM4620"
)

comp_use <- list(
  c("PBMCs + PHA", "PBMCs + PHA + BTP2"),
  c("PBMCs + PHA", "PBMCs + PHA + CM4620"),
  c("PBMCs + PHA + BTP2", "PBMCs + PHA + CM4620"),
  c("PBMCs (Control)", "PBMCs + BTP2"),
  c("PBMCs (Control)", "PBMCs + CM4620"),
  c("PBMCs (Control)", "PBMCs + PHA")
)

condition_colors <- c('#E78AC3','#FFD92F','#377EB8','#7570B3','#E41A1C','#4DAF4A')
names(condition_colors) <- condition_levels

library(scales)

# Choose your genes and order them
genes_keep <- c("CD3E", "CD4", "CD8A",
                "CD27","CD28",
                "ITGAE",
                "IL2","CD25",
                "IFNG","CXCL9","CXCL10","CXCR3",
                "PRF1","GZMB",
                "TGFB1","IL10","FOXP3","CTLA4")   # put them in desired order

# Make sure we use dplyr verbs explicitly (avoids masking)
library(dplyr)

# 0) Ensure factors align to your palette
condition_levels <- names(condition_colors)
df_sub <- df_long %>%
  dplyr::filter(Gene %in% genes_keep) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = genes_keep),
    Condition.of.Experiment = factor(Condition.of.Experiment, levels = condition_levels)
  )

# 1) Count finite rows per gene (avoid masked count(); do it explicitly)
df_gene_counts <- df_sub %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(n_finite = sum(is.finite(count)), .groups = "drop")

# 2) Drop genes with zero finite rows, then drop unused facet levels
empty_genes <- as.character(df_gene_counts$Gene[df_gene_counts$n_finite == 0])
if (length(empty_genes)) {
  warning("Dropping empty genes (no finite counts): ", paste(empty_genes, collapse = ", "))
  df_sub <- df_sub %>% dplyr::filter(!(Gene %in% empty_genes))
}
df_sub <- droplevels(df_sub)

# 3) Recompute y caps & limits only for genes that remain
y_caps <- df_sub %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(y_cap = suppressWarnings(max(count, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::mutate(y_cap = ifelse(is.finite(y_cap), y_cap, 0))

band_frac <- 0.50  # +50% over max
min_band  <- 25    # absolute min band height (tune to your scale)
pad_top_frac <- 0.05
pad_bot_frac <- 0.05

df_limits <- y_caps %>%
  dplyr::mutate(
    band    = pmax(y_cap * band_frac, min_band),
    y_upper = y_cap + band
  )

# 4) (Re)compute stats and bracket positions for the surviving genes
#    ? if you already have stats_all, filter it to kept genes; else compute now.
get_stats <- function(gene_df, comps) {
  gene_df <- gene_df %>% dplyr::filter(is.finite(count))
  out <- lapply(comps, function(cc) {
    sub <- gene_df %>% dplyr::filter(Condition.of.Experiment %in% cc)
    grp_present <- table(sub$Condition.of.Experiment)
    if (length(grp_present) < 2 || any(grp_present[cc] < 1)) {
      return(tibble(group1 = cc[1], group2 = cc[2], p = NA_real_))
    }
    pval <- tryCatch(
      stats::wilcox.test(count ~ Condition.of.Experiment, data = sub, paired = FALSE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    tibble(group1 = cc[1], group2 = cc[2], p = pval)
  })
  stats <- dplyr::bind_rows(out) %>%
    dplyr::mutate(
      p_adj   = ifelse(is.finite(p), p.adjust(p, method = "BH"), NA_real_),
      p_label = dplyr::case_when(
        is.na(p_adj)   ~ "n/a",
        p_adj < 1e-4   ~ "****",
        p_adj < 1e-3   ~ "***",
        p_adj < 1e-2   ~ "**",
        p_adj < 5e-2   ~ "*",
        TRUE           ~ "ns"
      )
    )
  stats
}

stats_all <- df_sub %>%
  dplyr::group_by(Gene) %>%
  dplyr::group_modify(~ get_stats(.x, comp_use)) %>%
  dplyr::ungroup()

valid_x <- condition_levels
stats_all_draw <- stats_all %>%
  dplyr::filter(is.finite(p_adj), group1 %in% valid_x, group2 %in% valid_x) %>%
  dplyr::left_join(df_limits, by = "Gene") %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(
    n_br  = dplyr::n(),
    idx   = dplyr::row_number(),
    band0 = y_cap + band * pad_bot_frac,
    band1 = y_cap + band * (1 - pad_top_frac),
    # ensure non-collapsing band even if y_cap == 0
    band0 = ifelse(band1 <= band0, y_cap + 0.4 * pmax(band, 1), band0),
    band1 = ifelse(band1 <= band0, y_cap + 0.6 * pmax(band, 1), band1),
    frac  = ifelse(n_br == 1, 0.5, (idx - 1) / pmax(n_br - 1, 1)),
    y.position = band0 + frac * (band1 - band0)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-dplyr::any_of("step.increase"))  # avoid ggpubr name collision

# 5) Build the plot (jitter-only, 6 columns, text size 6)
p <- ggplot(df_sub, aes(x = Condition.of.Experiment, y = count,
                        color = Condition.of.Experiment)) +
  # facet-specific top limits (no x mapping here)
  geom_blank(data = df_limits, aes(y = y_upper), inherit.aes = FALSE) +
  
  geom_point(aes(shape = DonorID),
             position = position_jitter(width = 0.12, height = 0),
             size = 1, alpha = 0.9) +
  
  facet_wrap(~ Gene, scales = "free_y", ncol = 6, drop = TRUE) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    labels = scales::number_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  { if (nrow(stats_all_draw) > 0)
    ggpubr::stat_pvalue_manual(
      stats_all_draw,
      label = "p_label",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0,
      size = 2,
      step.increase = 0
    )
    else
      NULL
  } +
  theme_classic(base_size = 6) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line.x   = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6),
    plot.margin = margin(4, 6, 4, 6)
  ) +
  labs(y = "qPCR count")

p

### ONLY PLOT/ NO STATS
p <- ggplot(df_sub, aes(x = Condition.of.Experiment, y = count,
                        color = Condition.of.Experiment)) +
  # facet-specific upper limits
  geom_blank(data = df_limits, aes(y = y_upper), inherit.aes = FALSE) +
  
  # points
  geom_point(aes(shape = DonorID),
             position = position_jitter(width = 0.12, height = 0),
             size = 1, alpha = 0.9) +
  
  facet_wrap(~ Gene, scales = "free_y", ncol = 6) +
  scale_color_manual(values = condition_colors, drop = FALSE) +
  coord_trans(y = "log10", clip = "off") +
  scale_y_continuous(
    breaks = c(10,50,100,500,1000,5000,10000,50000,100000,500000,1000000),
    labels = scales::label_number(accuracy = 1, big.mark = ","),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  theme_classic(base_size = 6) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.line.x   = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6),
    plot.margin = margin(4, 6, 4, 6)
  ) +
  labs(y = "")

p

pdf("./plots2/qPCR_all_NoStats.pdf", width = 6, height = 3.5)
(p)& theme(axis.text.x =  element_blank(), axis.title.x = element_blank()) 
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

