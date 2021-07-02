## 3.a. Heatmaps of differentially expressed genes
# Remove first five columns (chr, start, end, strand, length)
countdata_tg <- data.frame(supp[,7:65])  #66:127 for J20
rownames(countdata_tg)
colnames(countdata_tg)

rownames(countdata_tg) <- supp$V1

all(samples_tg$Sample_ID == colnames(countdata_tg))

row.names(samples_tg) <- colnames(countdata_tg)

# Convert to matrix
countdata_tg <- as.matrix(countdata_tg) # count matrix raw reads
head(countdata_tg)

countdata_tg
dds_sizefactors_tg <- estimateSizeFactors(dds_tg)
normalized_counts_tg <- counts(dds_sizefactors_tg, normalized=TRUE)
head(normalized_counts_tg)

countdata_j20
dds_sizefactors_j20<- estimateSizeFactors(dds_j20)
normalized_counts_j20 <- counts(dds_sizefactors_j20, normalized=TRUE)
head(normalized_counts_j20)

tg_diff_expr <- read.table("Tg4510_sig_genotype_13may2021.csv", sep=",", header=TRUE, row.names=1)
j20_diff_expr <- read.table("J20_sig_genotype_13may2021.csv", sep=",", header=TRUE, row.names=1)
tg_diff <- row.names(tg_diff_expr)
j20_diff <- row.names(j20_diff_expr)
tg <- as.data.frame (tg_diff)
j20 <- as.data.frame (j20_diff)
names(tg)[1] <- "genes"
tg$genes
names(j20)[1] <- "genes"
j20$genes

j <- which(rownames(rld_tg) %in%  tg$genes)
j
k <- which(rownames(rld_j20) %in%  j20$genes)
k

norm_counts_subset_tg<- as.matrix(normalized_counts_tg[j,])
norm_counts_subset_j20<- as.matrix(normalized_counts_j20[k,])

data_tg <- norm_counts_subset_tg[,order(samples_tg$Genotype, samples_tg$Age_months)]
data_j20 <- norm_counts_subset_j20[,order(samples_j20$Genotype, samples_j20$Age_months)]

symmetric_breaks <- seq(-max(abs(data)), max(abs(data)), length.out=101)

ann_colors = list(Genotype=c(rTg4510="#004680", WT_rTg4510="#faf0be", J20="#6d134c", WT_J20="#ff033e"),
                  Age_months=c("2"="#dbead5", "4"="#b7d5ac", "6"="#93bf85", "8"="#6eaa5e", "10"="#469536", "12"="#008000"))

pheatmap(
  data_tg,
  color = colorRampPalette(c("#0000FF", "white", "#FF0000"))(100), 
  breaks = seq(-max(abs(data)), max(abs(data)), length.out=101),
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  clustering_distance_rows="correlation",
  clustering_distance_cols="correlation",
  clustering_method="ward.D2",
  annotation_col=samples_tg_j20[order(samples_tg_j20$Genotype,samples_tg_j20$Age_months),c("Genotype", "Age_months")],
  treeheight_row=10,
  treeheight_column=0.5,
  show_colnames=FALSE,
  fontsize_row = 3,
  fontsize = 7, 
  scale="row",
)

pheatmap(
  data_j20,
  color = colorRampPalette(c("#0000FF", "white", "#FF0000"))(100), 
  breaks = seq(-max(abs(data)), max(abs(data)), length.out=101),
  border_color = NA,
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  clustering_distance_rows="correlation",
  clustering_distance_cols="correlation",
  clustering_method="ward.D2",
  annotation_col=samples_tg_j20[order(samples_tg_j20$Genotype,samples_tg_j20$Age_months),c("Genotype", "Age_months")],
  treeheight_row=15,
  treeheight_column=0.5,
  show_colnames=FALSE,
  fontsize_row = 6,
  scale="row" 
)

## 3.b. PCA plots - differentially expressed genes
tg_diff_expr <- read.table("Tg4510_sig_genotype_13may2021.csv", sep=",", header=TRUE, row.names=1)
row.names(tg_diff_expr)

dif <- which(rownames(rld_table_tg) %in%  row.names(tg_diff_expr))
dif

j20_diff_expr <- read.table("J20_sig_genotype_13may2021.csv", sep=",", header=TRUE, row.names=1)
row.names(j20_diff_expr)

dif_j20 <- which(rownames(rld_table_j20) %in%  row.names(j20_diff_expr))
dif_j20


plot1 = ggplot(pcaData_tg_dif, aes(PC1, PC2, color=Genotype, shape=Age_months)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(legend.key.size = unit(0.5, 'cm')) +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(values = c("#52854C", "#D16103")) +
  coord_fixed()
plot2 = ggplot(pcaData_j20_dif, aes(PC1, PC2, color=Genotype, shape=Age_months)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(legend.key.size = unit(0.5, 'cm')) +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  coord_fixed()
grid.arrange(plot1, plot2, ncol=2)

