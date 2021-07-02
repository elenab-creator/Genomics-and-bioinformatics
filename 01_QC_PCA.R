## 1. QC Samples - PCA all samples
# Read sample metadata
samples_tg <- read_excel("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/Paper_Supplementary_tables/mmc2.xlsx",
                         sheet = "rTg4510")

samples_j20 <- read_excel("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/Paper_Supplementary_tables/mmc2.xlsx",
                          sheet = "J20")
samples_tg <- as.data.table(samples_tg)
colnames(samples_tg) <- c("Animal_ID","Genotype","Age_months","Sample_ID","RIN_score",
                          "No_raw_reads","No_reads_after_filtering","Percent_Q30","Mean_Q",
                          "Unique_mapped_STAR_percent","Multiple_mapped_STAR_percent","Mapped_STARpercent")

samples_j20 <- as.data.table(samples_j20)
colnames(samples_j20) <- c("Animal_ID","Genotype","Age_months","Sample_ID","RIN_score",
                           "No_raw_reads","No_reads_after_filtering","Percent_Q30","Mean_Q",
                           "Unique_mapped_STAR_percent","Multiple_mapped_STAR_percent","Mapped_STARpercent")

samples_tg$Genotype[samples_tg$Genotype == "TG"] <- "rTg4510" 
samples_j20$Genotype[samples_j20$Genotype == "TG"] <- "J20" 
samples_tg$Genotype[samples_tg$Genotype == "WT"] <- "WT_rTg4510" 
samples_j20$Genotype[samples_j20$Genotype == "WT"] <- "WT_J20" 

all(colnames(samples_tg) == colnames(samples_j20))
unique(samples_tg$Genotype) 
unique(samples_j20$Genotype)

samples_tg_j20 <- rbind(samples_tg,samples_j20)

samples_tg_j20$Age_months <- as.factor(samples_tg_j20$Age_months)
samples_tg_j20$Genotype <- as.factor(samples_tg_j20$Genotype)
samples_tg_j20 <- as.data.frame(samples_tg_j20)

# Read raw counts
supp <- fread("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/Data/GSE125957_processed_data.csv")
countdata <- data.frame(supp[,7:127])
rownames(countdata)
colnames(countdata)
rownames(countdata) <- supp$V1

all(samples_tg_j20$Sample_ID == colnames(countdata))

row.names(samples_tg_j20) <- colnames(countdata)

countdata <- as.matrix(countdata)
head(countdata)

# Create data object for DESeq2
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = samples_tg_j20, design = ~1)

# Gene filtering
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ] # filtering for non-expressed and lowly expressed genes (sum of expression accross all samples at least 6 because final n = 6-8 mice per group)
nrow(dds)

# rlog transformation
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
rld_data <- assay(rld)

# PCA plot
pcaData = plotPCA(rld, intgroup = c("Genotype", "Age_months"), returnData=TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Genotype, shape=Age_months)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  # ggtitle("rTg4510 and J20 datasets combined") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  scale_color_manual(values = c("#E7B800", "#D16103", "#00AFBB", "#52854C")) +
  coord_fixed()