## Hypothesis testing for differential expression - Script from paper

#rTg4510

# Read samples metadata
coldata_tg <- read_excel("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/paper_Supplementary_tables/mmc2.xlsx",
                         sheet = "rTg4510")
coldata_tg <- as.data.table(coldata_tg)

rownames(coldata_tg) <- coldata_tg$`Animal ID`
colnames(coldata_tg) <- c("Animal_ID","Genotype","Age_months","Sample_ID","RIN_score",
                          "No_raw_reads","No_reads_after_filtering","Percent_Q30","Mean_Q",
                          "Unique_mapped_STAR_percent","Multiple_mapped_STAR_percent","Mapped_STARpercent")
coldata_tg$Genotype[coldata_tg$Genotype == "TG"] <- "rTg4510" 
coldata_tg$Age_months <- as.factor(coldata_tg$Age_months)
coldata_tg$Genotype <- as.factor(coldata_tg$Genotype)
coldata_tg$Genotype <- relevel(coldata_tg$Genotype, "WT")
levels(coldata_tg$Genotype)
head(coldata_tg)

#read counts
supp <- fread("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/Data/GSE125957_processed_data.csv")

# Remove first five columns (chr, start, end, strand, length)
countdata_tg <- data.frame(supp[,7:65])  #66:127 for J20
rownames(countdata_tg)
colnames(countdata_tg)

rownames(countdata_tg) <- supp$V1

all(coldata_tg$Sample_ID == colnames(countdata_tg))

row.names(coldata_tg) <- colnames(countdata_tg)

# Convert to matrix
countdata_tg <- as.matrix(countdata_tg) # count matrix raw reads
head(countdata_tg)

# design: ~condition + time + condition:time (time = age; condition = genotype)
dds_tg <- DESeqDataSetFromMatrix(countData = countdata_tg, colData = coldata_tg, design = ~Genotype + Age_months + Genotype:Age_months)

# Pre-filtering the data set
nrow(dds_tg)
dds_tg <- dds_tg[ rowSums(counts(dds_tg)) > 6, ] # filtering for non-expressed and lowly expressed genes (sum of expression accross all samples at least 6 because final n = 6-8 mice per group)
nrow(dds_tg)

rld_tg = rlog(dds_tg)

##### Differential expression analysis

dds_Wald_tg <- DESeq(dds_tg, test="Wald")

# Building the results table
res_Wald_tg <- results(dds_Wald_tg)
results_Wald_frame_tg <- data.frame(res_Wald_tg, stringsAsFactors=FALSE)
stats_Wald_tg <- mcols(dds_Wald_tg)
rownames(stats_Wald_tg) <- rownames(dds_Wald_tg)
stats_Wald_frame_tg <- data.frame(stats_Wald_tg, stringsAsFactors=FALSE)

# PROGRESSIVE CHANGES (genes which show a diff expression profile over time across genotype): LRT TEST ###
dds_LRT_tg <- DESeq(dds_tg, reduced=~Genotype + Age_months, test="LRT") # DESeq: standard differential expression analysis steps are wrapped in this single function
# LRT: Likelihood-ratio test
# The Likelihood-Ratio test (sometimes called the likelihood-ratio chi-squared test) is a hypothesis test that helps you choose the “best” model between two nested models.
# “Nested models” means that one is a special case of the other.
resultsNames(dds_LRT_tg)

# Building the results table
res_LRT_tg <- results(dds_LRT_tg)
res_LRT_tg # pvalue = interaction
results_LRT_frame_tg <- data.frame(res_LRT_tg, stringsAsFactors=FALSE)
write.csv(results_LRT_frame_tg, file = "Tg4510_res_LRT_13may2021.csv")
stats_LRT_tg <- mcols(dds_LRT_tg)
head(stats_LRT_tg)
stats_LRT_frame_tg <- data.frame(stats_LRT_tg, stringsAsFactors=FALSE)

# EFFECT OF AGE: MANUAL LRT
full_tg <- stats::model.matrix.default(~Genotype*Age_months, data = as.data.frame(colData(dds_tg)))
head(full_tg)
reduced_tg <- full_tg[,-c(3:5)]
head(reduced_tg)

dds_age_tg <- DESeq(dds_tg, full=full_tg, reduced=reduced_tg, test="LRT")
dds_age_tg 
resultsNames(dds_age_tg)

# Building the results table
res_age_tg <- results(dds_age_tg)
results_age_frame_tg <- data.frame(res_age_tg, stringsAsFactors=FALSE)
write.csv(results_age_frame_tg, file = "Tg4510_res_age_13may2021.csv")
stats_age_tg <- mcols(dds_age_tg)
head(stats_age_tg)
stats_age_frame_tg <- data.frame(stats_age_tg, stringsAsFactors=FALSE)

# Stats Tables (final)
FDR_adj_genotype_tg <- p.adjust(stats_Wald_tg[,"WaldPvalue_Genotype_rTg4510_vs_WT"], method = "fdr")
FDR_adj_age_tg <- p.adjust(stats_age_tg[,"LRTPvalue"], method = "fdr")
FDR_adj_LRT_tg <- p.adjust(stats_LRT_tg[,"LRTPvalue"], method = "fdr")

stats_table_tg <-cbind(FDR_adj_genotype_tg, as.data.frame(stats_Wald_tg[,c("WaldPvalue_Genotype_rTg4510_vs_WT", "Genotype_rTg4510_vs_WT")]), 
                       FDR_adj_age_tg,
                       res_age_tg[,"pvalue"],	
                       as.data.frame(stats_Wald_tg[,c("WaldPvalue_Age_months_4_vs_2","Age_months_4_vs_2","WaldPvalue_Age_months_6_vs_2","Age_months_6_vs_2",
                                                      "WaldPvalue_Age_months_8_vs_2","Age_months_8_vs_2")]), 
                       FDR_adj_LRT_tg,
                       res_LRT_tg[,"pvalue"],  
                       as.data.frame(stats_Wald_tg[,c("WaldPvalue_GenotyperTg4510.Age_months4","GenotyperTg4510.Age_months4","WaldPvalue_GenotyperTg4510.Age_months6","GenotyperTg4510.Age_months6",
                                                      "WaldPvalue_GenotyperTg4510.Age_months8","GenotyperTg4510.Age_months8")]))
rownames(stats_table_tg) <- rownames(res_LRT_tg)
write.csv(stats_table_tg, file = "Tg4510_stats_table_13may2021.csv")

sig_genotype_tg <- stats_table_tg[which(stats_table_tg[,"FDR_adj_genotype_tg"]<0.05),]
write.csv(sig_genotype_tg, file = "Tg4510_sig_genotype_13may2021.csv")


######################################

#J20

# Read samples metadata
coldata_j20 <- read_excel("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/paper_Supplementary_tables/mmc2.xlsx",
                          sheet = "J20")
coldata_j20 <- as.data.table(coldata_j20)
rownames(coldata_j20) <- coldata_j20$`Animal ID`
colnames(coldata_j20) <- c("Animal_ID","Genotype","Age_months","Sample_ID","RIN_score",
                           "No_raw_reads","No_reads_after_filtering","Percent_Q30","Mean_Q",
                           "Unique_mapped_STAR_percent","Multiple_mapped_STAR_percent","Mapped_STARpercent")
coldata_j20$Genotype[coldata_j20$Genotype == "TG"] <- "J20" 
coldata_j20$Age_months <- as.factor(coldata_j20$Age_months)
coldata_j20$Genotype <- as.factor(coldata_j20$Genotype)
coldata_j20$Genotype <- relevel(coldata_j20$Genotype, "WT")
levels(coldata_j20$Genotype)
head(coldata_j20)

# Read counts
supp <- fread("/Users/eleniballi/Documents/Bioinformatics/MSc_Courses/Genomics_And_Bioinformatics_MSc_EPFL/Project/Analysis/Data/GSE125957_processed_data.csv")

# Remove first five columns (chr, start, end, strand, length)
countdata_j20 <- data.frame(supp[,66:127])
rownames(countdata_j20)
colnames(countdata_j20)

rownames(countdata_j20) <- supp$V1

all(coldata_j20$Sample_ID == colnames(countdata_j20))

row.names(coldata_j20) <- colnames(countdata_j20)

# Convert to matrix
countdata_j20 <- as.matrix(countdata_j20) # count matrix raw reads
head(countdata_j20)

# design: ~condition + time + condition:time (time = age; condition = genotype)
dds_j20 <- DESeqDataSetFromMatrix(countData = countdata_j20, colData = coldata_j20, design = ~Genotype + Age_months + Genotype:Age_months)

# Pre-filtering the data set
# Remove the rows that have no or nearly no information about the amount of gene expression
# Removing rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds_j20)
dds_j20 <- dds_j20[ rowSums(counts(dds_j20)) > 6, ] # filtering for non-expressed and lowly expressed genes (sum of expression accross all samples at least 6 because final n = 6-8 mice per group)
nrow(dds_j20)

rld_j20 = rlog(dds_j20)

##### Differential expression analysis

dds_Wald_j20 <- DESeq(dds_j20, test="Wald")

resultsNames(dds_Wald_j20)

# Building the results table
res_Wald_j20 <- results(dds_Wald_j20)
res_Wald_j20
results_Wald_j20_frame <- data.frame(res_Wald_j20, stringsAsFactors=FALSE)

stats_Wald_j20 <- mcols(dds_Wald_j20)
rownames(stats_Wald_j20) <- rownames(dds_Wald_j20)
head(stats_Wald_j20)
stats_Wald_20_frame <- data.frame(stats_Wald_j20, stringsAsFactors=FALSE)

# PROGRESSIVE CHANGES (genes which show a diff expression profile over time across genotype): LRT TEST ###
dds_LRT_j20 <- DESeq(dds_j20, reduced=~Genotype + Age_months, test="LRT") # DESeq: standard differential expression analysis steps are wrapped in this single function
# LRT: Likelihood-ratio test
# The Likelihood-Ratio test (sometimes called the likelihood-ratio chi-squared test) is a hypothesis test that helps you choose the “best” model between two nested models.
# “Nested models” means that one is a special case of the other.

dds_LRT_j20 
resultsNames(dds_LRT_j20)

# Building the results table
res_LRT_j20 <- results(dds_LRT_j20)
res_LRT_j20 # pvalue = interaction
results_LRT_j20_frame <- data.frame(res_LRT_j20, stringsAsFactors=FALSE)
write.csv(results_LRT_j20_frame, file = "J20_res_LRT_24may2021.csv")
stats_LRT_j20 <- mcols(dds_LRT_j20)
head(stats_LRT_j20)
stats_LRT_j20_frame <- data.frame(stats_LRT_j20, stringsAsFactors=FALSE)

# EFFECT OF AGE: MANUAL LRT #
full_j20 <- stats::model.matrix.default(~Genotype*Age_months, data = as.data.frame(colData(dds_j20)))
head(full_j20)
reduced_j20 <- full_j20[,-c(3:5)]
head(reduced_j20)
dds_age_j20 <- DESeq(dds_j20, full=full_j20, reduced=reduced_j20, test="LRT")
dds_age_j20 
resultsNames(dds_age_j20)

# Building the results table
res_age_j20 <- results(dds_age_j20)
res_age_j20 
results_age_j20_frame <- data.frame(res_age_j20, stringsAsFactors=FALSE)
write.csv(results_age_j20_frame, file = "J20_res_age_24may2021.csv")
stats_age_j20 <- mcols(dds_age_j20)
head(stats_age_j20)
stats_age_j20_frame <- data.frame(stats_age_j20, stringsAsFactors=FALSE)

# Stats Tables (final)
FDR_adj_genotype_j20 <- p.adjust(stats_Wald_j20[,"WaldPvalue_Genotype_J20_vs_WT"], method = "fdr") # we have calculated the FDR-adj pvalue ourselves because 1) they were not available for all analysis,and 2) from DESeq some were "NA"
FDR_adj_age_j20 <- p.adjust(stats_Wald_j20[,"WaldPvalue_Genotype_J20_vs_WT"], method = "fdr")
FDR_adj_age_j20 <- p.adjust(stats_age_j20[,"LRTPvalue"], method = "fdr")
FDR_adj_LRT_j20 <- p.adjust(stats_LRT_j20[,"LRTPvalue"], method = "fdr")

stats_table_j20 <-cbind(FDR_adj_genotype_j20, as.data.frame(stats_Wald_j20[,c("WaldPvalue_Genotype_J20_vs_WT", "Genotype_J20_vs_WT")]), 
                        FDR_adj_age_j20,
                        res_age_j20[,"pvalue"],	
                        as.data.frame(stats_Wald_j20[,c("WaldPvalue_Age_months_8_vs_6","Age_months_8_vs_6","WaldPvalue_Age_months_10_vs_6","Age_months_10_vs_6",
                                                        "WaldPvalue_Age_months_12_vs_6","Age_months_12_vs_6")]), 
                        FDR_adj_LRT_j20,
                        res_LRT_j20[,"pvalue"],  
                        as.data.frame(stats_Wald_j20[,c("WaldPvalue_GenotypeJ20.Age_months8","GenotypeJ20.Age_months8","WaldPvalue_GenotypeJ20.Age_months10","GenotypeJ20.Age_months10",
                                                        "WaldPvalue_GenotypeJ20.Age_months12","GenotypeJ20.Age_months12")]))
rownames(stats_table_j20) <- rownames(res_LRT_j20)
write.csv(stats_table_j20, file = "J20_stats_table_24may2021.csv")

sig_genotype_j20 <- stats_table_j20[which(stats_table_j20[,"FDR_adj_genotype_j20"]<0.05),]
write.csv(sig_genotype_j20, file = "J20_sig_genotype_24may2021.csv")

