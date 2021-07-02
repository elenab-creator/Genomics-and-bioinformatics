## 4. Functional analysis, overrepresentation analysis - rTg4510
# 
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
unique(ah$species) %>% View()
ens <- query(ah, c("Mus musculus", "EnsDb"))
ens <- ens[["AH53222"]]
genes(ens, return.type = "data.frame") %>% View()
annotations_ahb <- genes(ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>% 
  dplyr::filter(gene_name %in% rownames(res_Wald_tg))
class(annotations_ahb$entrezid)
which(map(annotations_ahb$entrezid, length) > 1)
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()
which(is.na(annotations_ahb$gene_name)) %>% length()
which(duplicated(annotations_ahb$gene_name)) %>% length()
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)
annotations_ahb %>% nrow()
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]
annotations_ahb %>% nrow()
which(is.na(annotations_ahb$entrezid)) %>%  length()
res_Wald_tg_frame <- as.data.frame(res_Wald_tg)

# Add the rownames (gene names) as  column
res_Wald_tg_frame_genes <- cbind (gene_name = rownames(res_Wald_tg_frame), res_Wald_tg_frame)
head(res_Wald_tg_frame_genes)
rownames(res_Wald_tg_frame_genes) <- NULL
head(res_Wald_tg_frame_genes)

# Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_Wald_tg_frame_genes, annotations_ahb)
head(res_ids)
names(res_ids)[names(res_ids) == "gene_id"] <- "ENSEMBL"
head(res_ids)

# Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(res_ids$entrezid)

# Extract significant results
sig <- dplyr::filter(res_ids, padj < 0.05)
sig_genes <- as.character(sig$entrezid)

# Run GO enrichment analysis 
ego <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                # keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "clusterProfiler_Tg.csv")
save(ego, file="ego_Tg.rda")

# Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)

# Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, cex_line = 0.2, showCategory = 50)

