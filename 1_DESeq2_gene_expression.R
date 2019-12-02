library(readr)
library(DESeq2)
library(ggplot2)
library(rtracklayer)

#Takashima rna-seq data ribo-zero directional ----------------------------------

raw_counts <- read_delim("1_DESeq2_gene_expression/raw_read_counts_takashima_GRCh38.87_anno_opposing_strand_mrna.txt",
                         "\t", escape_double = FALSE, col_types = cols(Chromosome = col_character()),
                         trim_ws = TRUE) #done on mrna

#Rename duplicated probe names (should not remove anything)---------------------

raw_counts$ID <-  make.names(raw_counts$ID, unique=TRUE)
raw_counts <- as.data.frame(raw_counts) #data.frame from tibble
rownames(raw_counts) <- make.names(raw_counts$ID, unique=TRUE)
just_raw_counts <- raw_counts[,14:19]

#Run DESeq ---------------------------------------------------------------------

column_data <- data.frame(cell_type=as.factor(c("Naive","Naive","Naive",
                                                "Primed", "Primed", "Primed"))) 

# Make a DESeq data set from the counts and the design and specify which factors in
# the design to test
count_data_set <- DESeqDataSetFromMatrix(countData=just_raw_counts, colData=column_data, 
                                         design= ~ cell_type)

# Perform the analysis (pairwise comparison)
count_data_set <- DESeq(count_data_set)

# Retrieve the full set of results. Do independent filtering by default.
binomial_result <- results(count_data_set) #,independentFiltering=FALSE)

# Remove unmeasured results
binomial_result <- na.omit(binomial_result) 

de_results <- as.data.frame(binomial_result)



#convert ensembl ID to feature--------------------------------------------------

gtf <- readGFF("1_DESeq2_gene_expression/Homo_sapiens.GRCh38.87.gtf.gz", version=2L, 
               tags = c("gene_id", "gene_name", "gene_biotype"))

gtf_unq <- gtf[!duplicated(gtf$gene_id),]

gene_ids <- strsplit(rownames(de_results), "\\.")
gene_id <- plyr::ldply(gene_ids, rbind)

de_results$gene_name <- ifelse(is.na(gene_id$`2`), gtf_unq$gene_name[match(gene_id$`1`, gtf_unq$gene_id)],
                               paste(gtf_unq$gene_name[match(gene_id$`1`, gtf_unq$gene_id)], gene_id$`2`, sep="."))
de_results$gene_biotype <- gtf_unq$gene_biotype[match(gene_id$`1`, gtf_unq$gene_id)]



#merge with seqmonk log2 length corrected counts -------------------------------

fpkm_seqmonk<- read_delim("1_DESeq2_gene_expression/normalised_genes_takashima_GRCh38.87_anno_opposing_strand_fpkm.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols(Chromosome = col_character()))


fpkm_seqmonk$naive_mean <- rowMeans(fpkm_seqmonk[c("ERR590398_H9_reset_R2_1_val_1_GRCh38_hisat2.bam", 
                                                   "ERR590399_H9_reset_R3_1_val_1_GRCh38_hisat2.bam", 
                                                   "ERR590400_H9_reset_R1_1_val_1_GRCh38_hisat2.bam")], na.rm = TRUE)

fpkm_seqmonk$primed_mean <- rowMeans(fpkm_seqmonk[c("ERR590401_H9_R3_PE1_val_1_GRCh38_hisat2.bam",
                                                    "ERR590408_H9_R1_PE1_val_1_GRCh38_hisat2.bam",
                                                    "ERR590410_H9_R2_PE1_val_1_GRCh38_hisat2.bam")], na.rm = TRUE)


raw_counts_sub <- raw_counts[,c("Chromosome","Start","End", "ID")]
names(raw_counts_sub)[4] <- "ID_unq"
raw_counts_sub$ID <- as.character(plyr::ldply(strsplit(raw_counts_sub$ID, "\\."), rbind)[,1])
fpkm_seqmonk_id <- merge(fpkm_seqmonk[,c("Chromosome","Start","End", "naive_mean", "primed_mean", "ID")], 
                         raw_counts_sub, 
                         by=c("Chromosome","Start","End", "ID"), all.x=TRUE)



fpkm_mean_df <- data.frame(Naive_mean_fpkm_log2=fpkm_seqmonk_id$naive_mean, Primed_mean_fpkm_log2=fpkm_seqmonk_id$primed_mean)
rownames(fpkm_mean_df) <- fpkm_seqmonk_id$ID_unq


de_results_fpkm <- merge(de_results, fpkm_mean_df, by="row.names", all.x=TRUE)
#r doesn't like -
de_results_fpkm$gene_name <- gsub("-", ".", de_results_fpkm$gene_name)

#Get protain coding genes only
de_results_fpkm_pc <- de_results_fpkm[de_results_fpkm$gene_biotype %in% "protein_coding",]
de_results_fpkm_lincRNA <- de_results_fpkm[de_results_fpkm$gene_biotype %in% "lincRNA",]


# Write the hit names to a file
write.table(de_results_fpkm, file="1_DESeq2_gene_expression/de_genes_takashima_GRCh38.87_anno_opposing_strand_all.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(de_results_fpkm_pc, file="1_DESeq2_gene_expression/de_genes_takashima_GRCh38.87_anno_opposing_strand_prot_genes.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(de_results_fpkm_lincRNA, file="1_DESeq2_gene_expression/de_genes_takashima_GRCh38.87_anno_opposing_strand_lincRNA.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

