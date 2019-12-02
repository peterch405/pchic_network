
#Release 25 (GRCh38.p7)

library(rtracklayer)
library(GenomicRanges)
library(data.table)


source("annotate_promoters_functions.R")

#ensembl annotation
#ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/

gtf <- rtracklayer::import("annotate_promoters/Homo_sapiens.GRCh38.87.gtf.gz", 
                           colnames = c("gene_id", "gene_name", "gene_biotype",
                                        "transcript_id", "transcript_biotype", 
                                        "transcript_name", "type"), format="gff2")


gtf_transcript_all <- gtf[gtf$type == "transcript",] #annotate with transcript instead of gene to get TSS and hence promoters
gtf_gene_all <- gtf[gtf$type == "gene",]

#only keep certain biotypes
allowed_biotypes <- c("protein_coding", "non_coding", "antisense", "snRNA", "miRNA", "snoRNA", "lincRNA")
gtf_transcript_allowed <- gtf_transcript_all[gtf_transcript_all$transcript_biotype %in% allowed_biotypes]
gtf_transcript_allowed <- gtf_transcript_allowed[gtf_transcript_allowed$gene_biotype %in% allowed_biotypes]

gtf_gene_allowed <- gtf_gene_all[gtf_gene_all$gene_biotype %in% allowed_biotypes]


#add chr
gtf_transcript_allowed <- diffloop::addchr(gtf_transcript_allowed)
gtf_gene_allowed <- diffloop::addchr(gtf_gene_allowed)

baited_g38_promoter_HindIII_fragments_gr <- import("annotate_promoters/baited_g38_promoter_HindIII_fragments.txt",
                                                   format = "bed")
baited_g38_promoter_HindIII_fragments_gr <- diffloop::addchr(baited_g38_promoter_HindIII_fragments_gr)


#Need to create 1000bp upstream of gene and overlap with baited promoters
gtf_transcript <- promoter_start(gtf_transcript_allowed)
gtf_gene <- promoter_start(gtf_gene_allowed)

gtf_gene_prom <- gtf_gene
gtf_transcript_prom <- gtf_transcript
ranges(gtf_gene_prom)<- IRanges(gtf_gene$prom_start, gtf_gene$prom_end)
ranges(gtf_transcript_prom)<- IRanges(gtf_transcript$prom_start, gtf_transcript$prom_end)


#Annotate baits ----------------------------------------------------------------

# gtf_gr <- gtf_transcript
# fragments_gr <- baited_g38_promoter_HindIII_fragments_gr

bait_annotation_transcript_all <- annotate_fragments(gtf_transcript_prom, 
                                                     baited_g38_promoter_HindIII_fragments_gr, 
                                                     annotate_with="transcript")

bait_annotation_gene_name_all <- annotate_fragments(gtf_transcript_prom, 
                                               baited_g38_promoter_HindIII_fragments_gr, 
                                               annotate_with="gene")

bait_annotation_gene_all <- annotate_fragments(gtf_gene_prom, 
                                               baited_g38_promoter_HindIII_fragments_gr, 
                                               annotate_with="gene")


#1603 unnannotated baits!

write.table(bait_annotation_transcript_all, "annotate_promoters/baited_g38_promoter_HindIII_fragments_transcript_anno_20181207.txt",
            quote = FALSE, row.names = FALSE)
write.table(bait_annotation_gene_name_all, "annotate_promoters/baited_g38_promoter_HindIII_fragments_transcript_gene_anno_20181207.txt",
            quote = FALSE, row.names = FALSE)
write.table(bait_annotation_gene_all, "annotate_promoters/baited_g38_promoter_HindIII_fragments_gene_anno_20181207.txt",
            quote = FALSE, row.names = FALSE)

saveRDS(bait_annotation_transcript_all, "annotate_promoters/baited_g38_promoter_HindIII_fragments_transcript_anno_20181207.rds")
saveRDS(bait_annotation_gene_name_all, "annotate_promoters/baited_g38_promoter_HindIII_fragments_transcript_gene_anno_20181207.rds")
saveRDS(bait_annotation_gene_all, "annotate_promoters/baited_g38_promoter_HindIII_fragments_gene_anno_20181207.rds")

################################################################################
#######Annotate all restriction enzyme fragments with possible promoters########
################################################################################

library(readr)
DHs_GRCh38_HindIII <- read_delim("annotate_promoters/Digest_Homo_sapiens_GRCh38_HindIII_None_14-43-31_10-02-2016.txt.gz", 
                                                                              "\t", escape_double = FALSE, col_types = cols(Chromosome = col_character()), 
                                                                              trim_ws = TRUE, skip = 1)
#add chr to chroms
DHs_GRCh38_HindIII$Chromosome <- sub("^", "chr", DHs_GRCh38_HindIII$Chromosome)


names(DHs_GRCh38_HindIII)[1:3] <- c("chrom", "start", "end")

gr_DHs_GRCh38_HindIII <- makeGRangesFromDataFrame(DHs_GRCh38_HindIII, keep.extra.columns = TRUE)

#annotate with transcript gene names
annotation_gene_name_all <- annotate_fragments(gtf_transcript_prom, gr_DHs_GRCh38_HindIII, 
                                               annotate_with = "gene", 
                                               highlight_baited = baited_g38_promoter_HindIII_fragments_gr,
                                               keep_meta = TRUE, add_id=TRUE)


write.table(annotation_gene_name_all, "annotate_promoters/Digest_Homo_sapiens_GRCh38_HindIII_None_14-43-31_10-02-2016_anno_20181209.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

#annotate with gene's
annotation_gene_all <- annotate_fragments(gtf_gene_prom, gr_DHs_GRCh38_HindIII, 
                                          annotate_with = "gene", 
                                          highlight_baited = baited_g38_promoter_HindIII_fragments_gr,
                                          keep_meta = TRUE, add_id=TRUE)

write.table(annotation_gene_name_all, "annotate_promoters/Digest_Homo_sapiens_GRCh38_HindIII_None_14-43-31_10-02-2016_anno.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

