library(data.table)
library(igraph)
library(GenomicRanges)
library(readr)

source("annotate_promoters_functions.R")
source("network_enhancer_class.R")
source("network_enhancers_functions.R")


#' Create column of IDs matching overlaps between two gr's
#' essentailly giving an ID to intervals with multiple overlaps
#' 
#' @param gr GenomicRanges 1
#' @param gr2 GenomicRanges 2
#' 
#' gr2[1]
#' gr[queryHits(ovlp)[which(subjectHits(ovlp) == "1")]]
#' 
#' test <- aggregate_hits(Naive_enh, Primed_enh)
#' 
#' test$gr2[test$gr2$overlap %in% "q2"]
#' test$gr[test$gr$overlap %in% "q2"]
#' 
#' test$gr2[test$gr2$overlap %in% "s1"]
#' test$gr[test$gr$overlap %in% "s1"]
#' 
#' @return list with gr and gr2 with an overlap column that corresponds to
#' intervals that overlap
#' 
aggregate_hits <- function(gr, gr2){
  ovlp <- findOverlaps(gr, gr2)
  multi_query <- names(table(queryHits(ovlp))[table(queryHits(ovlp)) > 1]) 
  multi_q_pos <- queryHits(ovlp) %in% multi_query
  multi_subject <- names(table(subjectHits(ovlp))[table(subjectHits(ovlp)) > 1]) 
  multi_s_pos <- subjectHits(ovlp) %in% multi_subject
  
  stopifnot(!any(multi_s_pos & multi_q_pos))
  
  gr$overlap <- NA
  gr$overlap[queryHits(ovlp)] <- paste("q", queryHits(ovlp), sep = "")
  gr$overlap[queryHits(ovlp)[multi_s_pos]] <- paste("s", subjectHits(ovlp)[multi_s_pos], sep = "")
  
  gr2$overlap <- NA
  gr2$overlap[subjectHits(ovlp)] <- paste("q", queryHits(ovlp), sep = "")
  gr2$overlap[subjectHits(ovlp)[multi_s_pos]] <- paste("s", subjectHits(ovlp)[multi_s_pos], sep = "")
  
  
  return(list(gr=gr, gr2=gr2))
}



#' Make hindiii lookup table with center coordinates
#' @param df_to_add data.frame to convert to mid coordinates
#' @param start vector with start coordinates
#' @param end vector with end coordinates
center_coord <- function(df_to_add, start, end){
  
  stopifnot(length(start)==length(end))
  start <- as.numeric(as.character(start))
  end <- as.numeric(as.character(end))
  mid_start <- vector("numeric", length(start))
  mid_end <- vector("numeric", length(end))
  for(i in 1:length(start)){
    mid <- start[i] + (end[i]-start[i])/2
    if(mid%%1==0){
      mid_start[i] <- mid
      mid_end[i] <- mid +1
    }else{
      mid_start[i] <- floor(mid)
      mid_end[i] <- ceiling(mid)
    }
  }
  df_to_add$start <- mid_start
  df_to_add$end <- mid_end
  return(df_to_add)
}




#' Find the overlap of two genomic ranges
#' 
#' https://support.bioconductor.org/p/56880/
#' 
#' @param gr GenomicRanges to which overlap will be added
#' @param gr2 Second GenomicRanges with which overlap will be calculated
#' 
#' @return gr with overlap column
find_overlap <- function(gr, gr2){
  cover2.gr <- reduce(gr2)
  hits <- findOverlaps(gr,cover2.gr)
  gr.over <-pintersect(gr[queryHits(hits)],cover2.gr[subjectHits(hits)])
  gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x)
    sum(width(x)))
  gr$overlap<- 0
  gr$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
  
  return(gr)
}


#Load non-trans network --------------------------------------------------------
net_all_s <- readRDS("2_network_make/net_all_s_20200911.rds")



#HindIII fragment coordinates lookup genomicsranges -------------------------------------

#from Digest_Homo_sapiens_GRCh38_HindIII_None_14-43-31_10-02-2016_anno.txt in
#in network_visulisation.py
#should be 1-based coordinates
hindiii_lookup <- read_delim("/mnt/Projects/CHiC_naive_primed/hindiii_lookup.tab", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
names(hindiii_lookup) <- c("chrom", "start", "end", "ID")

hindiii_gr <- makeGRangesFromDataFrame(hindiii_lookup, keep.extra.columns = TRUE)

saveRDS(hindiii_gr, "lookup_files/hindiii_gr.rds")

#weird behaviour, if you convert gr straigt to as.data.table, it messes up the original table
hindiii_coord <- data.table(as(hindiii_gr,"data.frame"))
hindiii_coord[, c("width","strand"):=NULL ]
setkey(hindiii_coord, ID)

saveRDS(hindiii_coord, "lookup_files/hindiii_coord.rds")


hindiii_mid_lookup <- center_coord(hindiii_lookup, hindiii_lookup$start, hindiii_lookup$end)

saveRDS(hindiii_mid_lookup, "lookup_files/hindiii_mid_lookup.rds")





#Gene expression lookup table --------------------------------------------------
#annotate nodes with gene expression
links_nodes_cat_col_coord_deb2b <- readRDS("2_network_make/links_nodes_cat_col_coord_deb2b_20200911.rds")
nodes_all <- links_nodes_cat_col_coord_deb2b$nodes
de_genes <- readr::read_delim("/mnt/Projects/CHiC_naive_primed/RNA-seq/de_genes_takashima_GRCh38.87_anno_opposing_strand_prot_genes.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
names(de_genes)[8] <- "prot_genes"

#make gene transcript lookup table
nodes_sub <- nodes_all[,c("ID", "prot_genes", "baited")]
nodes_all_genes <- tidyr::separate_rows(nodes_sub, prot_genes, sep=",", convert = TRUE)
nodes_all_prot <- nodes_all_genes[!is.na(nodes_all_genes$prot_genes),]

#some are not in RNA-seq data (these are left as NA)
#don't care about nodes without genes, 
nodes_all_prot <- merge(nodes_all_prot, de_genes, by="prot_genes")#, all.x=TRUE)

#remove unbaited aswell as you want to annotate the baited genes with expression not the enhancers
nodes_all_prot <- nodes_all_prot[nodes_all_prot$baited == 1,]
#make into data.table for fast lookup
setDT(nodes_all_prot, key="ID")

saveRDS(nodes_all_prot, "lookup_files/nodes_all_prot.rds")

# Get gene coordinates from gtf ------------------------------------------------

gtf <- readGFF("/mnt/Projects/CHiC_naive_primed/Homo_sapiens.GRCh38.87.gtf.gz", version=2L, 
               tags = c("gene_id", "gene_name", "gene_biotype"))

# gtf_unq <- gtf[!duplicated(gtf$gene_id),]
gtf_unq <- gtf[!duplicated(gtf$gene_id),]

#GTF gene to granges object
#for karyoplotR gene plotting
gtf_unq_chr <- gtf_unq[gtf_unq$seqid %in% c(1:22, "X"),]
gtf_unq_chr$seqid <- paste("chr", gtf_unq_chr$seqid, sep="")

gtf_unq_genes_gr <- makeGRangesFromDataFrame(gtf_unq_chr[gtf_unq_chr$gene_biotype == "protein_coding",], keep.extra.columns = TRUE)
saveRDS(gtf_unq_genes_gr, "lookup_files/gtf_unq_genes_gr.rds")



#Gene level annotation (not transcript level) for enhancer expression analysis----
#1-based coordinates, will need to convert to 0-based
bait_annotation_gene_all <- readRDS("/mnt/Projects/CHiC_naive_primed/baited_g38_promoter_HindIII_fragments_gene_anno_20181207.rds")

#only interested in protein coding
bait_annotation_gene_all$prot_genes <- filter_biotype(bait_annotation_gene_all, 
                                                      "protein_coding")
#add hindiii ID
bait_annotation_gene_all$start <- bait_annotation_gene_all$start-1

baited_genes_hindiii <- merge(bait_annotation_gene_all, hindiii_lookup, 
                              by=c("chrom", "start", "end"), all.x=TRUE)



baited_genes_hindiii_split <- tidyr::separate_rows(baited_genes_hindiii[c("prot_genes", "ID")], 
                                                   prot_genes, sep=",", convert = TRUE)

baited_prot_genes <- baited_genes_hindiii_split[!is.na(baited_genes_hindiii_split$prot_genes),]

baited_prot_genes <- merge(baited_prot_genes, de_genes, by="prot_genes")#, all.x=TRUE)

#make into data.table for fast lookup
setDT(baited_prot_genes, key="ID")

saveRDS(baited_prot_genes, "lookup_files/baited_prot_genes.rds")






#Network simplified lookup table -----------------------------------------------

#Set up search data.table from simplified network
net_data_all <- igraph::as_long_data_frame(net_all_s)
setDT(net_data_all, key = c("from_name", "to_name"))
#reduce to only needed columns
cols <- c("from_name", "to_name", "origin", "b2b", "score_naive", "score_primed")
net_data <- net_data_all[, ..cols]

saveRDS(net_data_all, "lookup_files/net_data_all.rds")
saveRDS(net_data, "lookup_files/net_data.rds")


#Chromhmm lookup table ---------------------------------------------------------

naive_nodes <- read_delim("0_network_data_preparation/naive_nodes_wchip_20200911.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

primed_nodes <- read_delim("0_network_data_preparation/primed_nodes_wchip_20200911.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

chromhmm_lookup <- merge(naive_nodes[c("ID","chromhmm")], primed_nodes[c("ID","chromhmm")], by="ID", suffixes = c("_naive", "_primed"))

saveRDS(chromhmm_lookup, "lookup_files/chromhmm_lookup.rds")



chromhmm_coord <- merge(chromhmm_lookup, hindiii_coord, by="ID")
naive_chromhmm <- chromhmm_coord[c("seqnames", "start", "end", "chromhmm_naive", "ID")]
primed_chromhmm <- chromhmm_coord[c("seqnames", "start", "end", "chromhmm_primed", "ID")]

naive_chromhmm_gr <- makeGRangesFromDataFrame(naive_chromhmm, keep.extra.columns = TRUE)
primed_chromhmm_gr <- makeGRangesFromDataFrame(primed_chromhmm, keep.extra.columns = TRUE)





#Gene annotation for circos plots ----------------------------------------------

#GFF is 1-based coordinate system 
gtf <- rtracklayer::import("/mnt/Projects/CHiC_naive_primed/Homo_sapiens.GRCh38.87.gtf.gz", 
               colnames = c("gene_id", "transcript_id", "gene_name", "gene_biotype", "type"), format="gff2")

gtf_gene_all <- gtf[gtf$type == "gene",] 

#remove extra chromosomes and add 'chr'

gtf_gene <- gtf_gene_all[seqnames(gtf_gene_all) %in% c(seq(1,22,1), "X")]
seqlevels(gtf_gene) <- c(seq(1,22,1), "X")
#add chr to seqnames
gtf_gene_chr <- diffloop::addchr(gtf_gene)

#filter out unintresting gene_biotypes
gtf_gene_gr <- gtf_gene_chr[gtf_gene_chr$gene_biotype %in% c("protein_coding", "lincRNA", 
                                              "snoRNA", "antisense", "miRNA", 
                                              "snRNA", "processed transcript")]


#hg38 chromosome sizes lookup --------------------------------------------------

hg38_chrom <- readr::read_delim("/mnt/Projects/CHiC_naive_primed/scripts/kentUtils/hg38.chrom.sizes", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
names(hg38_chrom) <- c("chrom", "size")
setDT(hg38_chrom, key="chrom")

saveRDS(hg38_chrom, "lookup_files/hg38_chrom.rds")

# Enhancer lookup tables -------------------------------------------------------

# Load ROSE bed -----------------------------------------------------------------
Naive_enh <- rtracklayer::import("/mnt/Projects/CHiC_naive_primed/ROSE/Naive_1500/Naive_H3K27ac_peaks_narrowPeak_Gateway_Enhancers.bed",
                                 format="bed")
Naive_suprenh <- rtracklayer::import("/mnt/Projects/CHiC_naive_primed/ROSE/Naive_1500/Naive_H3K27ac_peaks_narrowPeak_Gateway_SuperEnhancers.bed",
                                     format="bed")

mcols(Naive_enh)$type <- "E"
mcols(Naive_enh)$type[mcols(Naive_enh)$name %in% mcols(Naive_suprenh)$name] <- "SE"

Primed_enh <- rtracklayer::import("/mnt/Projects/CHiC_naive_primed/ROSE/Primed_1500/Primed_H3K27ac_peaks_narrowPeak_Gateway_Enhancers.bed",
                                  format="bed")
Primed_suprenh <- rtracklayer::import("/mnt/Projects/CHiC_naive_primed/ROSE/Primed_1500/Primed_H3K27ac_peaks_narrowPeak_Gateway_SuperEnhancers.bed",
                                      format="bed")

mcols(Primed_enh)$type <- "E"
mcols(Primed_enh)$type[mcols(Primed_enh)$name %in% mcols(Primed_suprenh)$name] <- "SE"

# Make a lookup of shared primed specific naive specific enhancers
Naive_enh$name <- gsub("_naive_peak", "_np", gsub("_lociStitched", "", Naive_enh$name))
Primed_enh$name <- gsub("_Primed_H3K27ac_peak", "_pp", gsub("_lociStitched", "", Primed_enh$name))

np_ovlp <- findOverlaps(Naive_enh, Primed_enh)

# rename overlapping enhancers to shared
Naive_enh$name[queryHits(np_ovlp)] <- gsub("np", "nsp", Naive_enh$name[queryHits(np_ovlp)])
Primed_enh$name[subjectHits(np_ovlp)] <- gsub("pp", "psp", Primed_enh$name[subjectHits(np_ovlp)])


# add an identifier for overlapping enhancers
overlap_names <- aggregate_hits(Naive_enh, Primed_enh)
Naive_enh <- overlap_names$gr
Primed_enh <- overlap_names$gr2

Naive_enh$name <- paste(Naive_enh$name, Naive_enh$type, sep = "_")
Primed_enh$name <- paste(Primed_enh$name, Primed_enh$type, sep = "_")

Naive_enh$name <- paste(Naive_enh$name, Naive_enh$overlap, sep = "-")
Primed_enh$name <- paste(Primed_enh$name, Primed_enh$overlap, sep = "-")


#Overlap with OSN peaks (don't do this on a hindiii fragment, but on a peak resolution)
# 1 = naive; 2 = shared; 3 = primed
OSN_all <- rtracklayer::import("/mnt/Projects/CHiC_naive_primed/ChIP-seq/OSN/osn_num.bed",
                               format = "bed")

OSN_naive <- OSN_all[mcols(OSN_all)$name %in% c(1, 2)]
OSN_primed <- OSN_all[mcols(OSN_all)$name %in% c(2, 3)]


Primed_enh_hindiii_OSN <- ROSE_OSN_HindIII_table(Primed_enh, hindiii_gr, OSN_primed)
Naive_enh_hindiii_OSN <- ROSE_OSN_HindIII_table(Naive_enh, hindiii_gr, OSN_naive)

#enhancer lookup for OSN (has both naive and primed OSN overlaping naive primed enhancers)

enh_OSN_ovl <- findOverlaps(c(Primed_enh, Naive_enh), OSN_all)

OSN_enhancer_df <- data.frame(name=c(Primed_enh, Naive_enh)[queryHits(enh_OSN_ovl)]$name, 
                              OSN=OSN_all[subjectHits(enh_OSN_ovl)]$name)

OSN_enhancer_lookup <- aggregate(x = OSN_enhancer_df["OSN"], 
                                 by= OSN_enhancer_df["name"], 
                                 FUN = function(ids){
                                   paste(ids, collapse = ",")
                                 })

OSN_enhancer_lookup$origin <- NA
OSN_enhancer_lookup$origin[grepl("np|nsp", OSN_enhancer_lookup$name)] <- "naive"
OSN_enhancer_lookup$origin[grepl("pp|psp", OSN_enhancer_lookup$name)] <- "primed"

# OSN_enhancer_lookup$short_name <- ifelse(grepl("naive", OSN_enhancer_lookup$name), 
#                                          gsub("_naive_peak", "_np", gsub("_lociStitched", "", OSN_enhancer_lookup$name)),
#                                          gsub("_Primed_H3K27ac_peak", "_pp", gsub("_lociStitched", "", OSN_enhancer_lookup$name)))


setDT(OSN_enhancer_lookup, key="name")

#tidyr
naive_enhancer_lookup <- Naive_enh_hindiii_OSN %>% separate_rows(ID, OSN, sep=";")
primed_enhancer_lookup <- Primed_enh_hindiii_OSN %>% separate_rows(ID, OSN, sep=";")

np_enh <- merge(naive_enhancer_lookup[c("ID", "name", "score", "type", "OSN")],
                primed_enhancer_lookup[c("ID", "name", "score", "type", "OSN")],
                by = "ID", all = TRUE, suffixes = c("_naive","_primed"))

# np_enh$name_naive <- gsub("_naive_peak", "_np", gsub("_lociStitched", "", np_enh$name_naive))
# np_enh$name_primed <- gsub("_Primed_H3K27ac_peak", "_pp", gsub("_lociStitched", "", np_enh$name_primed))
#by hindIII fragments
np_enh_OSN_lookup <- aggregate(np_enh[,-1], by=list(np_enh$ID), 
                               FUN = function(x) {
                                 if(all(is.na(x))){
                                   return(NA)
                                 }else{
                                   paste(x, collapse = ";")
                                 }
                               })



np_enh_OSN_lookup <- deduplicated_collapsed(np_enh_OSN_lookup, "name_naive", c("name_naive", "score_naive", "type_naive", "OSN_naive"))      
np_enh_OSN_lookup <- deduplicated_collapsed(np_enh_OSN_lookup, "name_primed", c("name_primed", "score_primed", "type_primed", "OSN_primed"))      

names(np_enh_OSN_lookup)[1] <- "ID"

saveRDS(np_enh_OSN_lookup, "lookup_files/np_enh_OSN_lookup_20200915.rds")


#combine enhancer and chromhmm data
#capitalise X in chrx
np_enh_OSN_lookup$ID <- sapply(np_enh_OSN_lookup$ID, function(x) gsub("x", "X", x))
#convert NA to real NA
np_enh_OSN_lookup[np_enh_OSN_lookup == "NA"] <- NA

chromhmm_enh_OSN_lookup <- merge(chromhmm_lookup, np_enh_OSN_lookup, by="ID", all=TRUE)
setDT(chromhmm_enh_OSN_lookup, key = "ID")

saveRDS(chromhmm_enh_OSN_lookup, "lookup_files/chromhmm_enh_OSN_lookup_20200915.rds")



