library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(Sushi)
library(igraph)
library(data.table)
library(rtracklayer)
library(tcltk2)

source("network_dynamics_functions.R")

np_summary_deb2b <- readRDS("2_network_make/np_summary_deb2b_20200911.rds")


ggplot(np_summary_deb2b, aes(x = v, y = e)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  geom_point(data=np_summary_deb2b[np_summary_deb2b$.id %in% c("2301", "2609", "2387"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
ggsave("3_network_dynamics/break_networks_non_merged_ve_deb2b.pdf", 
       device="pdf", width = 10, height = 10, units="cm")

#Subnetwork plots --------------------------------------------------------------

net_all_s <- readRDS("2_network_make/net_all_s_20200911.rds")
E(net_all_s)$weight <- E(net_all_s)$distance
all_subgraphs <- decompose.graph(net_all_s)

hindiii_mid_lookup <- readRDS("3_network_dynamics/hindiii_mid_lookup.rds")

#HIST
primed_bedpe_2609 <- get_bedpe(all_subgraphs[[2609]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_2609 <- get_bedpe(all_subgraphs[[2609]], hindiii_mid_lookup, "primed", mark_dist = 2^20)

#Largest
primed_bedpe_2301 <- get_bedpe(all_subgraphs[[2301]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_2301 <- get_bedpe(all_subgraphs[[2301]], hindiii_mid_lookup, "primed", mark_dist = 2^20)

#HOXA
primed_bedpe_2650 <- get_bedpe(all_subgraphs[[2650]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_2650 <- get_bedpe(all_subgraphs[[2650]], hindiii_mid_lookup, "primed", mark_dist = 2^20)

#NKX
primed_bedpe_755 <- get_bedpe(all_subgraphs[[755]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_755 <- get_bedpe(all_subgraphs[[755]], hindiii_mid_lookup, "primed", mark_dist = 2^20)

#HOXD
primed_bedpe_1591 <- get_bedpe(all_subgraphs[[1591]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_1591 <- get_bedpe(all_subgraphs[[1591]], hindiii_mid_lookup, "primed", mark_dist = 2^20)



hg38_chrom <- readRDS("3_network_dynamics/hg38_chrom.rds")




plot_subnetwork_sushi(primed_bedpe_2609, naive_bedpe_2609, hg38_chrom,
                      "3_network_dynamics/2609_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_2609, naive_bedpe_2609, hg38_chrom,
                      "3_network_dynamics/2609_changing_sushi_2^20_sax_zoom.pdf",
                      zoom=TRUE)

plot_subnetwork_sushi(primed_bedpe_2301, naive_bedpe_2301, hg38_chrom,
                      "3_network_dynamics/2301_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_2650, naive_bedpe_2650, hg38_chrom,
                      "3_network_dynamics/2650_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_755, naive_bedpe_755, hg38_chrom,
                      "3_network_dynamics/755_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_1591, naive_bedpe_1591, hg38_chrom,
                      "3_network_dynamics/1591_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_1591, naive_bedpe_1591, hg38_chrom,
                      "3_network_dynamics/1591_changing_sushi_2^20_sax_zoom.pdf",
                      zoom=TRUE)

H3K27me3_naive <- rtracklayer::import.bw("3_network_dynamics/H3K27me3_Naive_SRR1515138_minus_input_chr_blacklist_noneg.bigwig")
H3K27me3_primed <- rtracklayer::import.bw("3_network_dynamics/H3K27me3_Primed_SRR1515137_minus_input_chr_blacklist_noneg.bigwig")

H3K27me3_naive_df <- data.frame(H3K27me3_naive[seqnames(H3K27me3_naive) == "chr5"])
H3K27me3_naive_df <- H3K27me3_naive_df[,!(names(H3K27me3_naive_df) %in% c("width", "strand"))]
H3K27me3_naive_df$start <- H3K27me3_naive_df$start-1
H3K27me3_primed_df <- data.frame(H3K27me3_primed[seqnames(H3K27me3_primed) == "chr5"])
H3K27me3_primed_df <- H3K27me3_primed_df[,!(names(H3K27me3_primed_df) %in% c("width", "strand"))]
H3K27me3_primed_df$start <- H3K27me3_primed_df$start-1



pdf("3_network_dynamics/H3K27me3_chr5.pdf", width=10, height=5)
par(mfrow = c(2,1))
chrom            = "chr5"
chromstart       = 1
chromend         = 181538259
plotBedgraph(H3K27me3_naive_df,chrom,chromstart,chromend,color = "#3C96EB", range = c(0,40))
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

plotBedgraph(H3K27me3_primed_df,chrom,chromstart,chromend,color = "#C00F14", range = c(0, 40))
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
dev.off()



pca_naive <- rtracklayer::import.bedGraph("3_network_dynamics/PCA_250kb.washu_edit_naive.bedgraph")
pca_primed <- rtracklayer::import.bedGraph("3_network_dynamics/PCA_250kb.washu_edit_primed.bedgraph")

pca_naive_df <- data.frame(pca_naive[seqnames(pca_naive)=="chr5"])
pca_naive_df <- pca_naive_df[,!(names(pca_naive_df) %in% c("width", "strand"))]
pca_naive_df$start <- pca_naive_df$start-1

pca_primed_df <- data.frame(pca_primed[seqnames(pca_primed)=="chr5"])
pca_primed_df <- pca_primed_df[,!(names(pca_primed_df) %in% c("width", "strand"))]
pca_primed_df$start <- pca_primed_df$start-1


pdf("3_network_dynamics/PCA_AB_chr5.pdf", width=10, height=5)
plotBedgraph(pca_naive_df,chrom,chromstart,chromend, color = "#3C96EB", range = c(-43,43), transparency=0.5, lwd=1)
plotBedgraph(pca_primed_df,chrom,chromstart,chromend, color = "#C00F14", overlay=TRUE, rescaleoverlay=TRUE, transparency=0.5, lwd=1)
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
dev.off()

#coregulation analysis based on chromhmm states --------------------------------
net_communities_dt <- community_level_dynamics(net_all_s)
net_communities_dt$v <- net_communities_dt$naive_n-net_communities_dt$primed_n
net_communities_dt$e <- net_communities_dt$naive_e-net_communities_dt$primed_e

prot_gene_lookup <- readRDS("lookup_files/nodes_all_prot.rds")
cols <- c("ID", "prot_genes")
prot_gene_lookup <- prot_gene_lookup[, ..cols]

net_nodes <- lapply(net_communities_dt$nodes, function(x) unlist(strsplit(x, ";")))
net_genes <- lapply(net_nodes, function(x) paste(unique(na.omit(prot_gene_lookup[x])$prot_genes), collapse = ";"))
net_communities_dt$prot_genes <- unlist(net_genes)

ggplot(net_communities_dt, aes(x = v, y = e)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  geom_point(data=net_communities_dt[net_communities_dt$.id %in% c("2609_3", "2387"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
ggsave("3_network_dynamics/break_networks_non_merged_ve_deb2b_community.pdf", 
       device="pdf", width = 10, height = 10, units="cm")

#Used in TAD community overlap analysis
saveRDS(net_communities_dt, "3_network_dynamics/net_communities_dt.rds")

# Bed files for pile-up plots --------------------------------------------------

#' Coarsen coordinates 
#' 
#' @param coordinate vector with start or end coordinate
#' @param resolution in bp
coarsen_coordinates <- function(coordinates, resolution){
  coord <- sapply(coordinates, function(c) as.integer(as.integer(c)/resolution)*resolution)
  return(coord)
}

c_cols <- c("start_1", "end_1", "start_2", "end_2")

primed_bedpe_2301_25k <- primed_bedpe_2301
primed_bedpe_2301_25k[c_cols] <- lapply(primed_bedpe_2301_25k[c_cols], coarsen_coordinates, resolution=25000)
primed_bedpe_2301_25k[c("chrom_1", "chrom_2")] <- lapply(primed_bedpe_2301_25k[c("chrom_1", "chrom_2")], function(x) gsub("chr", "", x))
# Export 6 column bed for pile-up plots
write.table(primed_bedpe_2301_25k[1:6], "3_network_dynamics/primed_2301_25k.bed", quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)

naive_bedpe_2301_25k <- naive_bedpe_2301
naive_bedpe_2301_25k[c_cols] <- lapply(naive_bedpe_2301_25k[c_cols], coarsen_coordinates, resolution=25000)
naive_bedpe_2301_25k[c("chrom_1", "chrom_2")] <- lapply(naive_bedpe_2301_25k[c("chrom_1", "chrom_2")], function(x) gsub("chr", "", x))
# Export 6 column bed for pile-up plots
write.table(naive_bedpe_2301_25k[1:6], "3_network_dynamics/naive_2301_25k.bed", quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)

