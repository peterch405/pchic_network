library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(Sushi)
library(igraph)
library(data.table)
library(rtracklayer)

source("network_dynamics_functions.R")

np_summary_deb2b <- readRDS("3_network_dynamics/np_summary_deb2b_20190829.rds")


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

net_all_s <- readRDS("3_network_dynamics/net_all_s_20190829.rds")
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
primed_bedpe_658 <- get_bedpe(all_subgraphs[[658]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_658 <- get_bedpe(all_subgraphs[[658]], hindiii_mid_lookup, "primed", mark_dist = 2^20)

#HOXD
primed_bedpe_1836 <- get_bedpe(all_subgraphs[[1836]], hindiii_mid_lookup, "naive", mark_dist = 2^20)
naive_bedpe_1836 <- get_bedpe(all_subgraphs[[1836]], hindiii_mid_lookup, "primed", mark_dist = 2^20)


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

plot_subnetwork_sushi(primed_bedpe_658, naive_bedpe_658, hg38_chrom,
                      "3_network_dynamics/658_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_1836, naive_bedpe_1836, hg38_chrom,
                      "3_network_dynamics/1836_changing_sushi_2^20_sax.pdf",
                      zoom=FALSE)

plot_subnetwork_sushi(primed_bedpe_1836, naive_bedpe_1836, hg38_chrom,
                      "3_network_dynamics/1836_changing_sushi_2^20_sax_zoom.pdf",
                      zoom=TRUE)

plot_subnetwork_sushi(primed_bedpe_2301, naive_bedpe_2301, hg38_chrom, zoom=FALSE)

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
#initialise progress bar
pb <- tkProgressBar("Progress bar", "Some information in %",
                    0, length(all_subgraphs), 0, width=300)


all_net_list <- list()
moduty <- list()
subnet_communities <- list()
subnet_comm <- list()
for(sg in seq_along(all_subgraphs)){
  sub_graph <- all_subgraphs[[sg]]
  comms <- multilevel.community(as.undirected(sub_graph), weights = E(sub_graph)$weight)
  
  moduty[[sg]] <- data.frame(mod=modularity(comms), groups=length(comms), e=length(E(sub_graph)), n=length(V(sub_graph)))
  
  #naive or primed subnetwork
  p_net <- cell_net(all_subgraphs[[sg]], "naive")
  n_net <- cell_net(all_subgraphs[[sg]], "primed")
  
  V(p_net)$btwn.cent <- betweenness(as.undirected(p_net))
  V(n_net)$btwn.cent <- betweenness(as.undirected(n_net))
  
  
  #add also naive and primed specific counts
  if(modularity(comms) >= 0.7){ #if modularity is high else just take entire netowork
    
    #loop through groups identified, save nodes within and give unique name to subnetwork
    for(c in seq_along(comms)){
      
      p_del <- delete.vertices(p_net, !(V(p_net)$name %in% comms[[c]]))
      n_del <- delete.vertices(n_net, !(V(n_net)$name %in% comms[[c]]))
      
      #include Robustness measure
      
      subnet_communities[[paste(sg, c, sep="_")]] <- data.frame(nodes=paste(comms[[c]], collapse = ";"),
                                                                naive_n=vcount(n_del),
                                                                primed_n=vcount(p_del),
                                                                naive_e=ecount(n_del),
                                                                primed_e=ecount(p_del),
                                                                modul=modularity(comms),
                                                                p_btwn_cent=paste(V(p_del)$btwn.cent, collapse = ";"),
                                                                n_btwn_cent=paste(V(n_del)$btwn.cent, collapse = ";"))
      #simple list for asigning in network
      subnet_comm[[paste(sg, c, sep="_")]] <- comms[[c]]
      
      
    }
  }else{
    subnet_communities[[as.character(sg)]] <- data.frame(nodes=paste(V(all_subgraphs[[sg]])$name, collapse = ";"),
                                                         naive_n=vcount(n_net),
                                                         primed_n=vcount(p_net),
                                                         naive_e=ecount(n_net),
                                                         primed_e=ecount(p_net),
                                                         modul=modularity(comms),
                                                         p_btwn_cent=paste(V(p_net)$btwn.cent, collapse = ";"),
                                                         n_btwn_cent=paste(V(n_net)$btwn.cent, collapse = ";"))
    
    subnet_comm[[as.character(sg)]] <- V(all_subgraphs[[sg]])$name
    
  }
  #progress bar
  info <- sprintf("%d%% done", (round(sg/length(all_subgraphs)*100)))
  setTkProgressBar(pb, sg, sprintf("Outer loop (%s)", info), info)
  
}

net_modularity <- plyr::ldply(moduty, rbind)
net_modularity$.id <- as.numeric(rownames(net_modularity))

net_communities <- plyr::ldply(subnet_communities, rbind)

#close progress bar
close(pb)


#remove nodes with only one central node (low robustness)
table(cut(table(as.numeric(unlist(strsplit(as.character(net_communities$p_btwn_cent), ";")))), seq(0,500, 10)))

net_communities$p_btwn_cent_cnt <- sapply(strsplit(as.character(net_communities$p_btwn_cent), ";"), function(x) sum((as.numeric(x) > 0)*1))
net_communities$n_btwn_cent_cnt <- sapply(strsplit(as.character(net_communities$n_btwn_cent), ";"), function(x) sum((as.numeric(x) > 0)*1))

#remove ones with btwn.cent 0 and 1

net_communities_sub <- net_communities[net_communities$p_btwn_cent_cnt > 1 & net_communities$n_btwn_cent_cnt > 1,]
net_communities_dt <- as.data.table(net_communities_sub) #lookup nodes
setkey(net_communities_dt, ".id")

saveRDS(net_communities_dt, "/media/chovanec/My_Passport/CHiC_naive_primed/TAD/net_communities_dt.rds")


