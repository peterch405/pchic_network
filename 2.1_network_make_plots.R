library(readr)
library(igraph)
library(ggplot2)
library(ggpubr)


source("network_dynamics_functions.R")



#Naive primed only -------------------------------------------------------------
net_all_s <- readRDS("2_network_make/net_all_s_20190829.rds")

net_all_naive <- isolate_cell_network(net_all_s, "naive")
net_all_primed <- isolate_cell_network(net_all_s, "primed")


#Make plots based on subnetwork id
# KRT cluster 1233
# Olfactory 749
# Hist 2609
# Pcdh 2387

#'Plot paired boxplot of the node degrees of subnetworks
#'
#'
degree_subnetwork <- function(net_all_naive, net_all_primed, net_number, gene_only=TRUE, gene_regex=NA){
  degree_naive <- data.frame(name=V(net_all_naive)$name[V(net_all_naive)$subnet == net_number],
                             gene=V(net_all_naive)$prot_genes[V(net_all_naive)$subnet == net_number],
                             count=unname(degree(net_all_naive, V(net_all_naive)$subnet == net_number)))
  degree_primed <- data.frame(name=V(net_all_primed)$name[V(net_all_primed)$subnet == net_number],
                              gene=V(net_all_primed)$prot_genes[V(net_all_primed)$subnet == net_number],
                              count=unname(degree(net_all_primed, V(net_all_primed)$subnet == net_number)))
  
  stopifnot(nrow(degree_naive) > 0 & nrow(degree_primed) > 0)
  
  degree_naive$origin <- "naive"
  degree_primed$origin <- "primed"
  
  if(gene_only){
    
    degree_naive <- degree_naive[!is.na(degree_naive$gene),]
    degree_primed <- degree_primed[!is.na(degree_primed$gene),]
  }
  
  if(!is.na(gene_regex)){
    degree_naive <- degree_naive[grepl(gene_regex, degree_naive$gene,  ignore.case = TRUE),]
    degree_primed <- degree_primed[grepl(gene_regex, degree_primed$gene,  ignore.case = TRUE),]
    
  }
  
  degree_all <- merge(degree_naive, degree_primed, by = "name", all=TRUE, suffixes = c("_naive","_primed"))
  degree_all$count_naive[is.na(degree_all$count_naive)] <- 0
  degree_all$count_primed[is.na(degree_all$count_primed)] <- 0
  
  degree_all$origin_naive <- "naive"
  degree_all$origin_primed <- "primed"

  
  degree_all <- data.frame(name=c(as.character(degree_all$name), as.character(degree_all$name)),
                           count=c(degree_all$count_naive, degree_all$count_primed),
                           origin=c(degree_all$origin_naive, degree_all$origin_primed),
                           gene=c(as.character(degree_all$gene_naive), as.character(degree_all$gene_primed)))
  
  
  
  ggpaired(degree_all, x = "origin", y = "count",
           color = "origin", ylab="gene degree", xlab=NULL,
           line.color = "gray",  line.size = 0.2, point.size = 0.5) +
           stat_compare_means(paired=TRUE)

}

degree_subnetwork(net_all_naive, net_all_primed, 2387, gene_only = FALSE, gene_regex = "PCDH")
ggsave("2_network_make/genes_degree_stat_pcdh_cluster.pdf", device="pdf",
       width = 6, height = 8, units="cm", useDingbats=FALSE)

degree_subnetwork(net_all_naive, net_all_primed, 2609, gene_only = FALSE, gene_regex = "HIST")
ggsave("2_network_make/genes_degree_stat_hist_cluster.pdf", device="pdf",
       width = 6, height = 8, units="cm", useDingbats=FALSE)

degree_subnetwork(net_all_naive, net_all_primed, 1233, gene_only = FALSE, gene_regex = "^KRT\\d")
ggsave("2_network_make/genes_degree_stat_krt_cluster.pdf", device="pdf",
       width = 6, height = 8, units="cm", useDingbats=FALSE)

degree_subnetwork(net_all_naive, net_all_primed, 749, gene_only = FALSE, gene_regex = "^OR\\d")
ggsave("2_network_make/genes_degree_stat_or_cluster.pdf", device="pdf",
       width = 6, height = 8, units="cm", useDingbats=FALSE)


