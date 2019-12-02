library(igraph)
library(ggplot2)
library(naniar)
library(ggpubr)

de_results <- readRDS("1_DESeq2_gene_expression/de_results.rds")


de_results$colour <- "NS"
de_results$colour[de_results$padj < 0.05] <- "P"
de_results$colour[de_results$log2FoldChange > 1.5 | de_results$log2FoldChange < -1.5] <- "FC"
de_results$colour[(de_results$log2FoldChange > 1.5 | de_results$log2FoldChange < -1.5) & de_results$padj < 0.05] <- "FC_P"

de_results$labels <- ""
de_results$labels[de_results$colour == "FC_P"] <- de_results$gene_name[de_results$colour == "FC_P"]


net_all_s <- readRDS("1_DESeq2_gene_expression/net_all_s_20190829.rds")


# KRT cluster 1233
# Olfactory 749
# Hist 2609
# Pcdh 2387

#' Plot expression of select genes in a subnetwork
#' 
plot_subnetwork_expression <- function(net_all_s, net_number, gene_regex){
  
  genes_net <- data.frame(name=V(net_all_s)$name[V(net_all_s)$subnet == net_number],
                          gene=V(net_all_s)$prot_genes[V(net_all_s)$subnet == net_number])
  net_genes <- unlist(strsplit(as.character(genes_net$gene), ","))
  select_net <-net_genes[grepl(gene_regex, net_genes,  ignore.case = TRUE)]
  
  
  to_plot <- de_results[de_results$gene_name %in% select_net,]
  #consider genes not present as not expressed
  to_plot <- merge(to_plot, data.frame(gene_name=select_net), by="gene_name", all=TRUE)
  to_plot$colour[is.na(to_plot$colour)] <- "NE"
  to_plot$log2FoldChange[is.na(to_plot$log2FoldChange)] <- runif(length(to_plot$log2FoldChange[is.na(to_plot$log2FoldChange)]), 
                                                                 min=-1, max=1)
  print(select_net)
  ggplot(to_plot, aes(x=log2FoldChange, y=-log10(padj), colour = colour)) + 
    geom_miss_point(alpha=0.5) +
    # geom_point(alpha=0.5) +
    geom_vline(xintercept = 1.5, linetype = 2)+
    geom_vline(xintercept = -1.5, linetype = 2)+
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    ggtitle("Volcano Plot Naive Primed via DESeq2") + 
    theme_pubr() +
    scale_color_manual(values=c(NS="grey30", FC="forestgreen", P="royalblue", FC_P="red2", NE="#EDEDED")) +
    xlim(-10,10)
}


plot_subnetwork_expression(net_all_s, 2609, "HIST")
ggsave("1_DESeq2_gene_expression/volcano_hist_10_subnet.pdf", width = 4, height = 4)
plot_subnetwork_expression(net_all_s, 2387, "PCDH")
ggsave("1_DESeq2_gene_expression/volcano_pcdh_10_subnet.pdf", width = 4, height = 4)

plot_subnetwork_expression(net_all_s, 1233, "^KRT\\d")
ggsave("1_DESeq2_gene_expression/volcano_krt_10_subnet.pdf", width = 4, height = 4)

plot_subnetwork_expression(net_all_s, 749, "^OR\\d")
ggsave("1_DESeq2_gene_expression/volcano_or_10_subnet.pdf", width = 4, height = 4)
