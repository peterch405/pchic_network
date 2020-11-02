library(readr)
require(igraph)
library(tcltk2)
library(plyr)
library(dplyr)
library(data.table)
library(Sushi)
library(rjson)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(plotly)

source("network_dynamics_functions.R")
source("network_make_functions.R")
source("ExportImportGraph.R")

#Load RNA-seq data -------------------------------------------------------------

de_genes <- read_delim("1_DESeq2_gene_expression/de_genes_takashima_GRCh38.87_anno_opposing_strand_prot_genes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)



#Load interaction data ---------------------------------------------------------

naive_nodes <- read_delim("0_network_data_preparation/naive_nodes_wchip_20181218.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

primed_nodes <- read_delim("0_network_data_preparation/primed_nodes_wchip_20181218.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

naive_intract_rep1 <- read_delim("/mnt/Projects/CHiC_naive_primed/Chicago_new/naive_rep1_intract_wchip.txt", 
                            "\t", escape_double = FALSE, col_names = TRUE, 
                            trim_ws = TRUE)

primed_intract_rep1 <- read_delim("/mnt/Projects/CHiC_naive_primed/Chicago_new/primed_rep1_intract_wchip.txt", 
                             "\t", escape_double = FALSE, col_names = TRUE, 
                             trim_ws = TRUE)

naive_intract_rep2 <- read_delim("/mnt/Projects/CHiC_naive_primed/Chicago_new/naive_rep2_intract_wchip.txt", 
                                 "\t", escape_double = FALSE, col_names = TRUE, 
                                 trim_ws = TRUE)

primed_intract_rep2 <- read_delim("/mnt/Projects/CHiC_naive_primed/Chicago_new/primed_rep2_intract_wchip.txt", 
                                  "\t", escape_double = FALSE, col_names = TRUE, 
                                  trim_ws = TRUE)


# Custom chromosome colour palette ---------------------------------------------

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,24), col=col_vector[c(5,6,10,11,13,14,15,16,17,18,19,20,21,22,25,26,29,49,54,55,56,57,71,74)])

custom_chrom_col <- apply(col2rgb(col_vector[c(5,6,10,11,13,14,15,16,17,18,19,20,21,22,25,26,29,49,54,55,56,57,71,74)]), 2, function(x) paste(x, collapse = ","))

#Make network Rep1 -------------------------------------------------------------

links_nodes_deb2b <- make_link_node_dfs(naive_nodes, primed_nodes, naive_intract_rep1, primed_intract_rep1, remove_b2b = TRUE)

#only protein coding genes
links_nodes_cat_deb2b <- categorise_add_data(links_nodes_deb2b, de_genes)

#Colour each tads with the same colour between naive and primed
links_nodes_cat_col_deb2b <- add_same_tad_colour(links_nodes_cat_deb2b)

#Enhancer categories
links_nodes_cal_col_bed2b_enh <- categorise_enhancers(links_nodes_cat_col_deb2b)

#Add Gephi layout coordinates (to create reproducible layouts, this has to be exported from gephi)
# links_nodes_cat_col_coord_deb2b <- add_layout_coordinate(links_nodes_cal_col_bed2b_enh, 
#                                                          "2_network_make/data.json")

#Make network 
net_out_deb2b <- make_network(links_nodes_cal_col_bed2b_enh, "2_network_make/np_summary_deb2b_rep1.rds",
                              directed = FALSE, save_nontrans_net = "2_network_make/net_all_s_rep1.rds",
                              chrom_colours=custom_chrom_col)

igraph::write_graph(net_out_deb2b, "2_network_make/net_all_trans_s_deb2b_rep1.graphml", 
                    format = "graphml")

saveRDS(links_nodes_cal_col_bed2b_enh,"2_network_make/links_nodes_cat_col_coord_deb2b_rep1.rds")

#Make network Rep2 -------------------------------------------------------------

links_nodes_deb2b <- make_link_node_dfs(naive_nodes, primed_nodes, naive_intract_rep2, primed_intract_rep2, remove_b2b = TRUE)

#only protein coding genes
links_nodes_cat_deb2b <- categorise_add_data(links_nodes_deb2b, de_genes)

#Colour each tads with the same colour between naive and primed
links_nodes_cat_col_deb2b <- add_same_tad_colour(links_nodes_cat_deb2b)

#Enhancer categories
links_nodes_cal_col_bed2b_enh <- categorise_enhancers(links_nodes_cat_col_deb2b)

#Add Gephi layout coordinates (to create reproducible layouts, this has to be exported from gephi)
# links_nodes_cat_col_coord_deb2b <- add_layout_coordinate(links_nodes_cal_col_bed2b_enh, 
#                                                          "2_network_make/data.json")

#Make network 
net_out_deb2b <- make_network(links_nodes_cal_col_bed2b_enh, "2_network_make/np_summary_deb2b_rep2.rds",
                              directed = FALSE, save_nontrans_net = "2_network_make/net_all_s_rep2.rds",
                              chrom_colours=custom_chrom_col)

igraph::write_graph(net_out_deb2b, "2_network_make/net_all_trans_s_deb2b_rep2.graphml", 
                    format = "graphml")

saveRDS(links_nodes_cal_col_bed2b_enh,"2_network_make/links_nodes_cat_col_coord_deb2b_rep2.rds")

# Reproducibility of network layouts between replicates ------------------------


rep1 <- read_layout_coordinate("2_network_make/Rep1.json")
rep2 <- read_layout_coordinate("2_network_make/Rep2.json")
rep1_cl2 <- read_layout_coordinate("2_network_make/Rep1_chrom_layout2.json")
rep1_l2 <- read_layout_coordinate("2_network_make/Rep1_layout2.json")

rep1_filt <- read_layout_coordinate("2_network_make/Rep1_chrom_layout_filt.json")
rep2_filt <- read_layout_coordinate("2_network_make/Rep2_chrom_layout_filt.json")


#' Allows to convert a json into an igraph object
#' It is essential that the json contains an igraph object
#' 
#' @param filename either the json or the filename of the file containing the json
importGraph <- function(filename){
  built.graph <-  rjson::fromJSON(file=filename)
  if("vertices" %in% names(built.graph)){
    built.g <- graph_from_data_frame(built.graph$edges, directed=built.graph$directed, 
                                     vertices=built.graph$vertices, x=build.graph$nodes$x,
                                     y=build.graph$nodes$y,
                                     name=build.graph$nodes$attributes$name)
  }else{
    built.g <- graph_from_data_frame(built.graph$edges, directed=built.graph$directed)
  }
  if("name" %in% names(built.graph)){
    built.g$name <- built.graph$name
  }
  return(built.g)
}

rep1_filt_net <- importGraph("2_network_make/Rep1_chrom_layout_filt.json")



# make coordinate lookup
rep1_list <- split(rep1[,2:3], seq(nrow(rep1)))
names(rep1_list) <- rep1$ID
# RStudio stack overflow error
rep1_lookup <- list2env(rep1_list)

rep2_list <- split(rep2[,2:3], seq(nrow(rep2)))
names(rep2_list) <- rep2$ID
# RStudio stack overflow error
rep2_lookup <- list2env(rep2_list)

rep1_cl2_list <- split(rep1_cl2[,2:3], seq(nrow(rep1_cl2)))
names(rep1_cl2_list) <- rep1_cl2$ID
# RStudio stack overflow error
rep1_cl2_lookup <- list2env(rep1_cl2_list)

rep1_l2_list <- split(rep1_l2[,2:3], seq(nrow(rep1_l2)))
names(rep1_l2_list) <- rep1_l2$ID
# RStudio stack overflow error
rep1_l2_lookup <- list2env(rep1_l2_list)


rep1_filt_list <- split(rep1_filt[,2:3], seq(nrow(rep1_filt)))
names(rep1_filt_list) <- rep1_filt$ID
# RStudio stack overflow error
rep1_filt_lookup <- list2env(rep1_filt_list)


rep2_filt_list <- split(rep2_filt[,2:3], seq(nrow(rep2_filt)))
names(rep2_filt_list) <- rep2_filt$ID
# RStudio stack overflow error
rep2_filt_lookup <- list2env(rep2_filt_list)

rep1_lookup <- readRDS("2_network_make/rep1_lookup.rds")
rep2_lookup <- readRDS("2_network_make/rep2_lookup.rds")


#' Read gephi layout coordinates from json file
#' 
#' @param gephi_json path to json file
#' @return data.frame with node name and x y coordinates
read_layout_coordinate <- function(gephi_json){
  
  graph_json <- rjson::fromJSON(file=gephi_json)
  
  graph_xy <- list()
  for(i in seq_along(graph_json$nodes)){
    graph_xy[[i]] <- data.frame(name=graph_json$nodes[[i]]$attributes$name,
                                x=graph_json$nodes[[i]]$x,
                                y=graph_json$nodes[[i]]$y)
    
  }
  
  graph_xy_df <- plyr::ldply(graph_xy, rbind)
  names(graph_xy_df)[1] <- "ID"
  
  return(graph_xy_df)
}

#' Calculate distance between x y coordinates
#' 
#' @param a x y coordinates of point a tuple
#' @param b x y coordinates of point b tuple
#' @return distance vector
coordinate_distance <- function(a, b){
  stopifnot(length(a) == length(b))
  s1 <- abs(a[1] - b[1])
  s2 <- abs(a[2] - b[2])

  return(sqrt(s1^2 + s2^2))
}

# Get shared links between rep1 and rep2 
net_rep1 <- readRDS("2_network_make/net_all_s_rep1.rds")
net_rep2 <- readRDS("2_network_make/net_all_s_rep2.rds")
net_data_rep1 <- igraph::as_long_data_frame(net_rep1)
net_data_rep2 <- igraph::as_long_data_frame(net_rep2)
setDT(net_data_rep1, key = c("from_name", "to_name"))
setDT(net_data_rep2, key = c("from_name", "to_name"))
#reduce to only needed columns
cols <- c("from_name", "to_name")
net_data_rep1 <- net_data_rep1[, ..cols]
net_data_rep2 <- net_data_rep2[, ..cols]


rep1_link_id <- paste(net_data_rep1$from_name, net_data_rep1$to_name, sep = "_")
rep2_link_id <- paste(net_data_rep2$from_name, net_data_rep2$to_name, sep = "_")
shared_id <- intersect(rep1_link_id, rep2_link_id)
shared_links <- net_data_rep1[which(rep1_link_id %in% shared_id),]


# rep1_lookup[[shared_links$from_name[1]]]
# rep1_lookup[[shared_links$to_name[1]]]
# 
# rep2_lookup[[shared_links$from_name[1]]]
# rep2_lookup[[shared_links$to_name[1]]]



rep1_dist <- apply(shared_links, 1, function(x) 
                   coordinate_distance(rep1_lookup[[x[1]]], rep1_lookup[[x[2]]]))
rep2_dist <- apply(shared_links, 1, function(x) 
                   coordinate_distance(rep2_lookup[[x[1]]], rep2_lookup[[x[2]]]))

rep1_cl2_dist <- apply(shared_links, 1, function(x) 
  coordinate_distance(rep1_cl2_lookup[[x[1]]], rep1_cl2_lookup[[x[2]]]))
rep1_l2_dist <- apply(shared_links, 1, function(x) 
  coordinate_distance(rep1_l2_lookup[[x[1]]], rep1_l2_lookup[[x[2]]]))

shared_filt <- shared_links[shared_links$from_name %in% rep1_filt$ID & shared_links$to_name %in% rep1_filt$ID] 
shared_filt <- shared_filt[shared_filt$from_name %in% rep2_filt$ID & shared_filt$to_name %in% rep2_filt$ID] 

# for(i in 1:nrow(shared_filt)){
#   print(i)
#   coordinate_distance(rep1_filt_lookup[[shared_filt[[3,1]]]], rep1_filt_lookup[[shared_filt[[3,2]]]])
# }
  
rep1_filt_dist <- apply(shared_filt, 1, function(x) 
  coordinate_distance(rep1_filt_lookup[[x[1]]], rep1_filt_lookup[[x[2]]]))
rep2_filt_dist <- apply(shared_filt, 1, function(x) 
  coordinate_distance(rep2_filt_lookup[[x[1]]], rep2_filt_lookup[[x[2]]]))

rep_dist <- data.frame(rep1=unlist(rep1_dist), rep2=unlist(rep2_dist))
rep_dist <- data.frame(rep1=unlist(rep1_dist), rep2=unlist(rep1_cl2_dist))

rep_dist <- data.frame(rep1=unlist(rep1_dist), rep2=unlist(rep1_l2_dist))

rep_dist <- data.frame(rep1=unlist(rep1_filt_dist), rep2=unlist(rep2_filt_dist))


ggplot(rep_dist, aes(rep1, rep2)) +
  geom_point(alpha=0.05) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    ggpmisc::stat_poly_eq(formula = y ~ x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +
  theme_bw() +
  coord_trans(x="log2", y="log2")

rep_dist <- rep_dist[!(rep_dist$rep1 > 200 | rep_dist$rep2 > 200),]

cor(unlist(rep1_dist), unlist(rep2_dist), method = c("pearson"))
cor(unlist(rep1_dist), unlist(rep1_cl2_dist), method = c("pearson"))
cor(unlist(rep1_dist), unlist(rep1_l2_dist), method = c("pearson"))
cor(unlist(rep1_filt_dist), unlist(rep2_filt_dist), method = c("pearson"))


# Subnetwork changes plots -----------------------------------------------------
np_summary_deb2b_rep1 <- readRDS("2_network_make/np_summary_deb2b_rep1.rds")
np_summary_deb2b_rep2 <- readRDS("2_network_make/np_summary_deb2b_rep2.rds")


ggplot(np_summary_deb2b_rep1, aes(x = v, y = e)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  geom_point(data=np_summary_deb2b_rep1[np_summary_deb2b_rep1$.id %in% c("5596", "5123"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
ggsave("2_network_make/Rep1_network_change.pdf", width = 5, height = 5)

ggplot(np_summary_deb2b_rep2, aes(x = v, y = e)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  geom_point(data=np_summary_deb2b_rep2[np_summary_deb2b_rep2$.id %in% c("5267", "4817"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
ggsave("2_network_make/Rep2_network_change.pdf", width = 5, height = 5)

# Community level changes plots ------------------------------------------------
net_all_s_rep1 <- readRDS("2_network_make/net_all_s_rep1.rds")
net_summary_rep1 <- community_level_dynamics(net_all_s_rep1)
net_summary_rep1$v <- net_summary_rep1$naive_n-net_summary_rep1$primed_n
net_summary_rep1$e <- net_summary_rep1$naive_e-net_summary_rep1$primed_e
net_summary_rep1$origin <- "rep1"

net_all_s_rep2 <- readRDS("2_network_make/net_all_s_rep2.rds")
net_summary_rep2 <- community_level_dynamics(net_all_s_rep2 )
net_summary_rep2$v <- net_summary_rep2$naive_n-net_summary_rep2$primed_n
net_summary_rep2$e <- net_summary_rep2$naive_e-net_summary_rep2$primed_e
net_summary_rep2$origin <- "rep2"

net_communities_dt <- readRDS("3_network_dynamics/net_communities_dt.rds")
net_communities_dt$origin <- "Combined"

prot_gene_lookup <- readRDS("lookup_files/nodes_all_prot.rds")
cols <- c("ID", "prot_genes")
prot_gene_lookup <- prot_gene_lookup[, ..cols]

net_nodes <- lapply(net_summary_rep1$nodes, function(x) unlist(strsplit(x, ";")))
net_genes <- lapply(net_nodes, function(x) paste(unique(na.omit(prot_gene_lookup[x])$prot_genes), collapse = ";"))
net_summary_rep1$prot_genes <- unlist(net_genes)

c("5596_44", "5123_35")

net_nodes <- lapply(net_summary_rep2$nodes, function(x) unlist(strsplit(x, ";")))
net_genes <- lapply(net_nodes, function(x) paste(unique(na.omit(prot_gene_lookup[x])$prot_genes), collapse = ";"))
net_summary_rep2$prot_genes <- unlist(net_genes)

c("5267_44", "4817_32")

cols <- c("v", "e", "origin", ".id", "prot_genes")
data_plot <- rbind(net_summary_rep1[, ..cols], net_summary_rep2[, ..cols], net_communities_dt[, ..cols])



p <- ggplot(data_plot, aes(x = v, y = e, label=prot_genes)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  geom_point(data=data_plot[data_plot$.id %in% c("5596_44", "5123_35", "5267_44", "4817_32", "2609_3", "2387"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  facet_wrap(~origin) + 
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p
ggsave("3_network_dynamics/break_networks_non_merged_ve_deb2b_community_rep.pdf", 
       device="pdf", width = 14, height = 5)

pp <- ggplotly(p)
htmlwidgets::saveWidget(as_widget(pp), "break_networks_non_merged_ve_deb2b_community_rep.html")

ggplot(net_summary_rep1, aes(x = v, y = e)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  # geom_point(data=np_summary_deb2b_rep2[np_summary_deb2b_rep2$.id %in% c("5267", "4817"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

ggplot(net_summary_rep2, aes(x = v, y = e)) +
  geom_point(size=1, pch=21, fill="grey", alpha=0.5) +
  # geom_point(data=np_summary_deb2b_rep2[np_summary_deb2b_rep2$.id %in% c("5267", "4817"),], fill="red", size=2, pch=21) +
  xlab("Vertices (Naive-primed)") +
  ylab("Edges (Naive-primed)") +
  theme_pubr() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))


# Long range interaction plots -------------------------------------------------

distance_boxplot <- function(links_nodes){
  plot_dist <- log10(links_nodes$links$distance)
  
  naive_dist <- data.frame(distance=plot_dist[links_nodes$links$origin == "naive"])
  naive_dist$range <- "other"
  #remove trans
  naive_dist <- naive_dist[!is.na(naive_dist$distance),]
  
  naive_dist$range[order(naive_dist$distance, decreasing = TRUE)][1:1000] <- "long"
  naive_dist$range[order(naive_dist$distance, decreasing = TRUE)][(round(length(naive_dist$distance)/2)-99):(round(length(naive_dist$distance)/2)+500)] <- "mid"
  naive_dist$range[order(naive_dist$distance, decreasing = FALSE)][1:1000] <- "short"

  naive_dist$origin <- "naive"
  
  #remove other
  naive_dist <- naive_dist[naive_dist$range != "other",]
  
  
  primed_dist <- data.frame(distance=plot_dist[links_nodes$links$origin == "primed"])
  primed_dist$range <- "other"
  #remove trans
  primed_dist <- primed_dist[!is.na(primed_dist$distance),]

  primed_dist$range[order(primed_dist$distance, decreasing = TRUE)][1:1000] <- "long"
  primed_dist$range[order(primed_dist$distance, decreasing = TRUE)][(round(length(primed_dist$distance)/2)-499):(round(length(primed_dist$distance)/2)+500)] <- "mid"
  primed_dist$range[order(primed_dist$distance, decreasing = FALSE)][1:1000] <- "short"

  primed_dist$origin <- "primed"
  
  #remove other
  primed_dist <- primed_dist[primed_dist$range != "other",]
  
  boxplot_dist  <- rbind(naive_dist, primed_dist)
  
  
  box_p <- ggplot() +
    # geom_jitter(aes(x = origin, y = distance),data=boxplot_dist,shape = 20,size = 0.5,alpha = 0.1) +
    geom_sina(aes(y = distance, x = origin, color=origin),data=boxplot_dist, size=0.5, alpha=0.2)+
    geom_boxplot(aes(y = distance, x = origin),data=boxplot_dist, outlier.size=0, size = 0.2) +
    scale_y_continuous(breaks = c(4,5,6,7,8),labels = c('10 kb', '100 kb', '1 Mb', '10 Mb', '100 Mb')) +
    ylab("log10 interaction distance") + theme_bw()
  # box_p + facet_grid(~range)
  box_p
}

links_nodes_cat_col_coord_deb2b_rep1 <- readRDS("2_network_make/links_nodes_cat_col_coord_deb2b_rep1.rds")
links_nodes_cat_col_coord_deb2b_rep2 <- readRDS("2_network_make/links_nodes_cat_col_coord_deb2b_rep2.rds")

distance_boxplot(links_nodes_cat_col_coord_deb2b_rep1)
ggsave("5_network_distance/distance_boxplot_log10_rep1.pdf", device="pdf",
       width = 9, height = 9, units="cm")

distance_boxplot(links_nodes_cat_col_coord_deb2b_rep2)
ggsave("5_network_distance/distance_boxplot_log10_rep2.pdf", device="pdf",
       width = 9, height = 9, units="cm")

