library(readr)
require(igraph)
library(tcltk2)
library(plyr)
library(dplyr)
library(data.table)
library(Sushi)

source("network_dynamics_functions.R")
source("network_make_functions.R")

#Load RNA-seq data -------------------------------------------------------------

de_genes <- read_delim("2_network_make/de_genes_takashima_GRCh38.87_anno_opposing_strand_prot_genes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)



#Load interaction data ---------------------------------------------------------

naive_nodes <- read_delim("2_network_make/naive_nodes_wchip_20181218.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

primed_nodes <- read_delim("2_network_make/primed_nodes_wchip_20181218.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

naive_intract <- read_delim("2_network_make/naive_intract_wchip_20181218.txt", 
                            "\t", escape_double = FALSE, col_names = TRUE, 
                            trim_ws = TRUE)

primed_intract <- read_delim("2_network_make/primed_intract_wchip_20181218.txt", 
                             "\t", escape_double = FALSE, col_names = TRUE, 
                             trim_ws = TRUE)

#with 3 score cutoff
naive_intract_3 <- read_delim("2_network_make/naive_3_intract_wchip_20181218.txt", 
                              "\t", escape_double = FALSE, col_names = TRUE, 
                              trim_ws = TRUE)

primed_intract_3 <- read_delim("2_network_make/primed_3_intract_wchip_20181218.txt", 
                               "\t", escape_double = FALSE, col_names = TRUE, 
                               trim_ws = TRUE)



#Make network without b2b duplicates for quantification of No. interactions ----

links_nodes_deb2b <- make_link_node_dfs(naive_nodes, primed_nodes, naive_intract, primed_intract, remove_b2b = TRUE)

#only protein coding genes
links_nodes_cat_deb2b <- categorise_add_data(links_nodes_deb2b, de_genes)

#Colour each tads with the same colour between naive and primed
links_nodes_cat_col_deb2b <- add_same_tad_colour(links_nodes_cat_deb2b)

#Enhancer categories
links_nodes_cal_col_bed2b_enh <- categorise_enhancers(links_nodes_cat_col_deb2b)

#Add Gephi layout coordinates (to create reproducible layouts, this has to be exported from gephi)
links_nodes_cat_col_coord_deb2b <- add_layout_coordinate(links_nodes_cal_col_bed2b_enh, 
                                                         "2_network_make/network_expression/data.json")

#Make network 
net_out_deb2b <- make_network(links_nodes_cat_col_coord_deb2b, "2_network_make/np_summary_deb2b_20190829.rds",
                              directed = FALSE, save_nontrans_net = "2_network_make/net_all_s_20190829.rds")

igraph::write_graph(net_out_deb2b, "2_network_make/net_all_trans_s_deb2b_20190829.graphml", 
                    format = "graphml")

saveRDS(links_nodes_cat_col_coord_deb2b,"2_network_make/links_nodes_cat_col_coord_deb2b_20190829.rds")


#Make 3 cutoff network without b2b ---------------------------------------------


links_nodes_3_deb2b <- make_link_node_dfs(naive_nodes, primed_nodes, naive_intract_3, primed_intract_3, remove_b2b = TRUE)

#only protein coding genes
links_nodes_3_cat_deb2b <- categorise_add_data(links_nodes_3_deb2b, de_genes)

#Colour each tads with the same colour between naive and primed
links_nodes_3_cat_col_deb2b <- add_same_tad_colour(links_nodes_3_cat_deb2b)

#Make network 
net_out_3_deb2b <- make_network(links_nodes_3_cat_col_deb2b, "2_network_make/np_summary_3_deb2b_20190829.rds",
                                directed = FALSE, save_nontrans_net = "2_network_make/net_all_3_s_20190829.rds")

igraph::write_graph(net_out_3_deb2b, "2_network_make/net_all_3_trans_s_deb2b_20190829.graphml", 
                    format = "graphml")

saveRDS(links_nodes_3_cat_col_deb2b,"2_network_make/links_nodes_3_cat_col_coord_deb2b_20190829.rds")

