library(readr)
require(igraph)
library(tcltk2)
library(plyr)
library(dplyr)
library(data.table)
library(Sushi)
library(rjson)
library(RColorBrewer)

source("network_dynamics_functions.R")
source("network_make_functions.R")

#Load RNA-seq data -------------------------------------------------------------

de_genes <- read_delim("1_DESeq2_gene_expression/de_genes_takashima_GRCh38.87_anno_opposing_strand_prot_genes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)



#Load interaction data ---------------------------------------------------------

# New ChromHMM
naive_nodes <- read_delim("0_network_data_preparation/naive_nodes_wchip_20200911.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

primed_nodes <- read_delim("0_network_data_preparation/primed_nodes_wchip_20200911.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

naive_intract <- read_delim('0_network_data_preparation/naive_intract_wchip_20200911.txt', 
                            "\t", escape_double = FALSE, col_names = TRUE, 
                            trim_ws = TRUE)

primed_intract <- read_delim('0_network_data_preparation/primed_intract_wchip_20200911.txt', 
                             "\t", escape_double = FALSE, col_names = TRUE, 
                             trim_ws = TRUE)

# with 3 score cutoff
naive_intract_3 <- read_delim("0_network_data_preparation/naive_3_intract_wchip_20200911.txt", 
                              "\t", escape_double = FALSE, col_names = TRUE, 
                              trim_ws = TRUE)

primed_intract_3 <- read_delim("0_network_data_preparation/primed_3_intract_wchip_20200911.txt", 
                               "\t", escape_double = FALSE, col_names = TRUE, 
                               trim_ws = TRUE)

# Custom chromosome colour palette ---------------------------------------------

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,24), col=col_vector[c(5,6,10,11,13,14,15,16,17,18,19,20,21,22,25,26,29,49,54,55,56,57,71,74)])

custom_chrom_col <- apply(col2rgb(col_vector[c(5,6,10,11,13,14,15,16,17,18,19,20,21,22,25,26,29,49,54,55,56,57,71,74)]), 2, function(x) paste(x, collapse = ","))


# Make network without b2b duplicates for quantification of No. interactions ----

links_nodes_deb2b <- make_link_node_dfs(naive_nodes, primed_nodes, naive_intract, primed_intract, remove_b2b = TRUE)

# only protein coding genes
links_nodes_cat_deb2b <- categorise_add_data(links_nodes_deb2b, de_genes)

# Colour each tads with the same colour between naive and primed
links_nodes_cat_col_deb2b <- add_same_tad_colour(links_nodes_cat_deb2b)

# Enhancer categories
links_nodes_cal_col_bed2b_enh <- categorise_enhancers(links_nodes_cat_col_deb2b)

# Add Gephi layout coordinates (to create reproducible layouts, this has to be exported from gephi)
links_nodes_cat_col_coord_deb2b <- add_layout_coordinate(links_nodes_cal_col_bed2b_enh, 
                                                         "2_network_make/data.json")

# Make network 
net_out_deb2b <- make_network(links_nodes_cat_col_coord_deb2b, "2_network_make/np_summary_deb2b_20200911.rds",
                              directed = FALSE, save_nontrans_net = "2_network_make/net_all_s_20200911.rds",
                              chrom_colours=custom_chrom_col)

igraph::write_graph(net_out_deb2b, "2_network_make/net_all_trans_s_deb2b_20200911.graphml", 
                    format = "graphml")

saveRDS(links_nodes_cat_col_coord_deb2b,"2_network_make/links_nodes_cat_col_coord_deb2b_20200911.rds")


# Make 3 cutoff network without b2b --------------------------------------------


links_nodes_3_deb2b <- make_link_node_dfs(naive_nodes, primed_nodes, naive_intract_3, primed_intract_3, remove_b2b = TRUE)

# only protein coding genes
links_nodes_3_cat_deb2b <- categorise_add_data(links_nodes_3_deb2b, de_genes)

# Colour each tads with the same colour between naive and primed
links_nodes_3_cat_col_deb2b <- add_same_tad_colour(links_nodes_3_cat_deb2b)

# Make network 
net_out_3_deb2b <- make_network(links_nodes_3_cat_col_deb2b, "2_network_make/np_summary_3_deb2b_20200911.rds",
                                directed = FALSE, save_nontrans_net = "2_network_make/net_all_3_s_20200911.rds")

igraph::write_graph(net_out_3_deb2b, "2_network_make/net_all_3_trans_s_deb2b_20200911.graphml", 
                    format = "graphml")

saveRDS(links_nodes_3_cat_col_deb2b,"2_network_make/links_nodes_3_cat_col_coord_deb2b_20200911.rds")


# Make supplemental table 1 ----------------------------------------------------

network_table <- as_long_data_frame(net_out_deb2b)
hindiii_lookup <- read_delim("7_enhancers/hindiii_lookup.tab", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
names(hindiii_lookup) <- c("chrom", "start", "end", "ID")

hindiii_lookup$coord <- apply(hindiii_lookup, 1, function(x) paste(x["chrom"], ":", x["start"], "-", x["end"], collapse=""))
hindiii_lookup$coord <- gsub(" ", "", hindiii_lookup$coord)

data.table::setnames(network_table, "from_name", "ID")
network_table_2 <- merge(network_table, hindiii_lookup[c("ID", "coord")], by="ID")
data.table::setnames(network_table_2, "coord", "from_coord")
data.table::setnames(network_table_2, "ID", "from_name")
data.table::setnames(network_table_2, "to_name", "ID")
network_table_2 <- merge(network_table_2, hindiii_lookup[c("ID", "coord")], by="ID")
data.table::setnames(network_table_2, "coord", "to_coord")
data.table::setnames(network_table_2, "ID", "to_name")
write.table(network_table_2, "2_network_make/S1_table.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
