library(igraph)
#rjson required

#'Replace dashes with dots
#'
remove_dash <- function(x) gsub("-", ".", x)


#'Make categories of naive primed SE and E to disply on network
#'
#'Notes: 
#'For each node (gene), an SE superseeds an Enhancer
#'
#'1 column
#'
#'B-backgroung
#'N-naive
#'P-primed
#'E-enhancer
#'SE-superenhancer
#'
#'N-B P-B = Bd
#'N-E P-B = Naive E
#'N-E P-E = Shared E
#'N-B P-E = Primed E
#'N-SE P-E = Naive SE
#'N-E P-SE = Primed SE
#'N-SE P-SE = Shared SE
#'N-B P-SE = Primed SE
#'N-SE P-B = Naive SE
#'
#'To make viewing easier
#'
#'nodes_all <- nodes_all[,c(setdiff(names(nodes_all),c("enhancer_naive", "enhancer_primed")), "enhancer_naive", "enhancer_primed")]
#'
categorise_enhancers <- function(links_nodes){

  nodes_all <- links_nodes$nodes
  nodes_all$enhancer_cat <- "B"
  
  
  nodes_all$enhancer_cat[nodes_all$enhancer_naive %in% c("S", "S_E", "S_S", "E_S") 
                         & nodes_all$enhancer_primed %in% c("S", "S_E", "S_S", "E_S")] <- "SE_shared"
  nodes_all$enhancer_cat[nodes_all$enhancer_naive %in% c("S", "S_E", "S_S", "E_S") 
                         & nodes_all$enhancer_primed %in% c("E", "E_E", "B")] <- "SE_naive"
  nodes_all$enhancer_cat[nodes_all$enhancer_naive %in% c("E", "E_E", "B") 
                         & nodes_all$enhancer_primed %in% c("S", "S_E", "S_S", "E_S")] <- "SE_primed"
  nodes_all$enhancer_cat[nodes_all$enhancer_naive %in% c("E", "E_E") 
                         & nodes_all$enhancer_primed %in% c("E", "E_E")] <- "E_shared"
  nodes_all$enhancer_cat[nodes_all$enhancer_naive %in% c("E", "E_E") 
                         & nodes_all$enhancer_primed %in% c("B")] <- "E_naive"
  nodes_all$enhancer_cat[nodes_all$enhancer_naive %in% c("B") 
                         & nodes_all$enhancer_primed %in% c("E", "E_E")] <- "E_primed"
  
  return(list(nodes=nodes_all, links=links_nodes$links))
}




#'simplify osn site overlap with hindiii fragments
#'hindiii containing both n and p are reclassified as mixed
#'
simplify_osn <- function(column_vector){
  
  #only process concatenated rows
  concat_rows <- grepl(";", column_vector)
  if(length(column_vector[grepl(";", column_vector)]) > 0){
    #returns a named vector
    osn_simp <- sapply(column_vector[grepl(";", column_vector)], 
                       function(x){
                         p_origin <- unique(unlist(strsplit(x, ";")))
                         if(length(p_origin) > 1){
                           return("m") #mixed
                         }else if(p_origin == "n"){
                           return("n")
                         }else if(p_origin == "p"){
                           return("p")
                         }else if(p_origin == "np"){
                           return("np")
                         }
                       })
    column_vector[concat_rows] <- unname(osn_simp)
    return(column_vector)
  }else{
    return(column_vector)
  }
}


#'
#'  interaction file from python script
#'  find duplicate b2b interactions and only keep highest score
#'  can lead to b2b to have different CHiCAGO score!
#'  naive b -> naive b 5.2 2.3
#'  primed b -> primed b 2.3 5.2
#'  primed b <- primed b 5.4 5.2
#'  naive b <- naive b 5.2 5.4
#'  so when mergining shared interactions take the max score to get 5.4 5.2
#'
remove_duplicate_b2b <- function(interact){

  
  original_nrow <- nrow(interact)
  #sort IDs by row to eliminate duplicate permutations
  b_oe_sorted <- t(apply(interact[c("b_ID", "oe_ID")], 1, sort))
  interact$b_ID_sorted <- b_oe_sorted[,1]
  interact$oe_ID_sorted <- b_oe_sorted[,2]
  #order by score to keep the highest score duplicate
  interact <- interact[order(-interact["score"]),]
  #mark duplicates
  interact$dups <- duplicated(interact[c("b_ID_sorted", "oe_ID_sorted")])
  
  #stop if not removing bait 2 bait only interactions
  stopifnot(interact$baited_b[interact$dups == "TRUE"] == interact$baited_oe[interact$dups == "TRUE"])
  
  interact <- interact[interact$dups == "FALSE",]
  
  print(paste(original_nrow - nrow(interact), "b2b duplicate interactions removed", 
              sep = " "))
  
  #drop added columns
  drop_columns <- c("b_ID_sorted", "oe_ID_sorted", "dups")
  interact_out <- interact[, !(names(interact) %in% drop_columns)]
  
  return(interact_out)
  
}


#'
#'  expected tables produced by network_visulisation.py
#'  include CHiCAGO scores for equivalent interactions ie score_naive score_primed
#'
make_link_node_dfs <- function(naive_nodes, primed_nodes, naive_interactions, primed_interactions, remove_b2b=FALSE){

  
  #remove dashes in gene names
  naive_intract <- dplyr::mutate_at(naive_interactions, c("b_genes", "oe_genes", 
                                                          "protein_b_genes", "protein_oe_genes"),
                                    remove_dash)
  
  primed_intract <- dplyr::mutate_at(primed_interactions, c("b_genes", "oe_genes", 
                                                            "protein_b_genes", "protein_oe_genes"),
                                     remove_dash)
  
  columns_sub <- c("b_ID","oe_ID","b_genes", "oe_genes","protein_b_genes", "protein_oe_genes",
                   "score", "score2", "distance")
  
  if(remove_b2b){
    naive_intract <- remove_duplicate_b2b(naive_intract)
    primed_intract <- remove_duplicate_b2b(primed_intract)
  }
  
  
  #merge naive and primed links into one
  naive_links <- naive_intract[,columns_sub]
  primed_links <- primed_intract[,columns_sub]
  
  #add bait2bait interaction column
  naive_links$b2b <- ifelse(naive_intract$baited_b == 1 & naive_intract$baited_oe == 1, 1, 0)
  primed_links$b2b <- ifelse(primed_intract$baited_b == 1 & primed_intract$baited_oe == 1, 1, 0)
  
  #add origin of interaction column
  naive_links$origin <- rep("naive", nrow(naive_links))
  primed_links$origin <- rep("primed", nrow(primed_links))
  
  links_all <- rbind(naive_links, primed_links)
  
  #make naive primed specific score columns
  links_all$score_naive <- ifelse(links_all$origin == "naive", links_all$score, links_all$score2)
  links_all$score_primed <- ifelse(links_all$origin == "primed", links_all$score, links_all$score2)
  #convert missing interactions to -1
  # links_all$score_naive[is.na(links_all$score_naive)] <- -1
  # links_all$score_primed[is.na(links_all$score_primed)] <- -1
  
  #subset all HindIII fragments to only those present within interaction data
  
  n <- data.frame(ID=c(naive_intract$b_ID, naive_intract$oe_ID))
  p <- data.frame(ID=c(primed_intract$b_ID, primed_intract$oe_ID))
  
  #remove duplicate nodes
  n <- unique(n)
  p <- unique(p)
  
  np_nodes <- merge(naive_nodes, primed_nodes, by=c("ID", "genes", "gtypes", "protein_genes", "OSN", "baited"), 
                    suffixes = c("_naive", "_primed"), all=TRUE)
  interacting_nodes <- merge(n, p, by=c("ID"), all=TRUE)
  #keep only nodes that interact
  nodes_all <- merge(interacting_nodes, np_nodes, by=c("ID"), all.x=TRUE)
  
  #remove - from gene names
  nodes_all <- dplyr::mutate_at(nodes_all, c("protein_genes", "genes"), remove_dash)
  #rename just keep compatibility with other functions
  names(nodes_all)[which(names(nodes_all) == "protein_genes")] <- "prot_genes"
  
  #Simplify OSN; resolve multiple HindIII OSN peak overlaps
  nodes_all$OSN <- simplify_osn(nodes_all$OSN)
  
  #change Inf value in insulation score (can't have them in network)
  in_naive_max <- max(nodes_all$in_naive[!is.infinite(nodes_all$in_naive)], na.rm = TRUE) + 1 
  in_primed_max <- max(nodes_all$in_primed[!is.infinite(nodes_all$in_primed)], na.rm = TRUE) + 1
  
  nodes_all$in_naive[is.infinite(nodes_all$in_naive)] <- in_naive_max
  nodes_all$in_primed[is.infinite(nodes_all$in_primed)] <- in_primed_max

  
  return(list(nodes=nodes_all, links=links_all))
}




#'  return the max fpkm value for a gene in a hindiii fragment
#'  fpkm values added from seqmonk onto the deseq results differential_analysis.R
#'
fpkm_max <- function(x, deseq_genes){

  if(!is.na(x)){
    
    # sapply(x, function(y) print(paste("test",y)))
    idx <- unlist(lapply(x, function(y){loc <- which(deseq_genes$gene_name %in% y)
    if(length(loc)==0){NA}else{loc}})) #resolve integer(0) returns
    #rough but works
    if(length(idx)>1){
      genes <- as.matrix(deseq_genes[idx,c("Naive_mean_fpkm_log2", "Primed_mean_fpkm_log2")])
    }else if(is.na(idx)){
      return(c(Naive_max_fpkm_log2=NA, Primed_max_fpkm_log2=NA))
    }else{
      genes <- as.matrix(deseq_genes[idx,c("Naive_mean_fpkm_log2", "Primed_mean_fpkm_log2")])
    }
    
    if(nrow(genes)> 1){
      naive_fpkm <- max(genes[,1], na.rm = TRUE)
      primed_fpkm <- max(genes[,2], na.rm = TRUE)
      return(c(Naive_max_fpkm_log2=naive_fpkm, Primed_max_fpkm_log2=primed_fpkm))
    }else{
      naive_fpkm <- unname(genes[,1])
      primed_fpkm <- unname(genes[,2])
      return(c(Naive_max_fpkm_log2=naive_fpkm, Primed_max_fpkm_log2=primed_fpkm))
    }
    
    
  }else{
    return(c(Naive_max_fpkm_log2=NA, Primed_max_fpkm_log2=NA))
  }
}


#'  Categorise expression into naive or primed (0 not significant, 1 primed significant,
#'  3 naive significant, 2 mixed)
#'  Categorise compartments (N - same in both, )
#'  Add FPKM max values per HindIII fragment
categorise_add_data <- function(links_nodes, deseq_genes){

  nodes_all <- links_nodes$nodes
  
  
  #Expression categorisation
  sig_de_genes <- deseq_genes[deseq_genes$padj < 0.05,]
  sig_naive <- sig_de_genes$gene_name[sig_de_genes$log2FoldChange < -0]
  sig_primed <- sig_de_genes$gene_name[sig_de_genes$log2FoldChange > 0]
  
  nodes_all$sig_prot_primed <- sapply(as.character(nodes_all$prot_genes), 
                                      function(x) sum((strsplit(x, ",")[[1]] %in% sig_primed)*1))
  nodes_all$sig_prot_naive <- sapply(as.character(nodes_all$prot_genes), 
                                     function(x) sum((strsplit(x, ",")[[1]] %in% sig_naive)*1))
  
  nodes_all$sign_prot <- 2 #everything else mixed
  nodes_all$sign_prot[nodes_all$sig_prot_primed == 0 & nodes_all$sig_prot_naive == 0] <- 0 #no sig
  nodes_all$sign_prot[nodes_all$sig_prot_primed > 0 & nodes_all$sig_prot_naive == 0] <- 1 #primed sig
  nodes_all$sign_prot[nodes_all$sig_prot_primed == 0 & nodes_all$sig_prot_naive > 0] <- 3 #naive sig
  
  
  #compartment categorisation
  nodes_all$compartments_n <- "N"
  nodes_all$compartments_p <- "N"
  n_a <- as.numeric(nodes_all$compartment_naive) > 0
  n_b <- as.numeric(nodes_all$compartment_naive) < 0
  p_a <- as.numeric(nodes_all$compartment_primed) > 0
  p_b <- as.numeric(nodes_all$compartment_primed) < 0
  nodes_all$compartments_n[n_a] <- "A"
  nodes_all$compartments_n[n_b] <- "B"
  nodes_all$compartments_p[p_a] <- "A"
  nodes_all$compartments_p[p_b] <- "B"
  
  #make AB BA and Same categories
  nodes_all$switched <- paste(nodes_all$compartments_n, nodes_all$compartments_p, sep = "")
  nodes_all$switched[(nodes_all$compartments_n == nodes_all$compartments_p)] <- "Same"
  
  #for compartments that switched get a difference value
  AB_switch <- nodes_all$switched == "AB"
  BA_switch <- nodes_all$switched == "BA"
  
  #filter low PC1 values
  AB_switched <- nodes_all$compartment_naive[AB_switch] > 5 &
                 nodes_all$compartment_primed[AB_switch] < -5
  BA_switched <- nodes_all$compartment_naive[BA_switch] < -5 &
                 nodes_all$compartment_primed[BA_switch] > 5
  
  nodes_all$switched_diff_np <- NA
  nodes_all$switched_diff_np[AB_switch][AB_switched] <- 
    nodes_all$compartment_naive[AB_switch][AB_switched]-nodes_all$compartment_primed[AB_switch][AB_switched]
  nodes_all$switched_diff_np[BA_switch][BA_switched] <- 
    nodes_all$compartment_naive[BA_switch][BA_switched]-nodes_all$compartment_primed[BA_switch][BA_switched]
  
  
    #add fpkm to nodes
    fpkm_list <- list()
    
    for(i in 1:nrow(nodes_all)){
      fpkm_list[[i]] <- fpkm_max(strsplit(as.character(nodes_all$prot_genes[i]), ","), deseq_genes)
    }
    
    
    rpkm_node_df <- plyr::ldply(fpkm_list, rbind)
    
    
    nodes_all$naive_log2FPKM <- rpkm_node_df$Naive_mean_fpkm_log2
    nodes_all$primed_log2FPKM <- rpkm_node_df$Primed_mean_fpkm_log2
  
  return(list(nodes=nodes_all, links=links_nodes$links))
}



get_spaced_colors <- function(n){
  #no R palletes contain enough colours for all TADs
  max_value <- 16581375 #255**3
  interval <- as.integer(max_value / n)
  hex_colours <- vector("character", n)
  x <- 1
  for(i in seq(0,max_value,interval)){
    hex_colours[x] <- sprintf("%06X",i)
    x <- x + 1
  }
  
  rgb_colours <- vector("character", n)
  for(i in seq_along(hex_colours)){
    rgb <- paste(strtoi(paste("0x", substr(hex_colours[i], 1, 2), sep="")),
                 strtoi(paste("0x", substr(hex_colours[i], 3, 4), sep="")),
                 strtoi(paste("0x", substr(hex_colours[i], 5, 6), sep="")),
                 sep = ",")
    
    # rgb_colours[i] <- paste("rgb(", rgb, ")", sep = "")
    rgb_colours[i] <- rgb
  }
  return(rgb_colours)
}


#'  Find equivalent tads and assign the same colour to all
add_same_tad_colour <- function(links_nodes){

  nodes_all <- links_nodes$nodes
  
  tad_id <- nodes_all[c("ID", "tad_primed", "tad_naive")]
  
  #create a list of equivalent tads
  tad_id_equival <- tad_id[!(is.na(tad_id$tad_primed) & is.na(tad_id$tad_naive)),]
  tad_id_equival$id <- paste(tad_id_equival$tad_primed, tad_id_equival$tad_naive, sep = "_")
  
  tad_id_eq <- tad_id_equival[!duplicated(tad_id_equival[c("tad_primed", "tad_naive")]),]
  tad_id_eq <- merge(plyr::count(tad_id_equival, c("id")), tad_id_eq, by="id", all = TRUE)
  
  #deduplicate
  #dup on one, freq then check if the other side has a dup, if not mark for unique numbering
  #create a colour table straight away
  # TAD ID              TAD ID freq HindIII ID  Primed TAD  Naive TAD
  # chr10_250_chr10_263 51          chr10_3552  chr10_250   chr10_263 <- colour by this
  # chr10_250_NA        11          chr10_3451  chr10_250   NA
  
  p_dup <- tad_id_eq$tad_primed[duplicated(tad_id_eq$tad_primed)]
  n_dup <- tad_id_eq$tad_naive[duplicated(tad_id_eq$tad_naive)]
  
  # test <- tad_id_eq[tad_id_eq$tad_primed %in% p_dup & tad_id_eq$tad_naive %in% n_dup,]
  
  for(i in seq_along(p_dup)){
    dup_idx <- which(tad_id_eq$tad_primed %in% p_dup[i])
    
    min_idx <- which.min(tad_id_eq$freq[dup_idx])
    
    tad_id_eq$tad_primed[dup_idx[min_idx]] <- NA
  }
  
  for(i in seq_along(n_dup)){
    dup_idx <- which(tad_id_eq$tad_naive %in% n_dup[i])
    
    min_idx <- which.min(tad_id_eq$freq[dup_idx])
    
    tad_id_eq$tad_naive[dup_idx[min_idx]] <- NA
  }
  
  tad_id_eq_col <- tad_id_eq[!(is.na(tad_id_eq$tad_primed) & is.na(tad_id_eq$tad_naive)),]
  tad_id_eq_col[is.na(tad_id_eq_col)] <- "empty"
  
  colours <- get_spaced_colors(nrow(tad_id_eq_col))
  
  tad_id_eq_col$colours <- colours[-1]
  
  tad_id_eq_col_n <- data.table::as.data.table(tad_id_eq_col)
  data.table::setkey(tad_id_eq_col_n, key="tad_naive")
  
  tad_id_eq_col_p <- data.table::as.data.table(tad_id_eq_col)
  data.table::setkey(tad_id_eq_col_p, key="tad_primed")
  
  
  #add colour to the nodes_all_df
  tads_n_find <- nodes_all$tad_naive[!is.na(nodes_all$tad_naive)] 
  tads_p_find <- nodes_all$tad_primed[!is.na(nodes_all$tad_primed)]
  
  
  nodes_all$color_tad_naive <- "255,255,255"
  nodes_all$color_tad_naive[!is.na(nodes_all$tad_naive)] <- tad_id_eq_col_n[.(tads_n_find), nomatch = 0L]$colours
  nodes_all$color_tad_primed <- "255,255,255"
  nodes_all$color_tad_primed[!is.na(nodes_all$tad_primed)] <- tad_id_eq_col_p[.(tads_p_find), nomatch = 0L]$colours
  
  
  # tads_n_find <- V(net_all_trans_s)$tad_naive[!is.na(V(net_all_trans_s)$tad_naive)]
  # tads_p_find <- V(net_all_trans_s)$tad_primed[!is.na(V(net_all_trans_s)$tad_primed)]
  # 
  # V(net_all_trans_s)$color_tad_naive <- "255,255,255"
  # V(net_all_trans_s)$color_tad_naive[!is.na(V(net_all_trans_s)$tad_naive)]  <- tad_id_eq_col_n[.(tads_n_find), nomatch = 0L]$colours
  # 
  # V(net_all_trans_s)$color_tad_primed <- "255,255,255"
  # V(net_all_trans_s)$color_tad_primed[!is.na(V(net_all_trans_s)$tad_primed)]  <- tad_id_eq_col_p[.(tads_p_find), nomatch = 0L]$colours
  
  # write_graph(net_all_trans_s, "/media/chovanec/My_Passport/CHiC_naive_primed/network/network_analysis/most_diff_networks/net_all_trans_s_cluster_reduced_tadcolor.graphml", 
  #             format = "graphml")
  
  return(list(nodes=nodes_all, links=links_nodes$links))
}


#'#Make an Igraph network (retaining trans interactions)
#'Add IDs for all the subnetworks so it can be linked back to dynimics analysis
#'
make_network <- function(links_nodes, np_summary_out, directed=TRUE, save_nontrans_net=NA){
  
  links_all_trans <- links_nodes$links
  links_all_trans$count <- 1
  #mark trans as 1 so distance can be used a weight
  links_all_trans$distance[is.na(links_all_trans$distance)] <- 1
  #drop score and score2 columns and only use score_naive score_primed
  links_all_trans <- links_all_trans[,!(names(links_all_trans) %in% c("score", "score2"))]
  
  links_all <- links_nodes$links[!is.na(links_nodes$links$distance),]
  links_all <- links_all[,!(names(links_all) %in% c("score", "score2"))]
  links_all$count <- 1
 
  #keep max CHiCAGO score because of b2b merging
  edge_simplify <- list(score_naive="max", score_primed="max", oe_genes="first", 
                        b_genes="first", protein_oe_genes="first", 
                        protein_b_genes="first", distance="first", count="sum", b2b="max",
                        origin=function(x) paste(x, collapse="_"))
  
  
  #make trans network
  net_all_trans <- igraph::graph_from_data_frame(links_all_trans, vertices = links_nodes$nodes, 
                                                 directed = directed)
  net_all_trans_s <- igraph::simplify(net_all_trans, edge.attr.comb = edge_simplify)
  
  #make network without trans for analysis of subnetworks
  net_all <- igraph::graph_from_data_frame(links_all, vertices = links_nodes$nodes, 
                                           directed = directed)
  net_all_s <- igraph::simplify(net_all, edge.attr.comb = edge_simplify)
  
  # decompose the graph
  all_subgraphs <- decompose.graph(net_all_s)
  
  #don't need to run this if file already exists
  if(file.exists(np_summary_out)){
    print(paste("File already exists and will be read from:", np_summary_out))
    # net_primed <- NA
    # net_naive <- NA
  }else{
    net_primed <- analyse_net(all_subgraphs, links_nodes$nodes, "naive") #deleting naive
    net_naive <- analyse_net(all_subgraphs, links_nodes$nodes, "primed") #deleting primed
  }
  
  np_summary <- summarise_net_analysis(net_naive, net_primed, np_summary_out)
  
  net_out <- add_subnet(net_all_trans_s, np_summary, all_subgraphs)
  net_out_notrans <- add_subnet(net_all_s, np_summary, all_subgraphs)
  
  #Make distance edge weights
  E(net_out)$weight <- E(net_out)$distance
  
  if(!is.na(save_nontrans_net)){
    saveRDS(net_out_notrans, save_nontrans_net)
  }
  
  return(net_out)
}


#' incorporate coordinates produced by gephi layout into igraph network
#' This requires that the network has been produced imported into gephi and
#' then its coordinates exported in json format. This json file can the be used here
add_layout_coordinate <- function(links_nodes, gephi_json){
  
  graph_json <- rjson::fromJSON(file=gephi_json)
  
  graph_xy <- list()
  for(i in seq_along(graph_json$nodes)){
    graph_xy[[i]] <- data.frame(name=graph_json$nodes[[i]]$attributes$name,
                                x=graph_json$nodes[[i]]$x,
                                y=graph_json$nodes[[i]]$y)
    
  }
  
  graph_xy_df <- plyr::ldply(graph_xy, rbind)
  names(graph_xy_df)[1] <- "ID"
  nodes_all <- merge(links_nodes$nodes, graph_xy_df, by="ID", all.x=TRUE)
  
  return(list(nodes=nodes_all, links=links_nodes$links))
}


