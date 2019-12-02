
#' the get_interaction_counts nodes will be baited genes
#' summaries which enhancers are contacting each gene - ID, type
#' get a summary of the chromatin state of all the interacting fragments
#' 
#' 
node_enhancer_summary <- function(root_node, nodes_dt, lookup_dt, score_diff=2){

  
  #TODO how do you export an s4 class for foreach?
  source("network_enhancer_class.R")
  
  lvls <- c("Active", "Polycomb Repressed", "Bivalent", 
            "Heterochromatin Repressed", "H3K4me1", "Mixed",
            "Background", "Unknown")
  
  if(nrow(nodes_dt) == 0){
    #make enhancer out table
    #everything should be NA if not interactions present
    out <- new("Interactions", node = root_node,
               interactions = setNames(data.table(matrix(nrow = 0, ncol = 7)), 
                                       c("from_name", "to_name", "origin", "b2b", 
                                         "score_naive", "score_primed", "new_origin")),
               naive_enhancers = setNames(data.table(matrix(nrow = 0, ncol = 6)), 
                                          c("ID", "chromhmm_naive", "name_naive", 
                                            "score_naive", "type_naive", "OSN_naive")),
               primed_enhancers = setNames(data.table(matrix(nrow = 0, ncol = 6)), 
                                           c("ID", "chromhmm_primed", "name_primed", 
                                             "score_primed", "type_primed", "OSN_primed")),
               naive_states = table(factor(NA, levels = lvls)),
               primed_states = table(factor(NA, levels = lvls))
               )
               
  }else{
    #seperate lookup table into naive and primed
    naive_cols <- c("ID", names(lookup_dt)[grepl("naive", names(lookup_dt))])
    primed_cols <- c("ID", names(lookup_dt)[grepl("primed", names(lookup_dt))])
    
    naive_lookup_dt <- lookup_dt[,..naive_cols]
    primed_lookup_dt <- lookup_dt[,..primed_cols]
    
    #when assigning origin also take score into account and make a 
    #new origin column 
    #if one enhancer below 5 score, record the difference in score
    #reasign origin only if the score difference is < 2 (or score_diff)
    less_than_five <- nodes_dt$score_naive < 5 | nodes_dt$score_primed < 5
    diff_less <- abs(nodes_dt$score_naive - nodes_dt$score_primed) < score_diff
    nodes_dt$new_origin <- nodes_dt$origin
    nodes_dt$new_origin[less_than_five & diff_less] <- "naive_primed"
    
    #seperate into naive and primed
    primed_orig <- nodes_dt$new_origin %in% c("primed", "naive_primed")
    intr_frag_primed <- c(nodes_dt$from_name[primed_orig], 
                          nodes_dt$to_name[primed_orig])
    
    naive_orig <- nodes_dt$new_origin %in% c("naive", "naive_primed")
    intr_frag_naive <- c(nodes_dt$from_name[naive_orig],
                         nodes_dt$to_name[naive_orig])
    
    #remove self node, only look for other end
    #will not tell you if the gene is itself an enhancer
    intr_frag_primed <- intr_frag_primed[!(intr_frag_primed %in% root_node)]
    intr_frag_naive <- intr_frag_naive[!(intr_frag_naive %in% root_node)]
    
    #lookup enhancers, OSN and chromhmm states
    primed_data <- primed_lookup_dt[intr_frag_primed]
    naive_data <- naive_lookup_dt[intr_frag_naive]
    
    
    #make enhancer out table
    out <- new("Interactions", node = root_node,
               interactions = nodes_dt,
               naive_enhancers = naive_data[!is.na(naive_data$name_naive)],
               primed_enhancers = primed_data[!is.na(primed_data$name_primed)],
               naive_states = table(factor(naive_data$chromhmm_naive, levels = lvls)),
               primed_states = table(factor(primed_data$chromhmm_primed, levels = lvls))
    )
  }
  return(out)
  
}




#'
#' output a summary list of expression values from a data.table lookup
#'
node_expression_summary <- function(search_nodes, lookup_dt, unique_name=NA){
  
  print(search_nodes)
  if(length(search_nodes) == 0){
    out_list <- list(intr_frag=NA, prot_genes=NA, Naive_mean_fpkm_log2=NA,
                     Primed_mean_fpkm_log2=NA, log2FoldChange=NA, Naive_expr_mean=NA,
                     Primed_expr_mean=NA, log2FoldChange_mean=NA)
    
  }else{
    
    out_list <- list()
    #only interested in these columns from gene_lookup_DT
    intr_frag_genes <- lookup_dt[search_nodes, c("prot_genes", 
                                                 "Naive_mean_fpkm_log2", 
                                                 "Primed_mean_fpkm_log2",
                                                 "log2FoldChange")]
    
    out_list[["intr_frag"]] <- paste(search_nodes, collapse = ";")
    out_list[["prot_genes"]] <- paste(intr_frag_genes$prot_genes, collapse = ";")
    out_list[["Naive_mean_fpkm_log2"]] <- paste(intr_frag_genes$Naive_mean_fpkm_log2, 
                                                collapse = ";")
    out_list[["Primed_mean_fpkm_log2"]] <- paste(intr_frag_genes$Primed_mean_fpkm_log2, 
                                                 collapse = ";")
    out_list[["log2FoldChange"]] <- paste(intr_frag_genes$log2FoldChange, 
                                          collapse = ";")
    out_list[["Naive_expr_mean"]] <- mean(intr_frag_genes$Naive_mean_fpkm_log2,
                                          na.rm = TRUE)
    out_list[["Primed_expr_mean"]] <- mean(intr_frag_genes$Primed_mean_fpkm_log2,
                                           na.rm = TRUE)
    out_list[["log2FoldChange_mean"]] <- mean(intr_frag_genes$log2FoldChange,
                                              na.rm = TRUE)
  }
  
  if(!is.na(unique_name)){
    names(out_list) <- paste(names(out_list), unique_name, sep = "_")
  }
  
  return(out_list)
}




#' Get count of interactions from a vector of nodes
#' 
#' @param node_IDs A vector or list of nodes
#' @param network An igraph network from which interaction data will be extracted
#' @param lookup_dt A data.table with lookup data for expression or enhancers
#' @param p Number of threads to use with foreach
#' @param type Add information about enhancers/OSN or about expression
#' @param score_diff If CHiCAGO score between naive primed is less then, reasign as shared
#' 
#' @return A data.frame with all values aggregated on node_IDs
#' @examples 
#' get_interaction_counts(gene_nodes_list, net_all_s, chromhmm_enh_OSN_lookup, p=7, count_without_b2b = FALSE, type = "enhancer", score_diff = 2, gene_lookup_dt)
#' get_interaction_counts(OSN_ID_primed, net_out_deb2b, nodes_all_prot, p=7, count_without_b2b = TRUE, type = "expression")
get_interaction_counts <- function(node_IDs, network, lookup_dt, p=7, 
                                   count_without_b2b=TRUE, 
                                   type = c("expression", "enhancer"),
                                   score_diff=2,
                                   gene_node_ID_dt=NULL){
  #find what a set of nodes is interacting with
  #lookup_dt = either gene expression lookup for node_expression_summary
  #or enhancer lookup
  #all genes in lookup_dt are baited protein coding genes (don't care about oe genes)
  
  #lookup_dt - for expression this will be the de_genes dt
  #          - for enhancers this will be the chromhmm and enhancer dt
  
  if(all(c("expression", "enhancer") %in% type)){
    stop("Specify type - either expression or enhancer")
  }
  
  #check the lookup table has expected columns
  if(type == "expression"){
    expected_cols <- c("prot_genes", "ID", "log2FoldChange", "padj", 
                       "Naive_mean_fpkm_log2", "Primed_mean_fpkm_log2")
    if(!all(expected_cols %in% names(lookup_dt))){
      missing_cols <- expected_cols[!(expected_cols %in% names(lookup_dt))]
      stop(paste("The lookup table is missing the following columns:", 
                 paste(missing_cols, collapse = ", "), sep = " "))
    }
  }
 
  if(type == "enhancer"){
    expected_cols <- c("ID", "chromhmm_naive", "chromhmm_primed", "name_naive", 
                       "score_naive", "type_naive", "OSN_naive", "name_primed", 
                       "score_primed", "type_primed","OSN_primed")
    if(!all(expected_cols %in% names(lookup_dt))){
      missing_cols <- expected_cols[!(expected_cols %in% names(lookup_dt))]
      stop(paste("The lookup table is missing the following columns:", 
                 paste(missing_cols, collapse = ", "), sep = " "))
    }
  } 
  
  #Set up search data.table from simplified network
  net_data_all <- igraph::as_long_data_frame(network)
  setDT(net_data_all, key = c("from_name", "to_name"))
  #reduce to only needed columns
  cols <- c("from_name", "to_name", "origin", "b2b", "score_naive", "score_primed")
  net_data <- net_data_all[, ..cols]

  #Set up parallel end
  cl <- makeCluster(p)
  registerDoSNOW(cl)

  #Set up progress bar
  pb <- txtProgressBar(max = length(node_IDs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  #return a list to keep it more flexible
  net_summary <-
    foreach(r=1:length(node_IDs), #.combine = "rbind",
            #.final = function(x) setNames(x, as.character(node_IDs)),
            .options.snow = opts, .packages = c("igraph", "data.table"),
            .inorder=TRUE, .export = c("node_expression_summary",
                                       "node_enhancer_summary")) %dopar% {

            # for(r in 1:length(node_IDs)){
              
              # #doing this from networks is very slow
              # #interaction ends
              # node_edges <- E(network)[from(V(network)[node_IDs[[r]]])]
              # int_ends <- ends(network, node_edges)
              # #origin of interaction
              # int_ends_origin <- cbind(int_ends, node_edges$origin, node_edges$b2b)
      
              #origin of interaction
              int_ends_origin <- rbind(net_data[CJ(unique(from_name), node_IDs[[r]]), 
                                                nomatch = 0L],
                                       net_data[node_IDs[[r]], nomatch = 0L])
              
              
              orig_freq <- table(factor(int_ends_origin$origin,
                                        levels = c("naive", "primed", "naive_primed")))
  
              if(count_without_b2b){
                #counts of only b to oe interactions (removing b2b interactions)
                int_b_oe <- int_ends_origin$b2b==0

                orig_freq_wo_b2b <- table(factor(int_ends_origin$origin[int_b_oe],
                                                 levels = c("naive", "primed", 
                                                            "naive_primed")))

                # out_df <- cbind(data.frame(as.list(orig_freq)),
                #                 data.frame(as.list(orig_freq_wo_b2b)))

                orig_freq_wo_b2b_list <- as.list(orig_freq_wo_b2b)
                names(orig_freq_wo_b2b_list) <- c("naive_wo_b2b", "primed_wo_b2b", 
                                                  "naive_primed_wo_b2b")

                freq_list <- c(as.list(orig_freq), orig_freq_wo_b2b_list)
              }else{

                freq_list <- as.list(orig_freq)
              }


              freq_list[["ID"]] <- paste(node_IDs[[r]], collapse = ";")


              #get gene and expression for each interacting hindiii
              if(type == "expression"){
                
                #check score and reassign origin based on difference
                #when assigning origin also take score into account and make a 
                #new origin column 
                #if one enhancer below 5 score, record the difference in score
                #reasign origin only if the score difference is < 2 (or score_diff)
                less_than_five <- int_ends_origin$score_naive < 5 | int_ends_origin$score_primed < 5
                diff_less <- abs(int_ends_origin$score_naive - int_ends_origin$score_primed) < score_diff
                int_ends_origin$new_origin <- int_ends_origin$origin
                int_ends_origin$new_origin[less_than_five & diff_less] <- "naive_primed"

                primed_orig <- int_ends_origin$new_origin %in% c("primed", "naive_primed")
                intr_frag_primed <- c(int_ends_origin$from_name[primed_orig],
                                      int_ends_origin$to_name[primed_orig])
                naive_orig <- int_ends_origin$new_origin %in% c("naive", "naive_primed")
                intr_frag_naive <- c(int_ends_origin$from_name[naive_orig],
                                     int_ends_origin$to_name[naive_orig])

                #remove self
                intr_frag_primed <- intr_frag_primed[!(intr_frag_primed %in% node_IDs[[r]])]
                intr_frag_naive <- intr_frag_naive[!(intr_frag_naive %in% node_IDs[[r]])]

                naive_list <- data.frame(node_expression_summary(intr_frag_naive,
                                                                 lookup_dt,
                                                                 "naive"))
                primed_list <- data.frame(node_expression_summary(intr_frag_primed,
                                                                  lookup_dt,
                                                                  "primed"))

                out_data <- cbind(data.frame(freq_list), data.frame(naive_list), 
                                  data.frame(primed_list))

              }else if(type == "enhancer"){
                #produces an Interaction class object
                enh_interactions <- node_enhancer_summary(node_IDs[[r]], int_ends_origin,
                                                          lookup_dt, score_diff)
                
                enh_sm_naive <- enhsummary(enh_interactions, "naive", gene_node_ID_dt)
                enh_sm_primed <- enhsummary(enh_interactions, "primed", gene_node_ID_dt)
                stat_naive <- statecounts(enh_interactions, "naive", gene_node_ID_dt)
                stat_primed <- statecounts(enh_interactions, "primed", gene_node_ID_dt)
                enh_cnt_naive <- enhcounts(enh_interactions, "naive", gene_node_ID_dt)
                enh_cnt_primed <- enhcounts(enh_interactions, "primed", gene_node_ID_dt)

                #don't want the chromhmm_naive and chromhmm_primed columns
                enh_sm_naive[, chromhmm_naive := NULL]
                enh_sm_primed[, chromhmm_primed := NULL]

                if(is.null(gene_node_ID_dt)){
                  out_enh <- plyr::join_all(list(enh_sm_naive, enh_sm_primed,
                                                 stat_naive, stat_primed,
                                                 enh_cnt_naive, enh_cnt_primed),
                                           by="int_node")
                }else{
                  out_enh <- plyr::join_all(list(enh_sm_naive, enh_sm_primed, 
                                                 stat_naive, stat_primed,
                                                 enh_cnt_naive, enh_cnt_primed),
                                           by="prot_genes")
                }
                
                #need to include new_origin and use those counts
                freq_list <- as.list(table(factor(interactions(enh_interactions)$new_origin,
                                          levels = c("naive", "primed", "naive_primed"))))
                freq_list[["ID"]] <- paste(node_IDs[[r]], collapse = ";")
                
                out_data <- cbind(data.frame(freq_list), out_enh)
              }else if(type == "raw"){
                #return S4 Interaction class for each node_ID
                #Useful for circos plots of each gene!
                out_data <- node_enhancer_summary(node_IDs[[r]], int_ends_origin,
                                                          lookup_dt, score_diff)
                #include gene name for node
                if(!is.null(gene_node_ID_dt)){
                  genes(out_data) <- paste(unique(gene_node_ID_dt[out_data@node]$prot_genes),
                                           collapse = ";")
                }
              }

              # print(out_data)
              return(out_data)

            }
  close(pb)
  stopCluster(cl)

  if(type %in% c("expression", "enhancer")){
    out <- do.call(rbind.data.frame, net_summary)
    #does something funky, converting string to 1 out <- plyr::ldply(net_summary, rbind)
  }else if(type == "raw"){
    
    #get a vector of gene names to use as list names
    gene_names <- vector("character")
    for(i in seq_along(net_summary)){
      gene_names[i] <- genes(net_summary[[i]])
    }
    
    names(net_summary) <- gene_names
    out <- net_summary
    
  }

  return(out)
}
