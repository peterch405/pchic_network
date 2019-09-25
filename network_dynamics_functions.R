library(igraph)
library(tcltk2)
library(Sushi)

#'Isolate naive or primed network from merged by edge origin attribute
#'
#'@param ig_network igraph network object
#'@param keep_type the cell type to keep
#'@param delete_shared interactions shared between cell types will be kept by default
#'
#'
#'
isolate_cell_network <- function(ig_network, keep_type=c("naive", "primed"), delete_shared=FALSE){
  
  if(keep_type == "naive"){
    type <- "primed"
  }else if(keep_type == "primed"){
    type <- "naive"
  }else{
    stop("keep_type can only be 'naive' or 'primed'")
  }
  
  isolated_net <- delete.edges(ig_network,E(ig_network)[E(ig_network)$origin == type]) #does not delete shared edges
  # plot(s_sub_graph, vertex.label=NA, vertex.size=4)
  # components(primed_sub_graph)$membership
  if(delete_shared){
    isolated_net <- delete.edges(isolated_net,E(isolated_net)[E(isolated_net)$origin == "naive_primed"])
  }
  
  #delete vertices with no connections
  isolated_net_out <- delete.vertices(isolated_net,which(degree(isolated_net)<1))
  
  return(isolated_net_out)
}







#functions from network dynamics for subnetwork analysis -----------------------

#'Sub-network analysis
#'
#'
#'
#'
#'
analyse_net <- function(all_subgraphs, nodes_all, type, del_shared=FALSE){
  baited <- nodes_all$ID[nodes_all$baited]
  #initialise progress bar
  pb <- tkProgressBar("Progress bar", "Some information in %",
                      0, length(all_subgraphs), 0, width=300)
  
  
  all_net_list <- list()
  for(sg in seq_along(all_subgraphs)){
    sub_graph <- all_subgraphs[[sg]]
    
    #which nodes are lost when remove naive or primed interactions (only have those genes listed)
    s_sub_graph <- delete.edges(sub_graph,E(sub_graph)[E(sub_graph)$origin == type]) #does not delete shared edges?
    # plot(s_sub_graph, vertex.label=NA, vertex.size=4)
    # components(primed_sub_graph)$membership
    if(del_shared){
      s_sub_graph <- delete.edges(s_sub_graph,E(s_sub_graph)[E(s_sub_graph)$origin == "naive_primed"])
    }
    
    #delete vertices with no connections
    sub_graph_net <- delete.vertices(s_sub_graph,which(degree(s_sub_graph)<1))
    # plot(s_sub_graph_net, vertex.label=NA, vertex.size=4)
    
    #catch empty networks
    #nodes that partake in trans interactions
    nets_list <- list()
    if(vcount(sub_graph_net)==0 & ecount(sub_graph_net)==0){
      nets_list[[1]] <- data.frame(sub_net=1, gene_names=NA, total_vertexes=NA, total_edges=NA,
                                  hindiii_stats=NA, frag_out=NA, distance=NA, avg_distance=NA)
      
    }else{
      #decompose and get information of individual networks (neighbourhoods)
      all_nets <- decompose.graph(sub_graph_net)
      for(g in seq_along(all_nets)){
        #reduce the complexity of large networks by looking only at communities 
        #...
        gene_names <- paste(unique(V(all_nets[[g]])$prot_genes), collapse=",") #protein coding genes only
        
        hindiii <- table(sapply(strsplit(unique(V(all_nets[[g]])$name), "_"), "[[", 1))
        hindiii_stats <- paste(names(hindiii), hindiii, collapse=",")
        
        #retain all fragment ids with connectivity info
        frags <- sapply(unique(V(all_nets[[g]])$name), function(x) length(E(all_nets[[g]])[from(x)]))
        
        #get distance info
        distance <- paste(E(all_nets[[g]])$distance, collapse = ",")
        avg_distance <- mean(E(all_nets[[g]])$distance, na.rm = TRUE)
        
        #are the genes within these fragments baited?
        frag_baited <- names(frags) %in% baited
        
        frag_anno <- paste0(names(frags), "-", frags)
        frag_b <- paste0(frag_anno[frag_baited], "-", rep("b", sum(frag_baited*1)))
        frag_oe <- paste0(frag_anno[!frag_baited], "-", rep("oe", sum(!frag_baited*1)))
        frag <- c(frag_b, frag_oe)
        
        frag_out <- paste(frag, collapse = ",")
        
        total_vertexes <- length(V(all_nets[[g]]))
        total_edges <- length(E(all_nets[[g]]))
        
        nets_list[[g]] <- data.frame(sub_net=g, gene_names, total_vertexes, 
                                     total_edges, hindiii_stats, frag_out, 
                                     distance, avg_distance)
        #other stats....
      }
      
    }
    
    #progress bar
    info <- sprintf("%d%% done", (round(sg/length(all_subgraphs)*100)))
    setTkProgressBar(pb, sg, sprintf("Analyse network (%s)", info), info)
    #net_data 
    all_net_list[[as.character(sg)]] <- plyr::ldply(nets_list, rbind)
    
  }
  #close progress bar
  close(pb)
  return(all_net_list)
}




#' Summarise the results of analyse_net function
#' 
#' @param ana_net_naive Analysed naive network
#' @param ana_net_primed Analysed primed network
#' @param out file path to which to save the summariesed object
#' 
#' @return data.frame
summarise_net_analysis <- function(ana_net_naive, ana_net_primed, out){
  #summarise the results of analyse_net  
  
  if(file.exists(out)){
    print(paste("File already exists and will be read from:", out))
  
    
    np_summary <- readRDS(out)
    return(np_summary)
  }
  
  
  net_primed_df <- plyr::ldply(ana_net_primed, rbind)
  net_naive_df <- plyr::ldply(ana_net_naive, rbind)
  
  
  #Aggregate results
  
  n_sum <- aggregate(x = net_naive_df[c("total_vertexes")], 
                     by = net_naive_df[c(".id")], 
                     FUN = function(subnet){
                       sum(subnet)})
  
  n_sum_e <- aggregate(x = net_naive_df[c("total_edges")], 
                       by = net_naive_df[c(".id")], 
                       FUN = function(subnet){
                         sum(subnet)})
  
  n_label <- aggregate(x = net_naive_df[c("gene_names")], 
                       by = net_naive_df[c(".id")], 
                       FUN = function(subnet){
    gsub("^,*|(?<=,),|,*$","", gsub("NA", "", paste(unlist(strsplit(as.character(subnet), ",")), collapse=",")), perl=TRUE)
  })
  
  p_label <- aggregate(x = net_primed_df[c("gene_names")], 
                       by = net_primed_df[c(".id")], 
                       FUN = function(subnet){
    gsub("^,*|(?<=,),|,*$","", gsub("NA", "", paste(unlist(strsplit(as.character(subnet), ",")), collapse=",")), perl=TRUE)
  })
  
  n_len <- aggregate(x = net_naive_df[c("sub_net")], 
                     by = net_naive_df[c(".id")], 
                     FUN = function(subnet){
                       length(subnet)})
  
  n_distance <- aggregate(x = net_naive_df[c("avg_distance")], 
                          by = net_naive_df[c(".id")], 
                          FUN = function(subnet){
                            mean(subnet)})
  
  p_sum_e <- aggregate(x = net_primed_df[c("total_edges")], 
                       by = net_primed_df[c(".id")], 
                       FUN = function(subnet){
                         sum(subnet)})
  
  p_sum <- aggregate(x = net_primed_df[c("total_vertexes")], 
                     by = net_primed_df[c(".id")], 
                     FUN = function(subnet){
                       sum(subnet)})
  
  p_len <- aggregate(x = net_primed_df[c("sub_net")], 
                     by = net_primed_df[c(".id")], 
                     FUN = function(subnet){
                       length(subnet)})
  
  p_distance <- aggregate(x = net_primed_df[c("avg_distance")], 
                          by = net_primed_df[c(".id")], 
                          FUN = function(subnet){
                            mean(subnet)})
  
  #how many are 2^20 long-range interactions
  p_distance_count <- aggregate(x = net_primed_df[c("distance")], 
                                by = net_primed_df[c(".id")], 
                                FUN = function(subnet){
    dists <- unlist(strsplit(as.character(subnet), ","))
    num_long <- 0
    for(x in seq_along(dists)){
      if(!is.na(dists[x])){
        if(as.double(dists[x]) > 2^20){
          num_long <- num_long + 1
        }
      }
    }
    return(num_long)  
  })
  
  n_distance_count <- aggregate(x = net_naive_df[c("distance")], by = net_naive_df[c(".id")], 
                                FUN = function(subnet){
    dists <- unlist(strsplit(as.character(subnet), ","))
    num_long <- 0
    for(x in seq_along(dists)){
      if(!is.na(dists[x])){
        if(as.double(dists[x]) > 2^20){
          num_long <- num_long + 1
        }
        
      }
      
    }
    return(num_long)  
  })
  
  
  n_summary <- plyr::join_all(list(n_sum, n_sum_e, n_len, n_distance, n_distance_count, n_label), 
                              by = ".id", type="left")
  p_summary <- plyr::join_all(list(p_sum, p_sum_e, p_len, p_distance, p_distance_count, p_label), 
                              by = ".id", type="left")
  
  np_summary <- merge(n_summary, p_summary, by=".id", suffixes = c("_n","_p"), all=TRUE)
  
  names(np_summary) <- c(".id", "n_total_vertexes", "n_total_edges", "n_sub_net", "n_avg_distance",
                         "n_distance_count", "n_label", "p_total_vertexes", "p_total_edges", "p_sub_net",
                         "p_avg_distance", "p_distance_count", "p_label")
  
  np_summary[is.na(np_summary)] <- 0
  
  np_summary$v <- np_summary$n_total_vertexes-np_summary$p_total_vertexes
  np_summary$n <- np_summary$n_sub_net-np_summary$p_sub_net
  np_summary$e <- np_summary$n_total_edges-np_summary$p_total_edges
  
  
  np_summary$dist_orig <- distance_origin(np_summary$v, np_summary$e)
  #need to add the positive negative to identify it as naive or primed
  np_summary$dist_orig_sign <- np_summary$dist_orig
  np_summary$dist_orig_sign[np_summary$v <= 0 & np_summary$e <= 0] <- -np_summary$dist_orig_sign[np_summary$v <= 0 & np_summary$e <= 0] 
  #remove points on the sides
  np_summary$dist_orig_sign[np_summary$v >= 0 & np_summary$e <= 0] <- NA
  np_summary$dist_orig_sign[np_summary$v <= 0 & np_summary$e >= 0] <- NA
  
  
  saveRDS(np_summary, out)
  
  return(np_summary)
  
}


distance_origin <- function(x,y){
  hypotenuse <- sqrt(x^2+y^2)
  return(hypotenuse)
}



#' Add subnetwork ID and distance to origin to network
#' 
#' @param network An igraph network
#' @param np_summary Summary of network produced by summarise_net_analysis function
#' @param all_subgraphs Output of decompose.graph performed on the input igraph network
#' 
#' @return igraph network
add_subnet <- function(network, np_summary, all_subgraphs){
  #Add the subnetwork ID produced in subnet_analysis to the network
  
  #initialise progress bar
  pb <- tkProgressBar("Progress bar", "Some information in %",
                      0, length(all_subgraphs), 0, width=300)
  
  # single_nodes <- 0
  dist_orig_sign_vector <- vector("numeric", length=length(V(network)))
  subnet_vector <- vector("character", length=length(V(network)))
  for(net in seq_along(all_subgraphs)){
    net_row <- which(np_summary$.id == net)
    #if net only contains a node and no edges
    #will correspond to vertices with trans interactions which have been removed
    # if(length(net_row) == 0){ 
    #   # print(vcount(all_subgraphs[[net]]), ecount(all_subgraphs[[net]]))
    #   single_nodes <- single_nodes +1
    # }else{
    select_v <- V(network)$name %in% V(all_subgraphs[[net]])$name
    # V(network)$dist_orig_sign[select_v] <- np_summary$dist_orig_sign[net_row]
    # V(network)$subnet[select_v] <- np_summary$.id[net_row]
    dist_orig_sign_vector[select_v] <- np_summary$dist_orig_sign[net_row]
    subnet_vector[select_v] <- np_summary$.id[net_row]
    # }
    #progress bar
    info <- sprintf("%d%% done", (round(net/length(all_subgraphs)*100)))
    setTkProgressBar(pb, net, sprintf("Adding subnet ID (%s)", info), info)
  }
  
  V(network)$dist_orig_sign <- dist_orig_sign_vector
  V(network)$subnet <- subnet_vector
  
  # print(paste("Single nodes:", single_nodes))
  #close progress bar
  close(pb)
  return(network)
}




#Ploting arcs for individual sub-networks --------------------------------------

get_link_coord <- function(sub_graph, hindiii_lookup, del_type, dist_filt=0){
  #remove primed or naive connections
  t_subgraph <- delete.edges(sub_graph, E(sub_graph)[E(sub_graph)$origin == del_type])
  #delete vertices with no connections
  t_graph <- delete.vertices(t_subgraph,which(igraph::degree(t_subgraph)<1))
  
  net_ends <- data.frame(ends(t_graph, E(t_graph)))
  #need to add index to preseve order
  net_start <- data.frame(ID=net_ends$X1, idx=seq(1:nrow(net_ends)))
  net_end <- data.frame(ID=net_ends$X2, idx=seq(1:nrow(net_ends)))
  net_starts_loc <- merge(net_start, hindiii_lookup, by="ID")
  net_ends_loc <- merge(net_end, hindiii_lookup, by="ID")
  #reorder
  net_starts_loc <- net_starts_loc[c("chrom", "start", "end", "ID", "idx")]
  net_ends_loc <- net_ends_loc[c("chrom", "start", "end", "ID", "idx")]
  
  net_starts_loc <- net_starts_loc[order(net_starts_loc$idx),]
  net_ends_loc <- net_ends_loc[order(net_ends_loc$idx),]
  
  #filter by distance
  net_dists <- abs((net_starts_loc$start+(net_starts_loc$end-net_starts_loc$start)/2)-
                     (net_ends_loc$start+(net_ends_loc$end-net_ends_loc$start)/2))
  net_starts_loc <- net_starts_loc[net_dists > dist_filt,]
  net_ends_loc <- net_ends_loc[net_dists > dist_filt,]
  
  return(list(net_starts_loc,net_ends_loc))
  
}

#' Get bedpe format for subgraph and its interacting hindiii fragments
#' bedpe format start 0-based coordinates end 1-based coordinates
#' 
#' @param sub_graph one of the network from decompose.graph function
#' @param hindiii_lookup hindiii coordinate lookup table
#' @param del_type Analysing naive, then delete primed and vice versa
#' @param mark_dist mark interaction longer than...
#' @param true_bedpe convert start coordinate to 0-based (not really needed for mid coordinates)
#' 
#' @return bedpe 
get_bedpe <- function(sub_graph, hindiii_lookup, del_type, mark_dist=0, true_bedpe=FALSE){
  #remove primed or naive connections
  t_subgraph <- delete.edges(sub_graph, E(sub_graph)[E(sub_graph)$origin == del_type])
  #delete vertices with no connections
  t_graph <- delete.vertices(t_subgraph,which(igraph::degree(t_subgraph)<1))
  
  net_ends <- data.frame(ends(t_graph, E(t_graph)))
  #include score
  net_ends$score <- edge_attr(t_graph, 
                              paste("score", setdiff(c("naive", "primed"), del_type), 
                                    sep="_"))
  # net_ends$score <- E(t_graph)$score_naive
  #need to add index to preseve order
  net_start <- data.frame(ID=net_ends$X1, score=net_ends$score, idx=seq(1:nrow(net_ends)))
  net_end <- data.frame(ID=net_ends$X2, score=net_ends$score, idx=seq(1:nrow(net_ends)))
  net_start_loc <- merge(net_start, hindiii_lookup, by="ID")
  net_end_loc <- merge(net_end, hindiii_lookup, by="ID")
  
  #bedpe
  net_bedpe <- merge(net_start_loc, net_end_loc, by = c("idx","score"), suffixes = c("_1", "_2"))
  #reorder
  net_bedpe <- net_bedpe[c("chrom_1", "start_1", "end_1", 
                           "chrom_2", "start_2", "end_2", 
                           "ID_1", "ID_2", "idx", "score")]
  
  if(true_bedpe){
    net_bedpe$start_1 <- net_bedpe$start_1-1 
    net_bedpe$start_2 <- net_bedpe$start_2-1 
  }
  
  #substitute any NA scores with 0
  net_bedpe$score[is.na(net_bedpe$score)] <- 0
  
  #add distance
  net_bedpe$distance <- abs((net_bedpe$start_1+(net_bedpe$end_1-net_bedpe$start_1)/2)-
                              (net_bedpe$start_2+(net_bedpe$end_2-net_bedpe$start_2)/2))
  net_bedpe$marked_dist <- 0
  net_bedpe$marked_dist[net_bedpe$distance > mark_dist] <- 1
  return(net_bedpe)
}





#'Plot subnetwork using the sushi package
#'
#'@param primed_bedpe BEDPE for primed interactions from get_bedpe function
#'@param naive_bedpe BEDPE for naive interactions from get_bedpe function
#'@param out_path where to save plot pdf
#'@param hg38_chrom_sizes chromosome sizes used for plot range
#'@param zoom either uses whole chromosome or +/- 100th of the bedpe range as a zoom in
#'@param keep_scales_equal ymax for naive and primed equal
#'
#'@return sushi plot
plot_subnetwork_sushi <- function(primed_bedpe, naive_bedpe, hg38_chrom_sizes, 
                                  out_path=NA, zoom=FALSE, keep_scales_equal=TRUE){
  
  #limit to coordinates in bedpe files
  if(keep_scales_equal){
    
    p_max <- max(primed_bedpe$score)
    n_max <- max(naive_bedpe$score)
    
    if(p_max > n_max){
     naive_sf <- (p_max/n_max) *1.04
     primed_sf <- 1.04
    }else{
      naive_sf <- 1.04
      primed_sf <- (n_max/p_max)*1.04
    }
    
  }else{
    naive_sf <- 1.04
    primed_sf <- 1.04
  }

  if(zoom){
    start_coord <- min(primed_bedpe$start_1, primed_bedpe$start_2, 
                       naive_bedpe$start_1, naive_bedpe$start_2)
    end_coord <- max(primed_bedpe$end_1, primed_bedpe$end_2, 
                     naive_bedpe$end_1, naive_bedpe$end_2)
    
    chrom            <- unique(primed_bedpe$chrom_1)
    chromstart       <- start_coord - as.integer((end_coord-start_coord)/100)
    chromend         <- end_coord + as.integer((end_coord-start_coord)/100)
  } else{
    chrom            <- unique(primed_bedpe$chrom_1)
    chromstart       <- 1
    chromend         <- hg38_chrom[chrom]$size  #181538259
  } 
  
  # pdf("/media/chovanec/My_Passport/CHiC_naive_primed/network/network_analysis/2301_changing_sushi_2^20.pdf", 
  #width=10, height=5)
  
  # cat(chrom, chromstart, chromend)
  if(!is.na(out_path)){
    pdf(out_path, width=10, height=5)  
  }
  
  
  layout(matrix(c(1,1,2,2), nrow = 4, ncol = 1, byrow = TRUE))
  pbpe <- plotBedpe(primed_bedpe,chrom,chromstart,chromend,
                    heights = primed_bedpe$score,plottype="loops",
                    colorby=primed_bedpe$marked_dist,
                    colorbycol=SushiColors(3), ymax = primed_sf)
  labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
  legend("topright",inset =0.01,legend=c(paste("Primed interaction distance", 
                                               chrom, chromstart, chromend, sep=" ")),
         col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
  axis(side=2,las=2,tcl=.2)
  mtext("CHiCAGO score",side=2,line=1.75,cex=.75,font=2)
  
  pbpe <- plotBedpe(naive_bedpe,chrom,chromstart,chromend,
                    heights = naive_bedpe$score,plottype="loops",
                    colorby=naive_bedpe$marked_dist,
                    colorbycol=SushiColors(3), ymax=naive_sf)
  labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
  legend("topright",inset =0.01,legend=c(paste("Naive interaction distance", 
                                               chrom, chromstart, chromend, sep=" ")),
         col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
  axis(side=2,las=2,tcl=.2)
  mtext("CHiCAGO score",side=2,line=1.75,cex=.75,font=2)
  
  if(!is.na(out_path)){
  dev.off()
  }
  
}







#' Plot subnetwork using the karyotype package
#' 
#' 
plot_subnetwork <- function(subnet, hindiii_lookup, long_range=2^20, long_only=FALSE, 
                            use_zoom=FALSE, protein_genes=NA, rna_genes=NA, title=""){
  
  n_links_subnet <- get_link_coord(subnet, hindiii_lookup, "primed")
  p_links_subnet <- get_link_coord(subnet, hindiii_lookup, "naive")
  
  n_links_subnet_long <- get_link_coord(subnet, hindiii_lookup, "primed", long_range)
  p_links_subnet_long <- get_link_coord(subnet, hindiii_lookup, "naive", long_range)
  
  n_starts <- makeGRangesFromDataFrame(n_links_subnet[[1]])
  n_ends <- makeGRangesFromDataFrame(n_links_subnet[[2]])
  
  p_starts <- makeGRangesFromDataFrame(p_links_subnet[[1]])
  p_ends <- makeGRangesFromDataFrame(p_links_subnet[[2]])
  
  if(nrow(n_links_subnet_long[[1]]) > 0){
    n_starts_long <- makeGRangesFromDataFrame(n_links_subnet_long[[1]])
    n_ends_long <- makeGRangesFromDataFrame(n_links_subnet_long[[2]])
    no_n_long <- FALSE
  }else{
    no_n_long <- TRUE
  }
  if(nrow(p_links_subnet_long[[1]]) > 0){
    p_starts_long <- makeGRangesFromDataFrame(p_links_subnet_long[[1]])
    p_ends_long <- makeGRangesFromDataFrame(p_links_subnet_long[[2]])
    no_p_long <- FALSE
  }else{
    no_p_long <- TRUE
  }
  
  chrom <- levels(seqnames(n_starts))
  
  stopifnot(chrom == levels(seqnames(p_starts))) #test both n and p same chr
  
  #zoom in on region +/-5kb otherwise whole chromosome ploted
  zoom_start = min(c(start(n_starts), start(p_starts), start(n_ends), start(p_ends)))-5000
  zoom_end = max(c(end(n_starts), end(p_starts), end(n_ends), end(p_ends)))+5000
  zoom_gr <- GRanges(seqnames = chrom, strand = "*", 
                     ranges = IRanges(start=zoom_start, width=zoom_end-zoom_start))
  
  pp <- karyoploteR::getDefaultPlotParams(plot.type = 2)
  pp$data1height <- 50
  pp$data2height <- 100
  pp$ideogramheight <- 10
  if(use_zoom){
    kp <- karyoploteR::plotKaryotype(chromosomes=c(chrom), plot.params = pp, 
                                     plot.type = 2, cex=0.6, zoom = zoom_gr , 
                                     main=title) 
  }else{
    kp <- karyoploteR::plotKaryotype(chromosomes=c(chrom), plot.params = pp, 
                                     plot.type = 2, cex=0.6, main=title) 
  }
  
  if(!is.na(protein_genes[1])){
    karyoploteR::kpPlotMarkers(kp, data=protein_genes, labels=protein_genes$gene_name, 
                               data.panel = 2, r0=0, r1=0.1, cex=0.5, 
                               adjust.label.position = FALSE)
    r0_pos = 0.5
    r1_pos = 1
  }
  if(!is.na(rna_genes[1])){
    karyoploteR::kpPlotMarkers(kp, data=rna_genes, labels=rna_genes$gene_name, 
                               data.panel = 2, r0=0, r1=0.1, cex=0.5, 
                               adjust.label.position = FALSE,
                  label.color = "grey")
    r0_pos = 0.5
    r1_pos = 1
  }
  
  if(is.na(protein_genes[1]) & is.na(rna_genes[1])){
    r0_pos = 0
    r1_pos = 0.5
  }
  
  
  #links
  if(long_only & !no_p_long){
    karyoploteR::kpPlotLinks(kp, data=p_starts_long, data2=p_ends_long, col="#FF2B2B", 
                             data.panel=2, r0=r0_pos, r1=r1_pos, lwd=0.3)
  }else if(long_only & !no_n_long){
    karyoploteR::kpPlotLinks(kp, data=n_starts_long, data2=n_ends_long, col="#3B85F5", 
                             data.panel=1, r0=0, r1=1, lwd=0.3)
  }else{
    karyoploteR::kpPlotLinks(kp, data=n_starts, data2=n_ends, col="#DBDBDB", 
                             data.panel=1, r0=0, r1=1, lwd=0.1)
    karyoploteR::kpPlotLinks(kp, data=p_starts, data2=p_ends, col="#DBDBDB", 
                             data.panel=2, r0=r0_pos, r1=r1_pos, lwd=0.1)
    if(!no_p_long){
      karyoploteR::kpPlotLinks(kp, data=p_starts_long, data2=p_ends_long, col="#FF2B2B", 
                               data.panel=2, r0=r0_pos, r1=r1_pos, lwd=0.3)
    }
    if(!no_n_long){
      karyoploteR::kpPlotLinks(kp, data=n_starts_long, data2=n_ends_long, col="#3B85F5", 
                               data.panel=1, r0=0, r1=1, lwd=0.3)
    }
  }
}







#Community analysis ------------------------------------------------------------


nodes2granges <- function(nodes, lookup_dt){
  if(grepl(",", nodes)){
    nodes <- unlist(strsplit(as.character(nodes), ","))
  }else if(grepl(";", nodes)){
    nodes <- unlist(strsplit(as.character(nodes), ";"))
  }
  
  grange_dt <- lookup_dt[nodes, nomatch=0L]
  gr <- GenomicRanges::makeGRangesFromDataFrame(grange_dt, keep.extra.columns=TRUE)
  
  return(gr)
}

