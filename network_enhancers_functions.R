




#'create a lookup table for enhancers, hindiii, OSN together
#'all inputs are GenomicRanges
#'
#'@param enhancer_gr enhancer GRanges object (contains mcols score and type i.e. E or SE)
#'@param hindiii_gr hindiii GRanges
#'@param OSN_gr OSN GRanges
#'
#'
ROSE_OSN_HindIII_table <- function(enhancer_gr, hindiii_gr, OSN_gr){

  #enhancer hindiii overlap
  enh_hindiii <- findOverlaps(enhancer_gr, hindiii_gr)
  
  #enhancer name + hindiii ID table
  enh_hindiii_df <- data.frame(name=enhancer_gr[queryHits(enh_hindiii)]$name, 
                                          ID=hindiii_gr[subjectHits(enh_hindiii)]$ID)
  
  #which OSN actually overlaps enhancer, don't care about other OSN's in hindiii if not overlapping enhancer
  
  #OSN overlap with hindiii fragments and then merge with enhancers base on hindiii ID
  #enhancer OSN overlap
  
  enh_OSN_ovl <- findOverlaps(enhancer_gr, OSN_gr)
  #subset OSN to only those overlapping enhancer
  
  OSN_enhancer_df <- data.frame(name=enhancer_gr[queryHits(enh_OSN_ovl)]$name, 
                                OSN=OSN_gr[subjectHits(enh_OSN_ovl)]$name)
  
  OSN_enhancer_lookup <- aggregate(x = OSN_enhancer_df["OSN"], 
                                  by= OSN_enhancer_df["name"], 
                                  FUN = function(ids){
                                    paste(ids, collapse = ",")
                                  })
  
  
  enh_hindiii_OSN <- merge(enh_hindiii_df, OSN_enhancer_lookup, 
                           by = "name", all = TRUE)
  
  
  
  #aggregate by enhancer (multiple hindiii fragments within enhancer)
  enh_hindiii_OSN_lookup <- aggregate(x = enh_hindiii_OSN[c("ID", "OSN")], 
                                  by=enh_hindiii_OSN["name"], 
                                  FUN = function(ids){
                                    paste(ids, collapse = ";")
                                  })

  #only want these columns
  enh_hindiii_OSN_all <- merge(data.frame(enhancer_gr)[c("name", "seqnames", "start", "end", 
                                                         "width", "strand", "score","type")], 
                               enh_hindiii_OSN_lookup, 
                               by = "name", all.x = TRUE)
  
  return(enh_hindiii_OSN_all)
}



#' Get mid of two coordinates start end
#' @param start
#' @param end
#' 
#' @example 
#' mid_coord(c(1,2,3), c(3,4,5))
#' 
#' returns 
#' c(2,3,4)
#' 
#' @return vector
mid_coord <- function(start, end){
  width <- end-start
  start+width/2
}




#'Deduplicated concatenated (aggregated) items by a defined column
#'
#'@param dedup_df data.frame with aggregated columns
#'@param dedup_with column to deduplicate with
#'@param dedup_column column(s) to deduplicated based on dedup_with column
#'
#'@return deduplicated data.frame
#'
#'If you want to count how many interactions are established with the same enhancer
#'do not deduplicate
#'
deduplicated_collapsed <- function(dedup_df, dedup_with="name_naive", dedup_column="type_naive"){
  
  unique_items <- sapply(strsplit(dedup_df[[dedup_with]], ";"), function(x) !duplicated(x))
  # return(unique_items)
  
  for(col in seq_along(dedup_column)){
    
    col_split <- strsplit(dedup_df[[dedup_column[col]]], ";")
    
    dedup_df[dedup_column[col]] <- sapply(seq(1,nrow(dedup_df), 1), function(x)
    {paste(col_split[[x]][unique_items[[x]]], collapse = ";")})
    
  }
  return(dedup_df)
}







#' Count number of OSN sites gained, lost
#' 
#' @param count_all when true, split multiple OSN in one hindiii and count all of them
#' 
#' "n"  1
#' "np"  2
#' "p"  3
count_osn_gain_loss <- function(main_osn_vector, other_osn_vector, 
                                type=c("gained", "lost"), return_vector=FALSE){
  #naive OSN can only have 1 or 2
  #primed OSN can only have 2 or 3
  
  if(!type %in% c("gained", "lost")){
    stop("Only gained or lost type accepted")
  }
  if(length(type) > 1){
    stop("Select only one type")
  }
  
  
  if(length(main_osn_vector) == 0){
    out <- list(gained=0, 
                lost=0,
                retained=0,
                none=0)
    return(out)
  }
  
  if(type == "gained"){
    poss_vals <- c("3", "2")
    poss_vals_other <- 1
  }
  if(type == "lost"){
    poss_vals <- c("1", "2")
    poss_vals_other <- 3
  }
  
  #function to resolve OSN  
    #for enhancers with multiple OSN sites:
      #if gained enhancer only count gained osn
      #if lost enhancer only count lost osn
      #if retained enhancer only count retained
  resolve_osn <- function(val, poss_vals){
    splt_val <- unique(unlist(strsplit(val, ","))) #return only a single value of multiple cell specific OSN
    if(length(splt_val) == 1){
      return(splt_val)
    }
    
    if(poss_vals[1] %in% splt_val & poss_vals[2] %in% splt_val){
      return(poss_vals[1]) #the cell specific type
    }else if(poss_vals[1] %in% splt_val){
      return(poss_vals[1])
    }else if(poss_vals[2] %in% splt_val){
      return(poss_vals[2])
    }
  }
 
  
  osn_counts <- table(sapply(main_osn_vector, 
                             function(x) ifelse(grepl(",", x), resolve_osn(x, poss_vals), x)), 
                      sapply(other_osn_vector, 
                             function(x) ifelse(grepl(",", x), resolve_osn(x, poss_vals_other), x)),
                      useNA="always")
  
  #return group categories instead of counts
  if(return_vector){
    osn_counts_main <- sapply(main_osn_vector, 
                               function(x) ifelse(grepl(",", x), resolve_osn(x, poss_vals), x))
    osn_counts_other <- sapply(other_osn_vector, 
                               function(x) ifelse(grepl(",", x), resolve_osn(x, poss_vals_other), x))
    
    stopifnot(length(osn_counts_main) == length(osn_counts_other))
    
    group_ids <- merge_vect(as.character(osn_counts_main), as.character(osn_counts_other))
    
    group_ids <- gsub("1", paste("enh", type, "osn_lost", sep = "_"), group_ids)
    group_ids <- gsub("2", paste("enh", type, "osn_retained", sep = "_"), group_ids)
    group_ids <- gsub("3", paste("enh", type, "osn_gained", sep = "_"), group_ids)
    
  
    return(group_ids)
    
  }
  
  if(type == "lost"){
    
    primed_g <- unname(osn_counts[is.na(rownames(osn_counts)), "3"])
    naive_l <- unname(osn_counts["1", is.na(colnames(osn_counts))])
    shared_r <- unname(osn_counts["2", is.na(colnames(osn_counts))])
    none <- unname(osn_counts[is.na(rownames(osn_counts)), is.na(colnames(osn_counts))])
    
  }else if(type == "gained"){
    
    primed_g <- unname(osn_counts["3", is.na(colnames(osn_counts))])
    naive_l <- unname(osn_counts[is.na(rownames(osn_counts)), "1"])
    shared_r <- unname(osn_counts["2", is.na(colnames(osn_counts))])
    none <- unname(osn_counts[is.na(rownames(osn_counts)), is.na(colnames(osn_counts))])
  }
  
  
  out <- list(gained=ifelse(is.na(primed_g), 0, primed_g), 
              lost=ifelse(is.na(naive_l), 0, naive_l),
              retained=ifelse(is.na(shared_r), 0, shared_r),
              none=ifelse(is.na(none), 0, none))
  
  return(out)
}
    



 
    
#'Count retained enhancers and their OSN 
#'
#'@param retained data.frame with enhancers in both naive and primed
#'@param Naive_enh GRanges object with naive enhancers
#'@param Primed_enh GRanges object with primed enhancers
#'@param OSN_enhancer_lookup data.table key short_name with overlaps of all OSN and all enhancers use to find lost_gained and gained_lost OSN's
#'
#'@note 
#' test <- apply(retained, 1, function(x){
#'   ovlp <- Naive_enh[Naive_enh$short_name %in% x["name_naive"]] %over% Primed_enh[Primed_enh$short_name %in% x["name_primed"]]
#'   return(ovlp)})
#'
count_retained <- function(retained, Naive_enh, Primed_enh, OSN_enhancer_lookup){
  
  # x <- c(ID="chr10_28855",chromhmm_naive="Active", name_naive="2_np_4353",
  #        score_naive="454", type_naive="SE", OSN_naive="1",
  #        chromhmm_primed="Active", name_primed="1_pp_2499", score_primed="2826" ,
  #        type_primed="E", OSN_primed="NA", b_ID="chr10_28932")
  
  retained_categories <- apply(retained, 1, function(x){
    #do the enhancers overlap?
    
    ovlp <- Naive_enh[Naive_enh$short_name %in% x["name_naive"]] %over% Primed_enh[Primed_enh$short_name %in% x["name_primed"]]
    
    #if yes then count only retained OSN
    if(ovlp){
      # ovlp_count %+=% 1
      if(grepl(2, x["OSN_naive"]) & grepl(2, x["OSN_primed"])){
        return("enh_retained_osn_retained")
        #if 1 or 2 and other is NA
      }else if((grepl(2, x["OSN_naive"]) | grepl(1, x["OSN_naive"])) & is.na(x["OSN_primed"])){
        return("enh_retained_osn_lost")
      }else if(is.na(x["OSN_naive"]) & (grepl(2, x["OSN_primed"]) | grepl(3, x["OSN_primed"]))){
        return("enh_retained_osn_gained")
      }else if(is.na(x["OSN_naive"]) & is.na(x["OSN_primed"])){
        return("enh_retained_osn_none")
      }else if(grepl(1, x["OSN_naive"]) & grepl(3, x["OSN_primed"])){
        return("enh_retained_osn_mixed")
      }else if(grepl(1, x["OSN_naive"]) & grepl(2, x["OSN_primed"])){
        return("enh_retained_osn_mixed")
      }else if(grepl(2, x["OSN_naive"]) & grepl(3, x["OSN_primed"])){
        return("enh_retained_osn_mixed")
      }
      print("ovlp")
      print(x)
      #not overlaping, treat as seperate enhancers
    }else{
      # no_ovlp_count %+=% 1
      
      #is not empty then get all OSN (naive and primed)
      n_osn <- OSN_enhancer_lookup[unname(x["name_naive"])]$OSN
      if(is.na(n_osn)){
        n_out <- "enh_lost_osn_none"
      }else if(grepl(2, x["OSN_naive"])){
        n_out <- "enh_lost_osn_retained"
      }else if(grepl(1, x["OSN_naive"])){
        n_out <- "enh_lost_osn_lost"
      }else if(grepl("3", n_osn)){
        n_out <- "enh_lost_osn_gained"
      }else{
        cat("naive", n_osn)
      }
      
      #for primed OSN
      p_osn <- OSN_enhancer_lookup[unname(x["name_primed"])]$OSN
      if(is.na(p_osn)){
        p_out <- "enh_gained_osn_none"
      }else if(grepl(2, x["OSN_primed"])){
        p_out <- "enh_gained_osn_retained"
      }else if(grepl("3", x["OSN_primed"])){
        p_out <- "enh_gained_osn_gained"
      }else if(grepl(1, p_osn)){
        p_out <- "enh_gained_osn_lost"
      }else{
        cat("primed", p_osn)
      }
      
      return(c(n_out, p_out))
      
      print("not_ovlp")
      print(x)
    }
  })
  
  return(retained_categories)
}




#'Calculate log2 odds ratio
#'function used by Paula in elife paper
#'
#'@param ctg count matrix
#'
#'Haldane-Anscombe correction - Just adding 0.5 to each of the cells and then 
#'cam ise loddsratio(ctg, correct = TRUE) for 0.5 correction
#'
logOddsRatio <- function(ctg){
  res3 = matrix(nrow=nrow(ctg), ncol=ncol(ctg))
  rownames(res3) <- rownames(ctg)
  colnames(res3) <- colnames(ctg)
  print(fisher.test(ctg, simulate.p.value = T, B=1e5))
  ## Now compute odds ratios
  for (i in 1:nrow(ctg)){
    for(j in 1:ncol(ctg)){
      res3[i,j] = ctg[i,j]/sum(ctg[i,-j])/(sum(ctg[-i,j])/sum(ctg[-i,-j]))
    }
  }
  
  library(gplots)
  heatmap.2(log(res3, 2), dendrogram = "none", Rowv = NA, Colv = NA, col=bluered(100), margins=c(15,15), breaks=seq(from = -1, to=1,by = 0.02), trace="none")
  return(res3)
}












#plot_range
#gene_to_plot_df
#chromhmm_to_plot_df
#b_bed_center_plot
#oe_bed_center_plot

#Circos plot data S4 class

setClass(Class = "CircosData",
         slots = c(plot_range = "data.frame",
                   genes = "data.frame",
                   chromhmm = "data.frame",
                   b_interac_center = "data.frame",
                   oe_interac_center = "data.frame",
                   b_interac = "data.frame",
                   oe_interac = "data.frame",
                   origin = "vector",
                   enhancers = "vector",
                   enhancers_bed = "data.frame",
                   osn_bed = "data.frame",
                   chromhmm_raw = "data.frame"),
         prototype = list(plot_range = data.frame(NULL),
                          genes = data.frame(NULL),
                          chromhmm = data.frame(NULL),
                          b_interac_center = data.frame(NULL),
                          oe_interac_center = data.frame(NULL),
                          b_interac = data.frame(NULL),
                          oe_interac = data.frame(NULL),
                          origin = vector("character"),
                          enhancers = vector("character"),
                          enhancers_bed = data.frame(NULL),
                          osn_bed = data.frame(NULL),
                          chromhmm_raw = data.frame(NULL)
                          )
         )

setGeneric("plotrange", function(object) standardGeneric("plotrange"))
setMethod("plotrange", "CircosData", function(object) object@plot_range)

#generic set in network_enhancer_class
# setGeneric("genes", function(object) standardGeneric("genes"))
setMethod("genes", "CircosData", function(object) object@genes)

setGeneric("chromhmm", function(object) standardGeneric("chromhmm"))
setMethod("chromhmm", "CircosData", function(object) object@chromhmm)

setGeneric("bc", function(object) standardGeneric("bc"))
setMethod("bc", "CircosData", function(object) object@b_interac_center)

setGeneric("oec", function(object) standardGeneric("oec"))
setMethod("oec", "CircosData", function(object) object@oe_interac_center)

setGeneric("b", function(object) standardGeneric("b"))
setMethod("b", "CircosData", function(object) object@b_interac)

setGeneric("oe", function(object) standardGeneric("oe"))
setMethod("oe", "CircosData", function(object) object@oe_interac)

setGeneric("origin", function(object) standardGeneric("origin"))
setMethod("origin", "CircosData", function(object) object@origin)

#generic set in network_enhancer_class
# setGeneric("enhancers", function(object) standardGeneric("enh"))
setMethod("enhancers", "CircosData", function(object, type) object@enhancers)

setGeneric("enhBED", function(object) standardGeneric("enhBED"))
setMethod("enhBED", "CircosData", function(object) object@enhancers_bed)

setGeneric("osnBED", function(object) standardGeneric("osnBED"))
setMethod("osnBED", "CircosData", function(object) object@osn_bed)

setGeneric("chromhmmraw", function(object) standardGeneric("chromhmmraw"))
setMethod("chromhmmraw", "CircosData", function(object) object@chromhmm_raw)



#create an enhancer/OSN BED track(s)
#Use ROSE BED and OSN BED



#' Plot a circos plot with gene and chromhmm tracks with interactions from a viewpoint
#' 
#' @param interactions_object Interactions class object from get_interaction_counts with "raw" setting
#' @param hindiii_coord
#' @param chromhmm_gr GenomicRanges object 
#' 
prepare_circlize_data <- function(interactions_object, hindiii_coord, chromhmm_gr,
                                  enhancer_gr, osn_gr,
                                  type=c("naive", "primed"),
                                  plot_type=c("enhancer", "significant", "full"),
                                  other_cell=TRUE,
                                  raw_chromhmm=NULL){
  
  if(!(type %in% c("naive", "primed"))){
    stop("Only 'naive' or 'primed' allowed in type")
  }
  
  if(!(plot_type %in% c("enhancer", "significant", "full"))){
    stop("Only 'enhancer', 'significant' or 'full' allowed in plot_type")
  }

  int_table <- interactions(interactions_object, type)
  # int_table <- interactions(gene_enhancer_raw[[100]], "naive")
  
  #when full only colour by significant possible
  
  if(nrow(int_table) == 0){
    cat("No interaction in interaction table")
    b_bed_center_plot <- data.table(NULL)
    oe_bed_center_plot <- data.table(NULL)
    b_bed_plot_df <- data.table(NULL)
    oe_bed_plot_df <- data.table(NULL)
    int_table[, enhancers:= "B"]
  }else{
    
    #include data for enhancer plots
    #mark which interactions contact enhancers
    enh_table <- enhancers(interactions_object, type)
    
    #in order of enh_table
    enh_rows <- c(which(int_table$from_name %in% enh_table$ID),
                  which(int_table$to_name %in% enh_table$ID))
    
    #add enhancer column
    int_table[, enhancers:= "B"]
    #categorise enhancers
    int_table$enhancers[enh_rows] <- sapply(enh_table[[paste("type", type, sep = "_")]], 
                                            function(x) ifelse("SE" %in% unlist(strsplit(x, ";")), 
                                                               "SE", "E"))
    
    
    
    #get coordinates from hindiii fragment IDs
    b_bed <- hindiii_coord[int_table$from_name]
    oe_bed <- hindiii_coord[int_table$to_name]
    
    b_bed_plot <- cbind(b_bed[,1:3], score=int_table[[paste("score", type, sep = "_")]], 
                        b_bed[,4])
    oe_bed_plot <- cbind(oe_bed[,1:3], score=int_table[[paste("score", type, sep = "_")]], 
                         oe_bed[,4])
    
    #which chromosome is to be plotted, make sure no trans present
    stopifnot(unique(as.character(b_bed_plot$seqnames))==unique(as.character(oe_bed_plot$seqnames)))
    
    chrom <- as.character(unique(b_bed_plot$seqnames))
    
    b_bed_plot_df <- cbind(name=b_bed_plot$seqnames, start=b_bed_plot$start,
                           end=b_bed_plot$end, b_bed_plot[,4:ncol(b_bed_plot)])
    oe_bed_plot_df <- cbind(name=oe_bed_plot$seqnames, start=oe_bed_plot$start,
                            end=oe_bed_plot$end, oe_bed_plot[,4:ncol(oe_bed_plot)])
    
    b_bed_center_plot <- cbind(name=b_bed_plot$seqnames, 
                               start=mid_coord(b_bed_plot$start,b_bed_plot$end),
                               end=mid_coord(b_bed_plot$start,b_bed_plot$end), 
                               b_bed_plot[,4:ncol(b_bed_plot)])
    oe_bed_center_plot <- cbind(name=oe_bed_plot$seqnames, 
                                start=mid_coord(oe_bed_plot$start,oe_bed_plot$end),
                                end=mid_coord(oe_bed_plot$start,oe_bed_plot$end), 
                                oe_bed_plot[,4:ncol(oe_bed_plot)])
  }
  
  
  #if other cell type will be plotted, make sure the plot range is the same
  if(other_cell){
    other_int_table <- interactions(interactions_object, setdiff(c("naive", "primed"), type))
    
    if(nrow(int_table) == 0){
      #if other cell type has no interactions and current has some, use only those
      if(nrow(other_int_table) == 0){
        stop("No interactions in either cell type, will not plot")
      }else{
        
        b_bed_other <- hindiii_coord[other_int_table$from_name]
        oe_bed_other <- hindiii_coord[other_int_table$to_name]
        
        chrom <- as.character(unique(b_bed_other$seqnames))
        
        stopifnot(length(chrom) == 1)
        
        start_other <- min(b_bed_other$start-1, oe_bed_other$start-1)
        end_other <- max(b_bed_other$end, oe_bed_other$end)
        #convert to 0-based bed data frame
        plot_range <- data.frame(name=chrom, 
                                 start=start_other,
                                 end=end_other, 
                                 stringsAsFactors = FALSE)
        
      }
      #else compare which coodrinates are larger and use those for plotting
    }else{
      b_bed_other <- hindiii_coord[other_int_table$from_name]
      oe_bed_other <- hindiii_coord[other_int_table$to_name]
      
      start_current <- min(b_bed_plot$start-1, oe_bed_plot$start-1)
      end_current <- max(b_bed_plot$end, oe_bed_plot$end)

      if(nrow(other_int_table) == 0){
        start <- start_current
        end <- end_current
      }else{
        #get coordinates from hindiii fragment IDs
        
        start_other <- min(b_bed_other$start-1, oe_bed_other$start-1)
        end_other <- max(b_bed_other$end, oe_bed_other$end)
        
        start <- ifelse(start_current <= start_other, start_current, start_other)
        end <- ifelse(end_current >= end_other, end_current, end_other)
        
      }
      #convert to 0-based bed data frame
      plot_range <- data.frame(name=chrom, 
                               start=start,
                               end=end, 
                               stringsAsFactors = FALSE)
    }
    
  }else{
    #convert to 0-based bed data frame
    plot_range <- data.frame(name=chrom, 
                             start=min(b_bed_plot$start-1, oe_bed_plot$start-1),
                             end=max(b_bed_plot$end, oe_bed_plot$end), 
                             stringsAsFactors = FALSE)
  }
  
  
  plotting_gr <- GRanges(seqnames = plot_range$name, 
                         ranges = IRanges(start = plot_range$start, 
                                          end = plot_range$end))
  
  #only keep fully overlapping hindiii intervals
  ovlp <- findOverlaps(chromhmm_gr, plotting_gr, type = "within")
  to_plot_gr <- chromhmm_gr[queryHits(ovlp)]
  #convert to 0-based bed data frame
  chromhmm_to_plot_df <- data.frame(name=seqnames(to_plot_gr), start=start(to_plot_gr)-1, 
                                    end=end(to_plot_gr), 
                                    chromhmm=mcols(to_plot_gr)[paste("chromhmm", type, sep="_")])
  
  #Get genes to plot
  gene_ovlp <- findOverlaps(gtf_gene_gr, plotting_gr)
  gene_to_plot_gr <- gtf_gene_gr[queryHits(gene_ovlp)]
  #convert to 0-based bed data frame
  gene_to_plot_df <- data.frame(name=seqnames(gene_to_plot_gr), start=start(gene_to_plot_gr)-1, 
                                end=end(gene_to_plot_gr), gene_biotype=gene_to_plot_gr$gene_biotype, 
                                gene_name=gene_to_plot_gr$gene_name, strand=strand(gene_to_plot_gr))
  #set direction of arrow
  gene_to_plot_df$strand <- ifelse(gene_to_plot_df$strand == "-", "start", "end")
  
  #if gene goes beyond plotting area, cut the end off
  gene_to_plot_df$start[gene_to_plot_df$start < plot_range$start] <- plot_range$start
  gene_to_plot_df$end[gene_to_plot_df$end > plot_range$end] <- plot_range$end
  
  
  #get enhancer and OSN peaks within plotting range
  enh_ovls <- findOverlaps(enhancer_gr, plotting_gr)
  enhancers_to_plot_gr <- enhancer_gr[queryHits(enh_ovls)]
  #convert to 0-based bed data frame
  enhancers_to_plot_df <- data.frame(name=seqnames(enhancers_to_plot_gr),
                                     start=start(enhancers_to_plot_gr)-1,
                                     end=end(enhancers_to_plot_gr),
                                     type=enhancers_to_plot_gr$type,
                                     ID=enhancers_to_plot_gr$ID,
                                     score=enhancers_to_plot_gr$score)
  
  #if enhancer goes beyond plotting area, cut the end off
  enhancers_to_plot_df$start[enhancers_to_plot_df$start < plot_range$start] <- plot_range$start
  enhancers_to_plot_df$end[enhancers_to_plot_df$end > plot_range$end] <- plot_range$end
  
  
  osn_ovls <- findOverlaps(osn_gr, plotting_gr)
  osn_to_plot_gr <- osn_gr[queryHits(osn_ovls)]
  #convert to 0-based bed data frame
  osn_to_plot_df <- data.frame(name=seqnames(osn_to_plot_gr),
                                     start=start(osn_to_plot_gr)-1,
                                     end=end(osn_to_plot_gr))
  
  #if osn goes beyond plotting area, cut the end off
  osn_to_plot_df$start[osn_to_plot_df$start < plot_range$start] <- plot_range$start
  osn_to_plot_df$end[osn_to_plot_df$end > plot_range$end] <- plot_range$end
  
  
  if(!is.null(raw_chromhmm)){
    #raw chromhmm (not aggregated into hindiii fragments)
    chromhmm_raw_ovls <- findOverlaps(raw_chromhmm, plotting_gr)
    chromhmm_raw_to_plot_gr <- raw_chromhmm[queryHits(chromhmm_raw_ovls)]
    
    #convert to 0-based bed data frame
    chromhmm_raw_to_plot_df <- data.frame(name=seqnames(chromhmm_raw_to_plot_gr),
                                          start=start(chromhmm_raw_to_plot_gr)-1,
                                          end=end(chromhmm_raw_to_plot_gr),
                                          state=mcols(chromhmm_raw_to_plot_gr)[["name"]])
    #if osn goes beyond plotting area, cut the end off
    chromhmm_raw_to_plot_df$start[chromhmm_raw_to_plot_df$start < plot_range$start] <- plot_range$start
    chromhmm_raw_to_plot_df$end[chromhmm_raw_to_plot_df$end > plot_range$end] <- plot_range$end
        
  }else{
    chromhmm_raw_to_plot_df <- data.frame(NULL)
  }
  
  
  
  #out CircosData object
  out <- new("CircosData", plot_range = plot_range,
             genes = gene_to_plot_df,
             chromhmm = chromhmm_to_plot_df,
             b_interac_center = b_bed_center_plot,
             oe_interac_center = oe_bed_center_plot,
             b_interac = b_bed_plot_df,
             oe_interac = oe_bed_plot_df, 
             origin = int_table$new_origin,
             enhancers = int_table$enhancers,
             enhancers_bed = enhancers_to_plot_df,
             osn_bed = osn_to_plot_df,
             chromhmm_raw = chromhmm_raw_to_plot_df)
 
  return(out)
}




#' Plot a circos plot with gene and chromhmm tracks with interactions from a viewpoint
#' 
#' @param interactions_object Interactions class object from get_interaction_counts with "raw" setting
#' @param hindiii_coord
#' @param chromhmm_gr GenomicRanges object 
#' 
plot_circlize <- function(CircosData_object,
                          type=c("naive", "primed"), 
                          plot_mid=TRUE, 
                          enhancer_osn=FALSE,
                          plot_type=c("enhancer", "significant", "full"),
                          plot_genes=c("inline", "stacked"),
                          plot_raw_chromhmm=FALSE){

  if(length(type)>1){
    cat("Using type:", type[1])
    type <- type[1]
  }
  
  if(!(type %in% c("naive", "primed"))){
    stop("Only 'naive' or 'primed' allowed in type")
  }
  
  if(!(plot_type %in% c("enhancer", "significant", "full"))){
    stop("Only 'enhancer', 'significant' or 'full' allowed in plot_type")
  }
  
  if(!(plot_genes %in% c("inline", "stacked"))){
    stop("Only 'inline' or 'stacked' allowed in plot_genes")
  }
  
  #Colours for different elements
  
  state_colours <- c(Active="#3aab04",
                     Background="#c2c2c2",
                     `Polycomb Repressed`="#d00101",
                     Bivalent="#ff9a01",
                     Unknown="#6e6e6e",
                     Mixed="#eedd9a",
                     `Heterochromatin Repressed`="#b340d5",
                     H3K4me1="#297800")
  
  state_number <- c(Active=1,
                    Background=2,
                    `Polycomb Repressed`=3,
                    Bivalent=4,
                    Unknown=5,
                    Mixed=6,
                    `Heterochromatin Repressed`=7,
                    H3K4me1=8)
  
  #colours from https://www.chicp.org/
  gene_colours <- c(protein_coding="#1f77b4",
                    lincRNA="#ff7f0e",
                    snoRNA="#d62728",
                    antisense="#9467bd",
                    miRNA="#e377c2",
                    snRNA="#8c564b",
                    `processed transcript`="#008000")
  
 raw_state_colours <- c("1"="#b340d5",
                        "2"="#b340d5",
                        "3"="#c2c2c2",
                        "4"="#eedd9a",
                        "5"="#3aab04",
                        "6"="#3aab04",
                        "7"="#3aab04",
                        "8"="#297800",
                        "9"="#297800",
                        "10"="#ff9a01",
                        "11"="#d00101",
                        "12"="#c2c2c2",
                        "13"="#c2c2c2",
                        "14"="#3aab04",
                        "15"="#3aab04",
                        "16"="#3aab04")
  
  
  # state_colours[chromhmm_enh_OSN_lookup$chromhmm_naive[1:10]]
  
  interaction_colours <- c(primed="#C00F14", naive="#3C96EB", naive_primed="#DBDBDB")
  enhancer_colours <- c(SE="#C00F14", E="#669900", B="#DBDBDB")
  
  plot_range <- plotrange(CircosData_object)
  gene_to_plot_df <- genes(CircosData_object)
  chromhmm_to_plot_df <- chromhmm(CircosData_object)
  chromhmm_raw <- chromhmmraw(CircosData_object)
  if(plot_mid){
    b_bed_plot <- bc(CircosData_object)
    oe_bed_plot <- oec(CircosData_object)
  }else{
    b_bed_plot <- b(CircosData_object)
    oe_bed_plot <- oe(CircosData_object)
  }
  
  
  if(plot_type == "enhancer"){
    int_enhs <-  enhancers(CircosData_object)
    int_colours <- unname(enhancer_colours[int_enhs])
  }else if(plot_type %in% c("significant", "full")){
    #colour by origin
    int_orig <- origin(CircosData_object)
    int_colours <- unname(interaction_colours[int_orig])
  }

  #Make circos plot
  
  par(lwd = 0.5)
  circos.par("cell.padding" = c(0, 0, 0, 0))
  circos.genomicInitialize(plot_range, tickLabelsStartFromZero = FALSE, plotType = NULL)
  
  cytoband <- read.cytoband()$df
  circos.genomicTrackPlotRegion(cytoband, ylim = c(0.5, 1), panel.fun = function(region, value, ...) {
    #draw cytoband, a lot of regions will not be plotted for zoomed in
    circos.genomicRect(region, value, col = cytoband.col(value[, 2]), border = NA, ...)
    #min and max of initialized track
    cell.xlim <- get.cell.meta.data("cell.xlim")
    cell.ylim <- get.cell.meta.data("cell.ylim")
    #draw black outline over plotted region
    circos.rect(cell.xlim[1], cell.ylim[1], cell.xlim[2], cell.ylim[2], border = "black")
    #major ticks
    step_size <- round((cell.xlim[2]-cell.xlim[1])/50)
    
    major.at <- seq(cell.xlim[1], cell.xlim[2], by = step_size)
    #convert label to MB
    major.labels <- round(major.at/1000000,2)
    #plot label every 5MB
    l <- (major.at-cell.xlim[1]) %% (step_size*10) == 0
    major.labels[l] = ""
    #plot breaks
    circos.axis("top", major.at = major.at, labels = major.labels,
                labels.facing = "bending.inside", labels.cex = 0.4, 
                major.tick.percentage = 0.5)
    #plot major breaks
    circos.text(major.at[l], rep(1.7, sum(l)), paste0(round(major.at[l]/1000000,2), "MB"), 
                cex = 0.8, facing = "clockwise", adj = c(0, 0.5), niceFacing = TRUE)
  }, bg.border = NA, track.height = 0.02)
  
  circos.par("cell.padding" = c(0.002, 0, 0.002, 0))
  
  #plot genes
  if(plot_genes == "stacked"){
    
    circos.genomicTrackPlotRegion(gene_to_plot_df, ylim=c(-0.5,0.5), panel.fun = function(region, value, ...) {
      for(i in seq_len(nrow(region))) {
        
        #if arrow head length if longer than the region
        if(ux(0.1, "cm") > (region[i, 2]-region[i, 1])){
          ahl <- ux(0.01, "cm")
        }else{
          ahl <- ux(0.1, "cm")
        }
        
        circos.arrow(region[i, 1], region[i, 2], ifelse(is.even(i),0.5, -0.5), width = 0.5,
                     arrow.head.width = 0.8, arrow.head.length = ahl,
                     col = unname(gene_colours[as.character(value[i, 1])]),
                     arrow.position = value[i,3])
        
        circos.genomicText(region[i,], value[i,], y = ifelse(is.even(i),0.5, -0.5), labels.column = 2, facing = "bending.inside",
                           adj = c(0, 0.5), cex = 0.5, posTransform = posTransform.fun,
                           niceFacing = TRUE)
      }
      
    }, bg.border = NA, track.height = 0.05)
    
  }else if(plot_genes == "inline"){
    
    #Plot genes inline
    circos.genomicTrackPlotRegion(gene_to_plot_df, stack = TRUE, panel.fun = function(region, value, ...) {
      for(i in seq_len(nrow(region))) {
        
        #if arrow head length if longer than the region
        if(ux(0.1, "cm") > (region[i, 2]-region[i, 1])){
          ahl <- ux(0.01, "cm")
        }else{
          ahl <- ux(0.1, "cm")
        }
        
        circos.arrow(region[i, 1], region[i, 2], width = 1, 
                     arrow.head.width = 0.8, arrow.head.length = ahl,
                     col = unname(gene_colours[as.character(value[i, 1])]), 
                     arrow.position = value[i,3])
        
        circos.genomicText(region[i,], value[i,], 1, labels.column = 2, facing = "bending.inside", 
                           adj = c(0, 0.5), cex = 0.5, posTransform = posTransform.fun,
                           niceFacing = TRUE)
      }
      
    }, bg.border = NA, track.height = 0.05)
  }
  
  circos.par("cell.padding" = c(0, 0, 0, 0))
  
  if(enhancer_osn){
    
    enh_plot_df <- enhBED(CircosData_object)
    osn_plot_df <- osnBED(CircosData_object)
    
    #plot enhancer OSN tracks
    circos.genomicTrackPlotRegion(osn_plot_df, stack=TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "black", border = "black", ...)
    }, bg.border = NA, track.height = 0.03)
    
    circos.genomicTrackPlotRegion(enh_plot_df, stack=TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = unname(enhancer_colours[as.character(value$type)]), 
                         border = "black", ...)
    }, bg.border = NA, track.height = 0.03)
    
  }
  
  #plot chromhmm track
  if(plot_raw_chromhmm){
    par(lwd = 1) #line thickness for hindiii fragment overlay
    circos.genomicTrackPlotRegion(chromhmm_to_plot_df, stack=TRUE, panel.fun = function(region, value, ...) {
      #plot raw chromhmm track
      circos.genomicRect(chromhmm_raw[c("start","end")], 
                         chromhmm_raw["state"], 
                         col = raw_state_colours[as.character(chromhmm_raw[["state"]])], 
                         border = NA, ...)
      #hindiii fragment overlay
      circos.genomicRect(region, value, col = NA, border = "white", ...)
      
    }, bg.border = NA, track.height = 0.05)
    
  }else{
    #plot hindiii chromhmm track
    circos.genomicTrackPlotRegion(chromhmm_to_plot_df, stack=TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = unname(state_colours[as.character(value$chromhmm)]),
                         border = "black", ...)
    }, bg.border = NA, track.height = 0.05)
  }

  
  if(nrow(b_bed_plot) > 0){
    #plot interactions
    circos.genomicLink(b_bed_plot, oe_bed_plot, col = int_colours,
                       border = "black", lwd = 1)#int_table$score_naive/2)
  }
  
  
  circos.clear()
  

}











#'Merge two vectors into one, collapsing values with ; and ignoring NA
#'
#'@param A first vector
#'@param B second vector
#'
#'@example 
#'A <- c(1, NA)
#'B <- c(NA, 2)
#'
#'@return
#'
#' "1","2"
merge_vect <- function(A, B){
  max_size <- max(length(A), length(B))
  
  sapply(seq(1,max_size), function(i){
    #if index out of bounds an NA is returned
    my_row <- c(A[i], B[i])
    if(all(is.na(my_row))){
      return(NA)
    }else{
      return(paste(my_row[!is.na(my_row)], collapse = ";"))
    }
  })
}




