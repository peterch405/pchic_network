

#' From a data.frame with collapsed genes and biotype "," extract a specific biotype
#' 
#' @param input_df Data.frame with gene_name and gene_type columns from which to isolate biotype
#' @param biotype Which biotype to isolate
#' 
#' @example 
#' 
#' filter_biotype(bait_annotation_gene_all, "protein_coding")
#' 
#' @return vector of collapsed gene names of specified biotype
filter_biotype <- function(input_df, biotype){
  biotype_genes <- list()
  for(line in 1:nrow(input_df)){
    genes <- unlist(strsplit(input_df$gene_name[line], ","))
    types <- unlist(strsplit(input_df$gene_type[line], ","))
    
    biotype_gene <- list()
    count <- 1
    for(g in seq_along(types)){
      if(is.na(types[g])){
        next
      }else if(types[g] == biotype){
        biotype_gene[count] <- genes[g]
        count <- count + 1
      }
    }
    #if list is empty
    if(length(biotype_gene) == 0){
      biotype_gene <- NA
    }
    biotype_genes[[line]] <- unique(unlist(biotype_gene))
  }
  
  out <- sapply(biotype_genes, function(x) {
    if(length(x)>1){
      paste(unlist(x), collapse=",")
    }else if(is.na(x)){
      return(NA)
    }else{
      paste(unlist(x), collapse=",")
    }
  })
  
  return(out)
}


#' Coordinates of the promoter upstream 1000bp of TSS
#' 
#' @param gtf_gr GenomicRanges object created from GFF annotation file
#' 
#' @return GenomicRanges object with additional prom_start and prom_end mcols
promoter_start <- function(gtf_gr){
  gtf_gr$prom_start <- 0
  gtf_gr$prom_end <- 0
  
  strand <- as.character(strand(gtf_gr))
  txStart <- start(gtf_gr)
  txEnd <- end(gtf_gr)
  
  gtf_gr$prom_start[strand == "-"] <- txEnd[strand == "-"]
  gtf_gr$prom_end[strand == "-"] <- txEnd[strand == "-"] + 1000
  gtf_gr$prom_start[strand == "+"] <- txStart[strand == "+"] - 1000
  gtf_gr$prom_end[strand == "+"] <- txStart[strand == "+"]
  return(gtf_gr)
}



#'Annotate HindIII fragments with gene promoters
#'
#'@param gtf_gr A GenomicRanges object from GFF file with gene promoter annotation
#'@param fragments_gr A GenomicRanges object with HindIII fragments
#'@param annotate_with Either 'gene' or 'transcript', annotate with promoters based on transcripts or genes
#'@param highlight_baited A GenomicRanges with baited HindIII fragments, if NULL will skip
#'@param keep_meta Keep the metadata within the fragments_gr object
#'@param add_id Add Ensembl gene or transcript ID
#'
#'@return data.frame (1-based coordinates)
annotate_fragments <- function(gtf_gr, fragments_gr, annotate_with=c("gene", "transcript"), 
                               highlight_baited=NULL, keep_meta=FALSE, add_id=FALSE){
  
  stopifnot("gene" %in% annotate_with | "transcript" %in% annotate_with)
  
  if(length(annotate_with)>1){
    print("Annotating with", annotate_with[1])
    annotate_with <- annotate_with[1]
  }

  
  type_col <- paste(annotate_with, "biotype", sep = "_")
  type_name <- paste(annotate_with, "type", sep = "_")
  name_col <- paste(annotate_with, "name", sep = "_")
  name_id <- paste(annotate_with, "id", sep = "_")
  
  hits <- findOverlaps(gtf_gr, fragments_gr)
  
  
  anno_df <- data.frame(gtf_hits=queryHits(hits), 
                        name=mcols(gtf_gr)[queryHits(hits),name_col], 
                        type=mcols(gtf_gr)[queryHits(hits),type_col], 
                        hindiii_hit=subjectHits(hits), 
                        chrom=seqnames(fragments_gr[subjectHits(hits)]),
                        start=start(fragments_gr[subjectHits(hits)]), 
                        end=end(fragments_gr[subjectHits(hits)]))
  
  setDT(anno_df)
  # anno_agg_df <- aggregate(anno_df[name_col], by=anno_df[c("chrom", "start", "end")],
  #                               FUN = function(x) {
  #                                 paste(unique(x), collapse = ",")
  #                               })
  
  #data.table faster than aggregate
  anno_agg_all <- anno_df[,list(name = paste(name, collapse = ",")), 
                          by = "chrom,start,end"]
  
  unq_names <- unname(sapply(anno_agg_all$name, function(x) paste(unique(unlist(strsplit(x, ","))), 
                                                                  collapse = ",")))
  anno_agg_all$name <- unq_names
  
  setnames(anno_agg_all, "name", name_col)
  
  # names(anno_df)[2:3] <- c(name_col, type_name)
  
  gene_type_lookup <- as.character(anno_df[["type"]])
  names(gene_type_lookup) <- as.character(anno_df[["name"]])
  
  anno_agg_all$type <- sapply(strsplit(anno_agg_all[[name_col]], ","), 
                              function(x) paste(unname(gene_type_lookup[unlist(x)]), 
                                                collapse = ","))
  
  setnames(anno_agg_all, "type", type_name)
  
  
  #Mark which hindIII fragments are baited
  if(!is.null(highlight_baited)){
    hits_bait <- findOverlaps(highlight_baited, fragments_gr)
    
    fragments_gr$baited <- 0
    fragments_gr$baited[subjectHits(hits_bait)] <- 1
    
    #add metadata columns
    if(keep_meta){
      frags <- data.frame(fragments_gr)
      names(frags)[1] <- "chrom"
    }else{
      frags <- data.frame(chrom=seqnames(fragments_gr),
                          start=start(fragments_gr),
                          end=end(fragments_gr),
                          baited=mcols(fragments_gr)[["baited"]])
    }
    
  }else{
    #include un-annotated fragments
    frags <- data.frame(chrom=seqnames(fragments_gr),
                        start=start(fragments_gr),
                        end=end(fragments_gr))
    #add metadata columns
    if(keep_meta){
      frags <- cbind(frags, data.frame(mcols(fragments_gr)))
    }
  }
  
  annotation_all <- merge(frags, anno_agg_all, by=c("chrom", "start", "end"), all=TRUE)
  
  #Add ensemble ID
  if(add_id){
    id_lookup <- as.character(mcols(gtf_gr)[[name_id]])
    names(id_lookup) <- as.character(mcols(gtf_gr)[[name_col]])
    
    annotation_all[,name_id] <- sapply(strsplit(annotation_all[[name_col]], ","), 
                                       function(x) paste(unname(id_lookup[unlist(x)]), 
                                                         collapse = ","))
  }
  
  if(is.null(highlight_baited)){
    annotation_out <- annotation_all
  }else{
    #move baited column to last
    annotation_out <- annotation_all[, c(setdiff(names(annotation_all), "baited"), "baited")]
  }
  
  return(annotation_out)
}

