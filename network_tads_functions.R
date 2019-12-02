

#'Percentage the top hit represents, how many nodes within a single top TAD
#'
#'@param item_vect hits from tad and community overlap
#'
#'@return %
percentage_of_top_hit <- function(item_vect){
  #if multiple tads overlap with a community hindiii they will be joined by ;
  item_vect <- unlist(strsplit(as.character(item_vect), ";"))
  #remove NA values
  item_vect <- item_vect[!is.na(item_vect)]
  if(length(item_vect) == 0){
    return(NA)
  }else{
    item_freq <- table(item_vect)
    pect <- (max(item_freq)/sum(item_freq))*100
    return(pect)
  }
}


#'Find the % of community nodes within the community dominant TAD
#'Can also return the number of TADs a community is overlapping with
#'
#'@param intervals community nodes
#'@param tads_gr GRanges object of TADs
#'@param hindiii_lookup_dt HindIII coordinate lookup data.table
#'@param percentage if TRUE return percentage of overlap, else return number of overlapping TADs
#'
tads_community_overlap <- function(intervals, tads_gr, hindiii_lookup_dt, percentage=TRUE){
  tads_overlaping_com <- vector("numeric", length=length(intervals))
  percentage_hindiii_overlap <- vector("numeric", length=length(intervals))
  for(i in seq_along(intervals)){
    #convert hindiii IDs to granges
    com_gr <- nodes2granges(intervals[i], hindiii_lookup_dt)
    #find overlaps with comminuty nodes and TADs
    overlap <- GenomicRanges::findOverlaps(com_gr, tads_gr)
    #how many TADs is a community overlapping with?
    tads_overlaping_com[i] <- length(unique(subjectHits(overlap)))
    #the % of the top TAD represented within overlaps
    #don't care about nodes not within a TAD
    percentage_hindiii_overlap[i] <- percentage_of_top_hit(subjectHits(overlap))
  }
  if(percentage){
    return(percentage_hindiii_overlap)
  }else{
    return(tads_overlaping_com)
  }
}

#'Evaluate function for permTest 
#'TAD community overlap (tco)
#'calculates the median of all overlaps and can finally plot a boxplot with these?
#'
#'@param B intervals
#'@param A tads_gr
tco_median <- function(A,B, ...){

  percentage_hindiii_overlap <- vector("numeric", length=length(B))
  #find overlaps with comminuty nodes and TADs
  overlap <- GenomicRanges::findOverlaps(B, A)
  
  mcols(B)$tad <- NA
  ovlp_tads <- aggregate(list(tads=subjectHits(overlap)), by=list(com=queryHits(overlap)), 
                         FUN = function(x) paste(x, collapse = ";"))
  mcols(B)$tad[ovlp_tads$com] <- ovlp_tads$tads

  
  #aggregate by communities
  com_agg <- aggregate(data.frame(mcols(B)["tad"]), 
                       by=data.frame(mcols(B)["com"]), 
                       FUN = percentage_of_top_hit)
  
  #the % of the top TAD represented within overlaps
  #don't care about nodes not within a TAD
  
  #Convert NA to 0
  com_agg$tad[is.na(com_agg$tad)] <- 0
  out <- median(com_agg$tad)
  
  return(out)
}


#'Only keep chromosomes listed
#'Taken from:
#'https://support.bioconductor.org/p/83588/
#'
#'@param genome BSgenome object
#'@param seqnames names of chromosomes to keep 
#'
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

