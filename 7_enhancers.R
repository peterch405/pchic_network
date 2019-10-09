library(cowplot)
library(foreach)
library(doSNOW)
library(readr)
library(data.table)
library(ggforce)
library(ggpubr)
library(GenomicRanges)
library(gplots)
library(UpSetR)
library(eulerr)

source("network_OSN_functions.R")
source("network_enhancer_class.R")
source("annotate_promoters_functions.R")
source("network_enhancers_functions.R")


#Determine distance of neighbouring fragments to merge by ROSE -----------------

#https://github.com/taoliu/MACS
#https://charlesjb.github.io/How_to_import_narrowPeak/

#MACS2 bed peaks 0-based; need to add +1 to start position
#https://www.biostars.org/p/84686/
#granges does this automatically and everything is 1-based


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

Naive_H3K27ac_peaks <- rtracklayer::import("7_enhancers/Naive_H3K27ac_peaks_chr.narrowPeak", 
                                           format="bed", extraCols = extraCols_narrowPeak)
Primed_H3K27ac_peaks <- rtracklayer::import("7_enhancers/Primed_H3K27ac_peaks_chr.narrowPeak",
                                            format="bed", extraCols = extraCols_narrowPeak)


naive_dtn <- distanceToNearest(Naive_H3K27ac_peaks)
primed_dtn <- distanceToNearest(Primed_H3K27ac_peaks)
p1 <- qplot(log2(mcols(naive_dtn)$distance), bins = 50, ylim = c(0, 4000)) + 
  scale_x_continuous(breaks = seq(5, 22, 1)) +
  theme_pubr() +
  geom_vline(xintercept = log2(1500), linetype="dashed", 
             color = "grey", size=1)
p2 <- qplot(log2(mcols(primed_dtn)$distance), bins = 50, ylim = c(0, 4000)) + 
  scale_x_continuous(breaks = seq(5, 22, 1)) +
  theme_pubr() +
  geom_vline(xintercept = log2(1500), linetype="dashed", 
             color = "grey", size=1)

cowplot::plot_grid(p1, p2)
ggsave("7_enhancers/ROSE_merge_size.pdf")




#Network enhancer viewpoint ----------------------------------------------------

#load non-trans network
net_all_s <- readRDS("2_network_make/net_all_s_20190829.rds")
chromhmm_enh_OSN_lookup <- readRDS("7_enhancers/chromhmm_enh_OSN_lookup_20190112.rds")

links_nodes_cat_col_coord_deb2b <- readRDS("2_network_make/links_nodes_cat_col_coord_deb2b_20190829.rds")
nodes_all <- links_nodes_cat_col_coord_deb2b$nodes

#Annotated by transcript gene names
nodes_all_prot <- readRDS("7_enhancers/nodes_all_prot.rds")

gene_lookup_dt <- nodes_all_prot[, c("ID","prot_genes"), with=FALSE]

#get hindiii fragments of enhancers and superenhancers
np_enh_OSN_lookup <- readRDS("7_enhancers/np_enh_OSN_lookup_20190709.rds")

nodes_all_prot <- readRDS("7_enhancers/nodes_all_prot.rds")


de_genes <- read_delim("1_DESeq2_gene_expression/de_genes_takashima_GRCh38.87_anno_opposing_strand_prot_genes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
names(de_genes)[8] <- "prot_genes"



#Make enhancer list
naive_enhancers<- aggregate(x = np_enh_OSN_lookup["ID"], 
                            by= np_enh_OSN_lookup[c("name_naive", "type_naive")], 
                            FUN = function(ids){
                              paste(ids, collapse = ",")
                            })

primed_enhancers<- aggregate(x = np_enh_OSN_lookup["ID"], 
                             by= np_enh_OSN_lookup[c("name_primed", "type_primed")], 
                             FUN = function(ids){
                               paste(ids, collapse = ",")
                             })


np_enhancers <- merge(naive_enhancers, primed_enhancers, by="ID", all=TRUE)

np_enhancers$np_name <- apply(np_enhancers, 1, function(x)
  gsub("(^_;|;_$)", "",
       paste(ifelse(is.na(x["name_naive"]), "_", x["name_naive"]), 
             ifelse(is.na(x["name_primed"]), "_", x["name_primed"]), 
             sep = ";"))
)

np_enhancers$np_type <- apply(np_enhancers, 1, function(x)
  gsub("(^_;|;_$)", "",
       paste(ifelse(is.na(x["type_naive"]), "_", x["type_naive"]), 
             ifelse(is.na(x["type_primed"]), "_", x["type_primed"]), 
             sep = ";"))
)


np_enhancer_list <- sapply(np_enhancers$ID, function(x) as.list(strsplit(x, ",")))



#find the interactions with enhancers ------------------------------------------


np_enhancer_summary_cnt <- get_interaction_counts(np_enhancer_list, net_all_s, 
                                                  nodes_all_prot, p=11, 
                                                  count_without_b2b = FALSE, 
                                                  type = "expression", score_diff = 2, 
                                                  gene_lookup_dt)


np_enhancers$ID <- unname(sapply(np_enhancers$ID, function(x) gsub(",", ";", x)))
np_enhancers_cnt <- merge(np_enhancer_summary_cnt, np_enhancers[c("ID", "np_name", "np_type")], by="ID", all = TRUE)


np_enhancers_cnt$prot_genes_naive_srt <- sapply(np_enhancers_cnt$prot_genes_naive, function(x) {
  g_vals <- unlist(strsplit(as.character(x), ";"))
  g_vals <- gsub("\\<NA\\>", NA, g_vals)
  g_vals <- g_vals[!is.na(g_vals)]
  if(length(g_vals) == 0){
    return(NA)
  }else{
    return(paste(g_vals, collapse = ";"))
  }
})

np_enhancers_cnt$prot_genes_primed_srt <- sapply(np_enhancers_cnt$prot_genes_primed, function(x) {
  g_vals <- unlist(strsplit(as.character(x), ";"))
  g_vals <- gsub("\\<NA\\>", NA, g_vals)
  g_vals <- g_vals[!is.na(g_vals)]
  if(length(g_vals) == 0){
    return(NA)
  }else{
    return(paste(g_vals, collapse = ";"))
  }
})

#column with shared interacting genes

np_enhancers_cnt$prot_genes_shared <- apply(np_enhancers_cnt, 1, function(x){
  n_genes <- unlist(strsplit(as.character(x["prot_genes_naive_srt"]), ";"))
  p_genes <- unlist(strsplit(as.character(x["prot_genes_primed_srt"]), ";"))
  out <- unique(n_genes[n_genes %in% p_genes])
  if(length(out[!is.na(out)]) == 0){
    return(NA)
  }else{
    return(paste(out, collapse = ";"))
  }
})

np_enhancers_cnt$prot_genes_naive_only <- apply(np_enhancers_cnt, 1, function(x){
  n_genes <- unlist(strsplit(as.character(x["prot_genes_naive_srt"]), ";"))
  p_genes <- unlist(strsplit(as.character(x["prot_genes_primed_srt"]), ";"))
  out <- unique(n_genes[!n_genes %in% p_genes])
  if(length(out[!is.na(out)]) == 0){
    return(NA)
  }else{
    return(paste(out, collapse = ";"))
  }
})


np_enhancers_cnt$prot_genes_primed_only <- apply(np_enhancers_cnt, 1, function(x){
  n_genes <- unlist(strsplit(as.character(x["prot_genes_naive_srt"]), ";"))
  p_genes <- unlist(strsplit(as.character(x["prot_genes_primed_srt"]), ";"))
  out <- unique(p_genes[!p_genes %in% n_genes])
  if(length(out[!is.na(out)]) == 0){
    return(NA)
  }else{
    return(paste(out, collapse = ";"))
  }
})

#gene of the enhancer viewpoint
np_enhancers_cnt$enh_prot_genes <- apply(np_enhancers_cnt, 1, function(x) {
  fnd_genes <- gene_lookup_dt[unlist(strsplit(as.character(x["ID"]), ";"))]$prot_genes
  out <- fnd_genes[!is.na(fnd_genes)]
  if(length(out) == 0){
    return(NA)
  }else{
    return(paste(out, collapse = ";"))
  }
})


np_enhancers_cnt$n_enh_type <- apply(np_enhancers_cnt, 1, function(x) {
  enh <- unlist(strsplit(as.character(x["np_name"]), ";"))
  enh_type <- unlist(strsplit(as.character(x["np_type"]), ";"))
  
  out <- enh_type[grepl("np", enh)]
  if(length(out) == 0){
    return(NA)
  }else{
    return(paste(out, collapse = ";"))
  }
})

np_enhancers_cnt$p_enh_type <- apply(np_enhancers_cnt, 1, function(x) {
  enh <- unlist(strsplit(as.character(x["np_name"]), ";"))
  enh_type <- unlist(strsplit(as.character(x["np_type"]), ";"))
  
  out <- enh_type[grepl("pp", enh)]
  if(length(out) == 0){
    return(NA)
  }else{
    return(paste(out, collapse = ";"))
  }
})



#seperate into SE and E

np_se_cnt <- np_enhancers_cnt[grepl("SE",np_enhancers_cnt$np_type),]

np_e_cnt <- np_enhancers_cnt[sapply(np_enhancers_cnt$np_type, function(x) "E" %in% unlist(strsplit(as.character(x), ";"))),]


#Make Venn diagram tables and plot boxplots for expression----------------------

enhs <- strsplit(np_enhancers_cnt$np_type, ";")
naive_only <- sapply(strsplit(np_enhancers_cnt$np_name, ";"), function(x) grepl("np", x))
primed_only <-  sapply(strsplit(np_enhancers_cnt$np_name, ";"), function(x) grepl("pp", x))

naive_list <- list()
for(i in seq_along(enhs)){
  naive_list[[i]] <- enhs[[i]][naive_only[[i]]]
}

primed_list <- list()
for(i in seq_along(enhs)){
  primed_list[[i]] <- enhs[[i]][primed_only[[i]]]
}



np_enhancers_cnt$naive_se <- sapply(naive_list, function(x) any(x=="SE"))

np_enhancers_cnt$primed_se <- sapply(primed_list, function(x) any(x=="SE"))

np_enhancers_cnt$naive_e <- sapply(naive_list, function(x) any(x=="E"))

np_enhancers_cnt$primed_e <- sapply(primed_list, function(x) any(x=="E"))


#Make new table with SE genes only and get expression ------------------------

SE_naive_prot_genes <- c(np_enhancers_cnt$prot_genes_naive_only[np_enhancers_cnt$naive_se], 
                         np_enhancers_cnt$prot_genes_shared[np_enhancers_cnt$naive_se])
SE_naive_prot_genes <- unique(unlist(strsplit(SE_naive_prot_genes, ";")))
SE_naive_prot_genes <- SE_naive_prot_genes[!(SE_naive_prot_genes == "NA")]
SE_naive_prot_genes <- SE_naive_prot_genes[!is.na(SE_naive_prot_genes)]

SE_primed_prot_genes <- c(np_enhancers_cnt$prot_genes_primed_only[np_enhancers_cnt$primed_se], 
                          np_enhancers_cnt$prot_genes_shared[np_enhancers_cnt$primed_se])
SE_primed_prot_genes <- unique(unlist(strsplit(SE_primed_prot_genes, ";")))
SE_primed_prot_genes <- SE_primed_prot_genes[!(SE_primed_prot_genes == "NA")]
SE_primed_prot_genes <- SE_primed_prot_genes[!is.na(SE_primed_prot_genes)]


SE_shared_prot_genes <- intersect(SE_naive_prot_genes, SE_primed_prot_genes)

SE_naive_prot_genes <- SE_naive_prot_genes[!(SE_naive_prot_genes %in% SE_shared_prot_genes)]

SE_primed_prot_genes <- SE_primed_prot_genes[!(SE_primed_prot_genes %in% SE_shared_prot_genes)]

SE_df <- data.frame(prot_genes=c(SE_shared_prot_genes, SE_naive_prot_genes, SE_primed_prot_genes), 
                    origin=(c(rep("shared", length(SE_shared_prot_genes)), 
                              rep("naive", length(SE_naive_prot_genes)), 
                              rep("primed", length(SE_primed_prot_genes)))))

#Add expression
SE_df_expr <- merge(SE_df, de_genes, by="prot_genes", all.x = TRUE)
SE_df_expr_plot <- SE_df_expr[c("prot_genes", "origin", "Naive_mean_fpkm_log2", "Primed_mean_fpkm_log2")]
SE_plot <- reshape2::melt(SE_df_expr_plot)


ggpaired(SE_plot, x = "variable", y = "value",facet.by = "origin",
         color = "variable", ylab="log2 FPKM", xlab=NULL,
         line.color = "grey", line.size = 0, label = NULL, point.size = 0,
         font.label = list(size = 5, color = "black"), width=0.2) +
  geom_violin() +
  stat_compare_means(paired = TRUE) 
ggsave("7_enhancers/SE_expression.pdf", width = 5, height = 5)

#Do the same for enhancers -----------------------------------------------------
#Remove SE from E list

E_naive_prot_genes <- c(np_enhancers_cnt$prot_genes_naive_only[np_enhancers_cnt$naive_e], np_enhancers_cnt$prot_genes_shared[np_enhancers_cnt$naive_e])
E_naive_prot_genes <- unique(unlist(strsplit(E_naive_prot_genes, ";")))
E_naive_prot_genes <- E_naive_prot_genes[!(E_naive_prot_genes == "NA")]
E_naive_prot_genes <- E_naive_prot_genes[!is.na(E_naive_prot_genes)]

E_primed_prot_genes <- c(np_enhancers_cnt$prot_genes_primed_only[np_enhancers_cnt$primed_e], np_enhancers_cnt$prot_genes_shared[np_enhancers_cnt$primed_e])
E_primed_prot_genes <- unique(unlist(strsplit(E_primed_prot_genes, ";")))
E_primed_prot_genes <- E_primed_prot_genes[!(E_primed_prot_genes == "NA")]
E_primed_prot_genes <- E_primed_prot_genes[!is.na(E_primed_prot_genes)]

#remove E that have SE's
remove_naive_from_E <- intersect(SE_naive_prot_genes, E_naive_prot_genes)
remove_primed_from_E <- intersect(SE_primed_prot_genes, E_primed_prot_genes)



E_shared_prot_genes <- intersect(E_naive_prot_genes, E_primed_prot_genes)
remove_shared_from_E <- intersect(SE_shared_prot_genes, c(E_shared_prot_genes,E_naive_prot_genes, E_primed_prot_genes))


E_naive_prot_genes_woSE <- E_naive_prot_genes[!(E_naive_prot_genes %in% c(remove_naive_from_E, remove_shared_from_E))]
E_primed_prot_genes_woSE <- E_primed_prot_genes[!(E_primed_prot_genes %in% c(remove_primed_from_E, remove_shared_from_E))]

E_shared_prot_genes_woSE <- E_shared_prot_genes[!(E_shared_prot_genes %in% c(remove_naive_from_E, remove_primed_from_E,remove_shared_from_E))]

E_naive_prot_genes_woSE <- E_naive_prot_genes_woSE[!(E_naive_prot_genes_woSE %in% E_shared_prot_genes)]

E_primed_prot_genes_woSE <- E_primed_prot_genes_woSE[!(E_primed_prot_genes_woSE %in% E_shared_prot_genes)]

E_df <- data.frame(prot_genes=c(E_shared_prot_genes_woSE, E_naive_prot_genes_woSE, E_primed_prot_genes_woSE), 
                   origin=(c(rep("shared", length(E_shared_prot_genes_woSE)), 
                             rep("naive", length(E_naive_prot_genes_woSE)), 
                             rep("primed", length(E_primed_prot_genes_woSE)))))

#Add expression
E_df_expr <- merge(E_df, de_genes, by="prot_genes", all.x = TRUE)
E_df_expr_plot <- E_df_expr[c("prot_genes", "origin", "Naive_mean_fpkm_log2", "Primed_mean_fpkm_log2")]
E_plot <- reshape2::melt(E_df_expr_plot)


ggpaired(E_plot, x = "variable", y = "value",facet.by = "origin",
         color = "variable", ylab="log2 FPKM", xlab=NULL,
         line.color = "grey", line.size = 0, label = NULL, point.size = 0,
         font.label = list(size = 5, color = "black"), width=0.2) +
  geom_violin() +
  stat_compare_means(paired = TRUE) 
ggsave("7_enhancers/E_expression.pdf", width = 5, height = 5)




write.csv(SE_df_expr_plot, "7_enhancers/SE_df.csv")
write.csv(E_df_expr_plot, "7_enhancers/E_df.csv")


#Plot venn diagrams ------------------------------------------------------------

pdf("7_enhancers/SE_E_venns.pdf", width = 4, height = 4)
e_fit <- euler(c(A = sum(E_df_expr_plot$origin == "naive"), 
                 B = sum(E_df_expr_plot$origin == "primed"), 
                 "A&B" = sum(E_df_expr_plot$origin == "shared")))
plot(e_fit, main="Enhancer", labels = c("Naive", "Primed", "Shared"),
     quantities = list(type = c("counts")))
grid.newpage()
se_fit <- euler(c(A = sum(SE_df_expr_plot$origin == "naive"), 
                  B = sum(SE_df_expr_plot$origin == "primed"), 
                  "A&B" = sum(SE_df_expr_plot$origin == "shared")))
plot(se_fit, main="SuperEnhancer", labels = c("Naive", "Primed", "Shared"),
     quantities = list(type = c("counts")))
dev.off()





#Load ROSE bed -----------------------------------------------------------------
hindiii_gr <- readRDS("6_chromhmm_plots/hindiii_gr.rds")

Naive_enh <- rtracklayer::import("7_enhancers/Naive_H3K27ac_peaks_narrowPeak_Gateway_Enhancers.bed",
                                 format="bed")
Naive_suprenh <- rtracklayer::import("7_enhancers/Naive_H3K27ac_peaks_narrowPeak_Gateway_SuperEnhancers.bed",
                                     format="bed")

mcols(Naive_enh)$type <- "E"
mcols(Naive_enh)$type[mcols(Naive_enh)$name %in% mcols(Naive_suprenh)$name] <- "SE"

Primed_enh <- rtracklayer::import("7_enhancers/Primed_H3K27ac_peaks_narrowPeak_Gateway_Enhancers.bed",
                                  format="bed")
Primed_suprenh <- rtracklayer::import("7_enhancers/Primed_H3K27ac_peaks_narrowPeak_Gateway_SuperEnhancers.bed",
                                      format="bed")

mcols(Primed_enh)$type <- "E"
mcols(Primed_enh)$type[mcols(Primed_enh)$name %in% mcols(Primed_suprenh)$name] <- "SE"



#Overlap with OSN peaks (don't do this on a hindiii fragment, but on a peak resolution)
# 1 = naive; 2 = shared; 3 = primed
OSN_all <- rtracklayer::import("7_enhancers/osn_num.bed",
                               format = "bed")

OSN_naive <- OSN_all[mcols(OSN_all)$name %in% c(1, 2)]
OSN_primed <- OSN_all[mcols(OSN_all)$name %in% c(2, 3)]


Primed_enh_hindiii_OSN <- ROSE_OSN_HindIII_table(Primed_enh, hindiii_gr, OSN_primed)
Naive_enh_hindiii_OSN <- ROSE_OSN_HindIII_table(Naive_enh, hindiii_gr, OSN_naive)

#enhancer lookup for OSN (has both naive and primed OSN overlaping naive primed enhancers)

enh_OSN_ovl <- findOverlaps(c(Primed_enh, Naive_enh), OSN_all)

OSN_enhancer_df <- data.frame(name=c(Primed_enh, Naive_enh)[queryHits(enh_OSN_ovl)]$name,
                              OSN=OSN_all[subjectHits(enh_OSN_ovl)]$name)

OSN_enhancer_lookup <- aggregate(x = OSN_enhancer_df["OSN"],
                                 by= OSN_enhancer_df["name"],
                                 FUN = function(ids){
                                   paste(ids, collapse = ",")
                                 })

OSN_enhancer_lookup$origin <- NA
OSN_enhancer_lookup$origin[grepl("naive", OSN_enhancer_lookup$name)] <- "naive"
OSN_enhancer_lookup$origin[grepl("Primed", OSN_enhancer_lookup$name)] <- "primed"

OSN_enhancer_lookup$short_name <- ifelse(grepl("naive", OSN_enhancer_lookup$name),
                                         gsub("_naive_peak", "_np", gsub("_lociStitched", "", OSN_enhancer_lookup$name)),
                                         gsub("_Primed_H3K27ac_peak", "_pp", gsub("_lociStitched", "", OSN_enhancer_lookup$name)))


setDT(OSN_enhancer_lookup, key="short_name")


#LogOdds ratio -----------------------------------------------------------------
#get interaction counts for all baits

bait_annotation_gene_all <- readRDS("7_enhancers/baited_g38_promoter_HindIII_fragments_gene_anno_20181207.rds")

hindiii_lookup <- read_delim("7_enhancers/hindiii_lookup.tab", 
                             "\t", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE)
names(hindiii_lookup) <- c("chrom", "start", "end", "ID")

#only interested in protein coding
bait_annotation_gene_all$prot_genes <- filter_biotype(bait_annotation_gene_all, 
                                                      "protein_coding")
#add hindiii ID
bait_annotation_gene_all$start <- bait_annotation_gene_all$start-1

baited_genes_hindiii <- merge(bait_annotation_gene_all, hindiii_lookup, 
                              by=c("chrom", "start", "end"), all.x=TRUE)


mcols(Naive_enh)$short_name <- gsub("_lociStitched", "", gsub("naive_peak", "np", Naive_enh$name))
mcols(Primed_enh)$short_name <- gsub("_lociStitched", "", gsub("Primed_H3K27ac_peak", "pp", Primed_enh$name))


baited_enhancer_summary_cnt <- get_interaction_counts(baited_genes_hindiii$ID, net_all_s, 
                                                      chromhmm_enh_OSN_lookup, p=7, 
                                                      count_without_b2b = FALSE, 
                                                      type = "raw", score_diff = 2, 
                                                      gene_lookup_dt)

#Simplify enhancer tables
np_seperate_enhancer <- list()
for(i in seq_along(baited_enhancer_summary_cnt)){
  
  #count two or more enhancers within interacting hindiii fragment seperatly
  
  n_enh <- enhancers(baited_enhancer_summary_cnt[[i]], "naive")
  p_enh <- enhancers(baited_enhancer_summary_cnt[[i]], "primed")
  
  
  
  
  if(nrow(n_enh) > 0 | nrow(p_enh) > 0){
    #add baited ID
    b_ID <- node(baited_enhancer_summary_cnt[[i]])
    np_enhs <- merge(n_enh, p_enh, by="ID", all=TRUE)
    np_enhs$b_ID <- b_ID
    
  }else{
    #if no enhancers present
    next 
  }
  
  # resolve multiple enhancers in a single hindiii fragment
  resolve <- grepl(";", np_enhs$name_naive) | grepl(";", np_enhs$name_primed)
  
  if(sum(resolve*1) > 0){
    
    np_enhs_resolved <- apply(np_enhs[resolve,], 1, function(x){
      
      #if both sides have multiple enhancers then resolve, else only split
      if(grepl(";", x[["name_naive"]]) & grepl(";",x[["name_primed"]])){
        n_x <- unlist(strsplit(x[["name_naive"]], ";"))
        p_x <- unlist(strsplit(x[["name_primed"]], ";"))
        
        
        
        naive_split <- data.frame(as.list(x[c("name_naive", "score_naive", "type_naive", "OSN_naive")])) %>% 
          separate_rows(name_naive, score_naive, type_naive, OSN_naive, sep=";")
        
        primed_split <- data.frame(as.list(x[c("name_primed", "score_primed", "type_primed", "OSN_primed")])) %>% 
          separate_rows(name_primed, score_primed, type_primed, OSN_primed, sep=";")
        
        #resolve which overlapping enhancers should be put together and which should be 
        #left with NA on other end
        ovlps <- findOverlaps(Naive_enh[Naive_enh$short_name %in% n_x],
                              Primed_enh[Primed_enh$short_name %in% p_x])
        
        
        naive_vect <- c(n_x[queryHits(ovlps)], 
                        setdiff(n_x, n_x[queryHits(ovlps)]), 
                        rep(NA, length(setdiff(p_x, p_x[subjectHits(ovlps)]))))
        
        primed_vect <- c(p_x[subjectHits(ovlps)], 
                         rep(NA, length(setdiff(n_x, n_x[queryHits(ovlps)]))),
                         setdiff(p_x, p_x[subjectHits(ovlps)]))
        
        stopifnot(length(naive_vect) == length(primed_vect))
        
        split_list <- list()
        #generate all the rows
        for(i in 1:length(naive_vect)){
          n_data_pos <- which(naive_vect[i] == naive_split$name_naive)
          p_data_pos <- which(primed_vect[i] == primed_split$name_primed)
          #generate an entire row, substituting name_naive name_primed and filling 
          #the rest from naive_split or primed_split or filling with NA
          split_list[[i]] <- data.frame(as.list(c(x[c("ID", "chromhmm_naive")], 
                                                  name_naive=naive_vect[i],
                                                  if(is.na(naive_vect[i])) 
                                                    c(score_naive=NA, type_naive=NA, OSN_naive=NA)  
                                                  else 
                                                    naive_split[n_data_pos,c("score_naive","type_naive", "OSN_naive")],
                                                  x["chromhmm_primed"],
                                                  name_primed=primed_vect[i], 
                                                  if(is.na(primed_vect[i])) 
                                                    c(score_primed=NA, type_primed=NA, OSN_primed=NA) 
                                                  else 
                                                    primed_split[p_data_pos,c("score_primed", "type_primed", "OSN_primed")],
                                                  x["b_ID"])))
        }
        #merge all lists into df
        resolved <- plyr::ldply(split_list, rbind)
        #convert to real NA
        resolved[resolved == "NA"] <- NA
        return(resolved)
        
        #if ; only on one side split that side, the other side will be filled with NA
      }else{
        if(grepl(";", x[["name_naive"]])){
          #need to convert named vector to data.frame for separate_rows 
          #(not sure if there is a method for named vectors)
          resolved <- data.frame(as.list(x)) %>% 
            separate_rows(name_naive, score_naive, type_naive, OSN_naive, sep=";")
        }
        if(grepl(";",x[["name_primed"]])){
          resolved <- data.frame(as.list(x)) %>% 
            separate_rows(name_primed, score_primed, type_primed, OSN_primed, sep=";")
        }
        #convert to real NA
        resolved[resolved == "NA"] <- NA
        return(resolved)
      }
    })
    
    np_enhs_out <- plyr::ldply(np_enhs_resolved, rbind)
    
    none_to_resolve <- FALSE
  }else{
    none_to_resolve <- TRUE
  }
  
  
  if(none_to_resolve){
    out <- np_enhs
  }else{
    out <- rbind(np_enhs[!resolve,], np_enhs_out)
  }
  
  np_seperate_enhancer[[i]] <- out
}


interacting_enhancer <- plyr::ldply(np_seperate_enhancer, rbind)

#gained enhancer interactions (primed)
gained <- interacting_enhancer[is.na(interacting_enhancer$name_naive),]
#add the OSN of the other cell type (despite loss of enhancer, these OSN sites are gained)
gained_both <- merge(gained, data.frame(name_primed=OSN_enhancer_lookup$short_name, OSN=OSN_enhancer_lookup$OSN), 
                     by="name_primed", all.x=TRUE)
gained_both$OSN_naive <- gained_both$OSN
gained_both$OSN_naive[!is.na(gained_both$OSN_primed)] <- NA
gained_both <- gained_both[,-which(names(gained_both) %in% "OSN")]


#lost enhancer interactions (naive)
lost <- interacting_enhancer[is.na(interacting_enhancer$name_primed),]
#add the OSN of the other cell type
lost_both <- merge(lost, data.frame(name_naive=OSN_enhancer_lookup$short_name, OSN=OSN_enhancer_lookup$OSN), 
                   by="name_naive", all.x=TRUE)
lost_both$OSN_primed <-  lost_both$OSN
lost_both$OSN_primed[!is.na(lost_both$OSN_naive)] <- NA
lost_both <- lost_both[,-which(names(lost_both) %in% "OSN")]

retained <- interacting_enhancer[!is.na(interacting_enhancer$name_primed) &
                                   !is.na(interacting_enhancer$name_naive),]
retained[retained == "NA"] <- NA
#if the two enhancers overlap, count as retained, otherwise count as lost and gained taking OSN into consideration


lost[lost == "NA"] <- NA
gained[gained == "NA"] <- NA

gained_osn_counts <- count_osn_gain_loss(gained_both$OSN_primed, as.character(gained_both$OSN_naive), "gained")
lost_osn_counts <- count_osn_gain_loss(lost_both$OSN_naive, as.character(lost_both$OSN_primed), "lost")


gained_osn_groups <- count_osn_gain_loss(gained_both$OSN_primed, as.character(gained_both$OSN_naive), "gained", return_vector = TRUE)
gained_both$groups <- gained_osn_groups

lost_osn_groups <- count_osn_gain_loss(lost_both$OSN_naive, as.character(lost_both$OSN_primed), "lost", return_vector = TRUE)
lost_both$groups <- lost_osn_groups



retained_osn_counts <- table(unlist(count_retained(retained, Naive_enh, Primed_enh, OSN_enhancer_lookup)))
retained_osn_groups  <- count_retained(retained, Naive_enh, Primed_enh, OSN_enhancer_lookup)

retained$groups <- sapply(retained_osn_groups, paste, collapse = ";")

#' incremental operator
#' https://stackoverflow.com/questions/5738831/r-plus-equals-and-plus-plus-equivalent-from-c-c-java-etc
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

counts_matrix <- matrix(0, nrow=5, ncol=3, 
                        dimnames = list(c("OSN_lost", "OSN_gained", "OSN_retained", "OSN_none", "OSN_mixed"),
                                        c("enhancer_interaction_lost", "enhancer_interaction_gained", "enhancer_interaction_retained")))

counts_matrix["OSN_gained","enhancer_interaction_gained"] <- gained_osn_counts$gained
counts_matrix["OSN_lost","enhancer_interaction_gained"] <- gained_osn_counts$lost
counts_matrix["OSN_retained","enhancer_interaction_gained"] <- gained_osn_counts$retained
counts_matrix["OSN_none","enhancer_interaction_gained"] <- gained_osn_counts$none


counts_matrix["OSN_gained","enhancer_interaction_lost"] <- lost_osn_counts$gained
counts_matrix["OSN_lost","enhancer_interaction_lost"] <- lost_osn_counts$lost
counts_matrix["OSN_retained","enhancer_interaction_lost"] <- lost_osn_counts$retained
counts_matrix["OSN_none","enhancer_interaction_lost"] <- lost_osn_counts$none

counts_matrix["OSN_gained", "enhancer_interaction_lost"] %+=% unname(retained_osn_counts["enh_lost_osn_gained"])
counts_matrix["OSN_lost", "enhancer_interaction_lost"] %+=% unname(retained_osn_counts["enh_lost_osn_lost"])
counts_matrix["OSN_none", "enhancer_interaction_lost"] %+=% unname(retained_osn_counts["enh_lost_osn_none"])
counts_matrix["OSN_retained", "enhancer_interaction_lost"] %+=% unname(retained_osn_counts["enh_lost_osn_retained"])
counts_matrix["OSN_mixed", "enhancer_interaction_lost"] %+=% unname(ifelse(is.na(retained_osn_counts["enh_lost_osn_mixed"]), 0, 
                                                                           unlist(retained_osn_counts["enh_lost_osn_mixed"])))

counts_matrix["OSN_gained", "enhancer_interaction_gained"] %+=% unname(retained_osn_counts["enh_gained_osn_gained"])
counts_matrix["OSN_lost", "enhancer_interaction_gained"] %+=% unname(retained_osn_counts["enh_gained_osn_lost"])
counts_matrix["OSN_none", "enhancer_interaction_gained"] %+=% unname(retained_osn_counts["enh_gained_osn_none"])
counts_matrix["OSN_retained", "enhancer_interaction_gained"] %+=% unname(retained_osn_counts["enh_gained_osn_retained"])
counts_matrix["OSN_mixed", "enhancer_interaction_gained"] %+=% unname(ifelse(is.na(retained_osn_counts["enh_gained_osn_mixed"]), 0, 
                                                                             unlist(retained_osn_counts["enh_gained_osn_mixed"])))

counts_matrix["OSN_gained", "enhancer_interaction_retained"] %+=% unname(retained_osn_counts["enh_retained_osn_gained"])
counts_matrix["OSN_lost", "enhancer_interaction_retained"] %+=% unname(retained_osn_counts["enh_retained_osn_lost"])
counts_matrix["OSN_none", "enhancer_interaction_retained"] %+=% unname(retained_osn_counts["enh_retained_osn_none"])
counts_matrix["OSN_retained", "enhancer_interaction_retained"] %+=% unname(retained_osn_counts["enh_retained_osn_retained"])
counts_matrix["OSN_mixed", "enhancer_interaction_retained"] %+=% unname(ifelse(is.na(retained_osn_counts["enh_retained_osn_mixed"]), 0, 
                                                                               unlist(retained_osn_counts["enh_retained_osn_mixed"])))



pdf("7_enhancers/log2oddsratio_enhancer_interaction.pdf", height = 6, width = 6)
res3 <- logOddsRatio(counts_matrix)
grid.newpage()
gridExtra::grid.table(data.frame(counts_matrix), theme=gridExtra::ttheme_default(base_size=7))
res <- logOddsRatio(counts_matrix[1:4,1:3])
dev.off()
  

#Barplots / Upset plots --------------------------------------------------------
#only interacting enhancers

links_nodes_cat_col_coord_deb2b <- readRDS("/media/chovanec/My_Passport/CHiC_naive_primed/network/network_analysis/links_nodes_cat_col_coord_deb2b_20190829.rds")
links_all <- links_nodes_cat_col_coord_deb2b$links
nodes_all <- links_nodes_cat_col_coord_deb2b$nodes

interacting_nodes <- unique(as.character(nodes_all$ID))

naive_interacting <- unique(c(links_all$b_ID[links_all$origin == "naive"], 
                              links_all$oe_ID[links_all$origin == "naive"]))
primed_interacting <- unique(c(links_all$b_ID[links_all$origin == "primed"], 
                               links_all$oe_ID[links_all$origin == "primed"]))



#need coordinates to overlap with enhancers
interacting_hindiii_gr <- hindiii_gr[hindiii_gr$ID %in% interacting_nodes]

interacting_naive_enh_ovlp <- findOverlaps(Naive_enh, interacting_hindiii_gr)
interacting_primed_enh_ovlp <- findOverlaps(Primed_enh, interacting_hindiii_gr)
Naive_enh_int <- Naive_enh[queryHits(interacting_naive_enh_ovlp)]
Primed_enh_int <- Primed_enh[queryHits(interacting_primed_enh_ovlp)]



expressionInput_int_se <- c(Naive = table(Naive_enh_int $ID[Naive_enh_int$type == "SE"])[["n"]], 
                            Primed = table(Primed_enh_int$ID[Primed_enh_int$type == "SE"])[["p"]], 
                            Shared = 0,
                            `Naive&Shared` = table(Naive_enh_int $ID[Naive_enh_int$type == "SE"])[["s"]], 
                            `Primed&Shared` = table(Primed_enh_int$ID[Primed_enh_int$type == "SE"])[["s"]])

plot_int_se_df <- fromExpression(expressionInput_int_se)
pdf("7_enhancers/SE_np_overlap_interacting_only.pdf", 
    onefile=FALSE, width = 4, height = 4)
upset(plot_int_se_df, nsets = 3, keep.order = TRUE, #order.by = "freq",
      point.size = 2.5, line.size = 1.5, 
      mainbar.y.label = "SE Intersections")
dev.off()

expressionInput_int_enh <- c(Naive = sum((mcols(Naive_enh_int)$ID == "n")*1), 
                             Primed = sum((mcols(Primed_enh_int)$ID == "p")*1), 
                             Shared = 0,
                             `Naive&Shared` = sum((mcols(Naive_enh_int)$ID == "s")*1), 
                             `Primed&Shared` = sum((mcols(Primed_enh_int)$ID == "s")*1))

plot_int_enh_df <- fromExpression(expressionInput_int_enh)

pdf("7_enhancers/enhancer_overlap_interacting_only.pdf", 
    onefile=FALSE, width = 4, height = 4)
upset(plot_int_enh_df, nsets = 6, keep.order = TRUE,order.by = "freq",
      point.size = 2.5, line.size = 1.5, 
      mainbar.y.label = "Enhancer Intersections")
dev.off()







