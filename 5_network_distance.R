library(readr)
library(tidyr)
library(plyr)
library(ComplexHeatmap)
library(ggplot2)

#interactions with chip-seq data -----------------------------------------------

naive_intract <- read_delim("0_network_data_preparation/naive_intract_wchip_20200911.txt", 
                            "\t", escape_double = FALSE, col_names = TRUE, 
                            trim_ws = TRUE)

primed_intract <- read_delim("0_network_data_preparation/primed_intract_wchip_20200911.txt", 
                             "\t", escape_double = FALSE, col_names = TRUE, 
                             trim_ws = TRUE)

#need to fill chip-seq information from nodes not from interactions as nodes may still have mark even without interaction!

naive_nodes <- read_delim("0_network_data_preparation/naive_nodes_wchip_20200911.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

primed_nodes <- read_delim("0_network_data_preparation/primed_nodes_wchip_20200911.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)



#Heatmap of 1000 interactions by distance---------------------------------------


primed_fheat <- primed_intract[c("b_ID", "oe_ID", "distance")]
naive_fheat <- naive_intract[c("b_ID", "oe_ID", "distance")]


primed_fheat$ID <- paste(primed_fheat$b_ID, primed_fheat$oe_ID, sep = ":")
naive_fheat$ID <- paste(naive_fheat$b_ID, naive_fheat$oe_ID, sep = ":")


np_fheat <- merge(primed_fheat, naive_fheat, by="ID", all=TRUE, suffixes = c("_primed","_naive"))

np_fheat <- np_fheat[c("ID", "distance_primed", "distance_naive")]
#seperate ID into seperate columns
np_fheat <- np_fheat %>% separate("ID", c("b_ID", "oe_ID"),  sep = ":")



primed_nodes_sub <- primed_nodes[c("ID", "chromhmm", "compartment")]

naive_nodes_sub <- naive_nodes[c("ID", "chromhmm", "compartment")]

#primed b
orig_names <- colnames(primed_nodes_sub)
colnames(primed_nodes_sub)[1] <- paste("b", colnames(primed_nodes_sub)[1], sep = "_")
colnames(primed_nodes_sub)[2:length(orig_names)] <- paste(colnames(primed_nodes_sub)[2:length(orig_names)], "b_primed", sep = "_")
np_fheat <- merge(np_fheat, primed_nodes_sub, by="b_ID", all.x=TRUE)
#primed oe
colnames(primed_nodes_sub) <- orig_names
colnames(primed_nodes_sub)[1] <- paste("oe", colnames(primed_nodes_sub)[1], sep = "_")
colnames(primed_nodes_sub)[2:length(orig_names)] <- paste(colnames(primed_nodes_sub)[2:length(orig_names)], "oe_primed", sep = "_")
np_fheat <- merge(np_fheat, primed_nodes_sub, by="oe_ID", all.x=TRUE)


#naive b
orig_names <- colnames(naive_nodes_sub)
colnames(naive_nodes_sub)[1] <- paste("b", colnames(naive_nodes_sub)[1], sep = "_")
colnames(naive_nodes_sub)[2:length(orig_names)] <- paste(colnames(naive_nodes_sub)[2:length(orig_names)], "b_naive", sep = "_")
np_fheat <- merge(np_fheat, naive_nodes_sub, by="b_ID", all.x=TRUE)
#naive oe
colnames(naive_nodes_sub) <- orig_names
colnames(naive_nodes_sub)[1] <- paste("oe", colnames(naive_nodes_sub)[1], sep = "_")
colnames(naive_nodes_sub)[2:length(orig_names)] <- paste(colnames(naive_nodes_sub)[2:length(orig_names)], "oe_naive", sep = "_")
np_fheat <- merge(np_fheat, naive_nodes_sub, by="oe_ID", all.x=TRUE)

#Make clustering possible ------------------------------------------------------

np_fheat$chromhmm_b_naive_name <- np_fheat$chromhmm_b_naive
np_fheat$chromhmm_oe_naive_name <- np_fheat$chromhmm_oe_naive
np_fheat$chromhmm_b_primed_name <- np_fheat$chromhmm_b_primed
np_fheat$chromhmm_oe_primed_name <- np_fheat$chromhmm_oe_primed



chromhmm_state2_num <- function(column){
  column %<>%
    gsub("Active", 1, .) %>%
    gsub("Background", 2, .) %>%
    gsub("Bivalent", 3, .) %>%
    gsub("Heterochromatin Repressed", 4, .) %>%
    gsub("Mixed", 5, .) %>%
    gsub("Polycomb Repressed", 6, .) %>%
    gsub("Unknown", 7, .) %>%
    gsub("Unclassified", 7, .) %>%
    gsub("H3K4me1", 8, .) 
  return(column)
}

np_fheat$chromhmm_b_primed <- chromhmm_state2_num(np_fheat$chromhmm_b_primed)
np_fheat$chromhmm_b_naive <- chromhmm_state2_num(np_fheat$chromhmm_b_naive)
np_fheat$chromhmm_oe_primed <- chromhmm_state2_num(np_fheat$chromhmm_oe_primed)
np_fheat$chromhmm_oe_naive <- chromhmm_state2_num(np_fheat$chromhmm_oe_naive)



#order on primed dist ----------------------------------------------------------


p_fheat_1k_long <- np_fheat[order(np_fheat$distance_primed, decreasing = TRUE),][1:1000,]

long <- as.data.frame(cbind(log2(p_fheat_1k_long[c("distance_primed","distance_naive")]), 
                            p_fheat_1k_long[c("chromhmm_b_primed", "chromhmm_oe_primed",
                                              "chromhmm_b_naive", "chromhmm_oe_naive", 
                                              "b_ID", "oe_ID")]))




n_fheat_1k_long <- np_fheat[order(np_fheat$distance_naive, decreasing = TRUE),][1:1000,]

n_long <- as.data.frame(cbind(log2(n_fheat_1k_long[c("distance_primed","distance_naive")]), 
                              n_fheat_1k_long[c("chromhmm_b_primed", "chromhmm_oe_primed",
                                                "chromhmm_b_naive", "chromhmm_oe_naive", 
                                                "b_ID", "oe_ID")]))

#Colours have been changed for the final paper figure!

# active 1 - #00AE0F
# background 2 - #c9c9c9
# bivalent 3: #FF8300
# heterochromatin 4: #C300FF
# Mixed 5- #6e6e6e
# polycomb 6: #FF0100
# unknown 7 - #eedd9a
# h3k4me1 8 - #297800


set.seed(3)
km <- kmeans(as.matrix(long[,3:4]), centers = 6)$cluster

hc1 <- hclust(dist(long[,3:4]))
od1 <- hc1$order


#         1         2           3         4           5          6            8
col <- c('#00AE0F', '#c9c9c9', '#FF8300', '#C300FF', '#6e6e6e', '#FF0100', '#297800')

h1 <- Heatmap(as.matrix(long[,3:4]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, 
              split=long$chromhmm_b_primed, col=col, km=km, row_order = od1) +
  Heatmap(as.matrix(long[,1]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, row_order = od1,
          circlize::colorRamp2(c(-1, 0, 30), c("grey", "blue", "red")))

#                      1         2           3         4           5          6         7          8
# col_long_naive <- c('#00AE0F', '#c9c9c9', '#FF8300', '#C300FF', '#6e6e6e', '#FF0100', '#eedd9a', '#297800')

h2 <- Heatmap(as.matrix(long[,5:6]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, 
              split=long$chromhmm_b_primed, col=col, km=km, row_order = od1) +
  Heatmap(as.matrix(long[,2]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, row_order = od1,
          circlize::colorRamp2(c(-1, 0, 30), c("grey", "blue", "red")))


pdf("5_network_distance/chrommhmm_hist_all_heatmap_primed_20200912.pdf")
h1+h2
dev.off()

# Get numbers for long interaction heatmap -------------------------------------
# table(is.na(long$distance_primed), is.na(long$distance_naive))
# 
# table(long$chromhmm_b_primed, long$chromhmm_b_naive)
# table(long$chromhmm_b_primed)
# table(long$chromhmm_b_naive)

#order on naive dist -----------------------------------------------------------


n_fheat_1k_long <- np_fheat[order(np_fheat$distance_naive, decreasing = TRUE),][1:1000,]

n_long <- as.data.frame(cbind(log2(n_fheat_1k_long[c("distance_primed","distance_naive")]), 
                              n_fheat_1k_long[c("chromhmm_b_primed", "chromhmm_oe_primed",
                                                "chromhmm_b_naive", "chromhmm_oe_naive", 
                                                "b_ID", "oe_ID")]))


#Colours have been changed for the final paper figure!

# active 1 - #00AE0F
# background 2 - #c9c9c9
# bivalent 3: #FF8300
# heterochromatin 4: #C300FF
# Mixed 5- #6e6e6e
# polycomb 6: #FF0100
# unknown 7 - #eedd9a
# h3k4me1 8 - #297800
set.seed(3)
nkm <- kmeans(as.matrix(n_long[,5:6]), centers = 6)$cluster

nhc2 <- hclust(dist(n_long[,5:6]))
nod2 <- nhc2$order


#           1         2          3           4          5          6         8 
ncol <-  c('#00AE0F', '#c9c9c9', '#FF8300', '#C300FF', "#6e6e6e", '#FF0100', '#297800')
#           1           2         3           4          5           6          8
ncol2 <-  c('#00AE0F', '#c9c9c9', '#FF8300', '#C300FF', "#6e6e6e", '#FF0100', '#297800')


nh1 <- Heatmap(as.matrix(n_long[,3:4]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, 
               split=n_long$chromhmm_b_naive, col=ncol, km=nkm, row_order = nod2) +
  Heatmap(as.matrix(n_long[,1]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, row_order = nod2,
          circlize::colorRamp2(c(0, 30), c("blue", "red")))

nh2 <- Heatmap(as.matrix(n_long[,5:6]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, 
               split=n_long$chromhmm_b_naive, col=ncol2, km=nkm, row_order = nod2) +
  Heatmap(as.matrix(n_long[,2]), show_row_names = FALSE, gap = unit(2, "mm"), cluster_columns = FALSE, row_order = nod2,
          circlize::colorRamp2(c(0, 30), c("blue", "red")))


pdf("5_network_distance/chrommhmm_hist_all_heatmap_naive_20200912.pdf")
nh1+nh2
dev.off()




#ChromHMM states by distance heatmap -------------------------------------------

links_nodes_cat_col_coord_deb2b <- readRDS("2_network_make/links_nodes_cat_col_coord_deb2b_20200911.rds")


allu_df <- data.frame(b_ID=as.character(links_nodes_cat_col_coord_deb2b$links$b_ID), 
                      oe_ID=as.character(links_nodes_cat_col_coord_deb2b$links$oe_ID),
                      origin=links_nodes_cat_col_coord_deb2b$links$origin, 
                      distance=log2(links_nodes_cat_col_coord_deb2b$links$distance))

#add chromhmm state from nodes_all df
naive_chromhmm_lkup <- data.frame(b_ID=naive_nodes$ID, chromhmm_naive_b=naive_nodes$chromhmm)
allu_df <- merge(allu_df, naive_chromhmm_lkup, by="b_ID", all.x=TRUE)
data.table::setnames(naive_chromhmm_lkup, c("b_ID", "chromhmm_naive_b"), c("oe_ID", "chromhmm_naive_oe"))
allu_df <- merge(allu_df, naive_chromhmm_lkup, by="oe_ID", all.x=TRUE)

primed_chromhmm_lkup <- data.frame(b_ID=primed_nodes$ID, chromhmm_primed_b=primed_nodes$chromhmm)
allu_df <- merge(allu_df, primed_chromhmm_lkup, by="b_ID", all.x=TRUE)
data.table::setnames(primed_chromhmm_lkup, c("b_ID", "chromhmm_primed_b"), c("oe_ID", "chromhmm_primed_oe"))
allu_df <- merge(allu_df, primed_chromhmm_lkup, by="oe_ID", all.x=TRUE)

allu_df$chromhmm_naive <- paste(allu_df$chromhmm_naive_b, allu_df$chromhmm_naive_oe, sep="_")
allu_df$chromhmm_primed <- paste(allu_df$chromhmm_primed_b, allu_df$chromhmm_primed_oe, sep="_")

allu_df$chromhmm <- ifelse(allu_df$origin %in% "naive", allu_df$chromhmm_naive, allu_df$chromhmm_primed)
allu_df$chromhmm_origin <- paste(allu_df$chromhmm, allu_df$origin, sep="_")


#' remove outliers in data using Turky method
#' 
#' @param dataf data.frame of value with chromhmm_origin
#' @param var the chromhmm_origin variable for which to remove outliers
#' 
#' 
outlier_dist <- function(dataf, var) {
  var_data <- dataf$distance[dataf$chromhmm_origin %in% var]
  
  #get outliers using Turkey method (1.5*IQR)
  outlier <- boxplot.stats(var_data)$out
  
  # mo <- mean(outlier)
  
  var_na_outlier <- dataf[dataf$chromhmm_origin %in% var,][!(var_data %in% outlier),]
  
  return(var_na_outlier)
}



out_df <- list()
for(i in names(table(allu_df$chromhmm_origin))){
  out_df[[i]] <- outlier_dist(allu_df, i)
}

allu_outlier <- plyr::rbind.fill(out_df)


#split into cell type
allu_outlier_n <- allu_outlier[allu_outlier$origin %in% "naive",]
allu_outlier_p <- allu_outlier[allu_outlier$origin %in% "primed",]

#distance intervals
allu_outlier_p$interval <- cut(allu_outlier_p$distance, breaks=(seq(10,28,1.5)), include.lowest = TRUE)
#list of tables of for each interval and chromhmm state
allu_out_int_p <- tapply(allu_outlier_p$interval, allu_outlier_p$chromhmm, function(x) as.data.frame(table(x)))

#rename freq col by file name
for(i in seq_along(allu_out_int_p)){
  colnames(allu_out_int_p[[i]]) <- c("interval", names(allu_out_int_p)[i])
}
#reduce into a single data.frame
allu_out_int_p_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "interval", all = TRUE), allu_out_int_p)

#do the same for primed
allu_outlier_n$interval <- cut(allu_outlier_n$distance, breaks=(seq(10,28,1.5)), include.lowest = TRUE)
allu_out_int_n <- tapply(allu_outlier_n$interval, allu_outlier_n$chromhmm, function(x) as.data.frame(table(x)))

#rename freq col by file name
for(i in seq_along(allu_out_int_n)){
  colnames(allu_out_int_n[[i]]) <- c("interval", names(allu_out_int_n)[i])
}

allu_out_int_n_df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "interval", all = TRUE), allu_out_int_n)


#add in missing rows to make heatmap
missing_cols <- colnames(allu_out_int_n_df)[!(colnames(allu_out_int_n_df) %in% colnames(allu_out_int_p_df))]
for(i in 1:length(missing_cols)){
  allu_out_int_p_df[missing_cols[i]] <- rep(0, 12)
}


#order column (same order in both)
allu_out_int_p_df <- allu_out_int_p_df[colnames(allu_out_int_n_df)]


allu_out_int_p_df_log2 <- data.frame(interval=allu_out_int_p_df[1], log2(allu_out_int_p_df[2:ncol(allu_out_int_p_df)]+1))
allu_out_int_p_df_t <- t(allu_out_int_p_df_log2[2:ncol(allu_out_int_p_df_log2)])
colnames(allu_out_int_p_df_t) <- as.character(allu_out_int_p_df$interval)

# Heatmap(allu_int_p_df_t, cluster_columns = FALSE, col = rev(RColorBrewer::brewer.pal(11, "Spectral")), km=5)

allu_out_int_n_df_log2 <- data.frame(interval=allu_out_int_n_df[1], log2(allu_out_int_n_df[2:ncol(allu_out_int_n_df)]+1))
allu_out_int_n_df_t <- t(allu_out_int_n_df_log2[2:ncol(allu_out_int_n_df_log2)])
colnames(allu_out_int_n_df_t) <- as.character(allu_out_int_n_df$interval)


row_ord_p <- hclust(dist(allu_out_int_p_df_t))
row_ord_n <- hclust(dist(allu_out_int_n_df_t))


allu_out_int_p_df_t[allu_out_int_p_df_t == 0] <- NA
allu_out_int_n_df_t[allu_out_int_n_df_t == 0] <- NA

p_n_mat <- ifelse(is.na(allu_out_int_p_df_t), 0, allu_out_int_p_df_t) - ifelse(is.na(allu_out_int_n_df_t), 0, allu_out_int_n_df_t)
set.seed(2)
p_n_mat_km <- kmeans(p_n_mat, centers = 5)$cluster



pdf("5_network_distance/chromhmm_heatmap_without_outliers_p-n_20200912.pdf", width=5.5, height=7)
set.seed(2)
Heatmap(p_n_mat, cluster_columns = FALSE, cluster_rows=FALSE, col = rev(c("#ef8a62", "#f7f7f7", "#67a9cf")), 
        km=5, row_order=row_ord_p$order, gap = unit(2, "mm"))
grid.newpage()
#NA as grey
p_n_mat[ifelse(is.na(allu_out_int_p_df_t) & is.na(allu_out_int_n_df_t), TRUE, FALSE)] <- NA
Heatmap(p_n_mat, cluster_columns = FALSE, cluster_rows=FALSE, col = rev(c("#ef8a62", "#f7f7f7", "#67a9cf")), 
        split=paste0("km", p_n_mat_km), row_order=row_ord_p$order, gap = unit(2, "mm"))
dev.off()




#make boxplot of distances------------------------------------------------------

plot_dist <- log10(links_nodes_cat_col_coord_deb2b$links$distance)

naive_dist <- data.frame(distance=plot_dist[links_nodes_cat_col_coord_deb2b$links$origin == "naive"])
naive_dist$range <- "other"
#remove trans
naive_dist <- naive_dist[!is.na(naive_dist$distance),]

naive_dist$range[order(naive_dist$distance, decreasing = TRUE)][1:1000] <- "long"
naive_dist$range[order(naive_dist$distance, decreasing = TRUE)][(round(length(naive_dist$distance)/2)-499):(round(length(naive_dist$distance)/2)+500)] <- "mid"
naive_dist$range[order(naive_dist$distance, decreasing = FALSE)][1:1000] <- "short"

naive_dist$origin <- "naive"

#remove other
naive_dist <- naive_dist[naive_dist$range != "other",]


primed_dist <- data.frame(distance=plot_dist[links_nodes_cat_col_coord_deb2b$links$origin == "primed"])
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
  # geom_sina(aes(y = distance, x = origin, color=origin),data=boxplot_dist, size=0.5, alpha=0.2)+
  geom_boxplot(aes(y = distance, x = origin),data=boxplot_dist, outlier.size=0, size = 0.2) +
  scale_y_continuous(breaks = c(4,5,6,7,8),labels = c('10 kb', '100 kb', '1 Mb', '10 Mb', '100 Mb')) +
  ylab("log10 interaction distance")
box_p + facet_grid(~range) + theme_bw()
  
ggsave("5_network_distance/distance_boxplot_log10.pdf", device="pdf",
       width = 9, height = 9, units="cm")




  