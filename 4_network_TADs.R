library(readr)
library(rtracklayer)
library(regioneR)
library(data.table)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(grid)
library(igraph)
library(dplyr)
library(dabestr)
library(ggpmisc)


source("network_tads_functions.R")
source("network_dynamics_functions.R")

#requires "BSgenome.Hsapiens.UCSC.hg38" to be installed

#needed to edit plot.permTestResultsList (submitted pull request, newest version probably fine)
source("plot.permTestResults.R")
source("plot.permTestResultsList.R")

#TAD and community correlation -------------------------------------------------

naive_tads <- import("4_network_TADs/naive_hESC-1_exact_TADs_5000_25000.washu_sorted.bed")
primed_tads <- import("4_network_TADs/primed_hESC-1_exact_TADs_5000_25000.washu_sorted.bed")


net_communities_dt <- readRDS("3_network_dynamics/net_communities_dt.rds")
hindiii_lookup_dt <- readRDS("4_network_TADs/hindiii_coord.rds")

#Make a gr with community data
com_gr_list <- apply(net_communities_dt, 1, function(x){
  com_gr <- nodes2granges(x["nodes"], hindiii_lookup_dt)
  com_gr$com <- x[".id"]
  return(com_gr)
})

#Remove alt chromosomes
com_gr_all <- do.call("c", com_gr_list)
seqlevels(com_gr_all) <- setdiff(seqlevels(com_gr_all), c("chrY", "chrMT"))

hg38_ref <- BSgenome.Hsapiens.UCSC.hg38
sequences_to_keep <- paste0("chr", c(1:22, "X"))
hg38 <- keepBSgenomeSequences(hg38_ref, sequences_to_keep)


set.seed(3)
naive_tad_pt <- permTest(A=naive_tads, B=com_gr_all, ntimes=500, randomize.function=randomizeRegions, 
                         genome="hg38", allow.overlaps=TRUE, per.chromosome=TRUE,
                         evaluate.function=tco_median)

primed_tad_pt <- permTest(A=primed_tads, B=com_gr_all, ntimes=500, randomize.function=randomizeRegions, 
                          genome="hg38", allow.overlaps=TRUE, per.chromosome=TRUE,
                          evaluate.function=tco_median)


#make plots
pdf("4_network_TADs/primed_tad_permtest.pdf", width = 5, height = 5)
plot(primed_tad_pt, xlim=c(55, 85))
grid.newpage()
plot(naive_tad_pt, xlim=c(55, 85))
dev.off()

#Single permutation plot -------------------------------------------------------

naive_percentage <- tads_community_overlap(net_communities_dt$nodes, naive_tads, 
                                           hindiii_lookup_dt, percentage=TRUE)
primed_percentage <- tads_community_overlap(net_communities_dt$nodes, primed_tads, 
                                            hindiii_lookup_dt, percentage=TRUE)

#need to randomly shuffle tad coordinates
set.seed(3)
mc.set.seed <- FALSE
naive_tads_random <- regioneR::randomizeRegions(naive_tads, genome=hg38, allow.overlaps=TRUE, per.chromosome=TRUE)
primed_tads_random <- regioneR::randomizeRegions(primed_tads, genome=hg38, allow.overlaps=TRUE, per.chromosome=TRUE)

naive_percentage_random <- tads_community_overlap(net_communities_dt$nodes, naive_tads_random, hindiii_lookup_dt, TRUE)
primed_percentage_random <- tads_community_overlap(net_communities_dt$nodes, primed_tads_random, hindiii_lookup_dt, TRUE)


tad_community_bxplt <- rbind(data.frame(origin=rep("naive", length(naive_percentage)), percentage=naive_percentage),
                             data.frame(origin=rep("naive_random", length(naive_percentage_random)), percentage=naive_percentage_random),
                             data.frame(origin=rep("primed", length(primed_percentage)), percentage=primed_percentage),
                             data.frame(origin=rep("primed_random", length(primed_percentage_random)), percentage=primed_percentage_random))

tad_community_bxplt$percentage[is.na(tad_community_bxplt$percentage)] <- 0

my_comparisons_tad <- list(c("naive", "naive_random"), c("primed", "primed_random"))


ggviolin(tad_community_bxplt, x = "origin", y = "percentage",
         fill = "origin", trim=TRUE,
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons_tad, method= "wilcox.test") +
  scale_y_continuous(breaks=c(0,25, 50, 75,100))
ggsave("4_network_TADs/tad_overlap_community_hindiii_wilcox.pdf", 
       width = 5, height = 5, device = "pdf")



# Insulation scores of community assigned TAD(s) with centrality score ---------



naive_tads_boundaries <- import("4_network_TADs/naive_hESC-1_25000rounded_TADBs_5000_25000.washu.bed")
primed_tads_boundaries <- import("4_network_TADs/primed_hESC-1_25000rounded_TADBs_5000_25000.washu.bed")

# for each hindiii, assign a tad 
com_gr_all

net_all_s <- readRDS("2_network_make/net_all_s_20200911.rds")


in_score <- data.frame(ID=V(net_all_s)$name,
                   in_primed=V(net_all_s)$in_primed,
                   in_naive=V(net_all_s)$in_naive,
                   tad_primed=V(net_all_s)$tad_primed_col,
                   tad_naive=V(net_all_s)$tad_naive_col)

# Add insulation scores
mcols(com_gr_all) <- merge(data.frame(mcols(com_gr_all)), in_score, by="ID", 
                           all.x=TRUE, sort=FALSE)

# get TAD boundary HindIII fragments and plot their insulation score
nb <- findOverlaps(com_gr_all, naive_tads_boundaries)
pb <- findOverlaps(com_gr_all, primed_tads_boundaries)

mcols(com_gr_all)$n_boundary_in <- NA
mcols(com_gr_all)$p_boundary_in <- NA

com_gr_all$n_boundary_in[queryHits(nb)] <- com_gr_all$in_naive[queryHits(nb)]
com_gr_all$p_boundary_in[queryHits(pb)] <- com_gr_all$in_primed[queryHits(pb)]

boundaries <- com_gr_all[!is.na(com_gr_all$n_boundary_in) & !is.na(com_gr_all$p_boundary_in)]

community_in_mean <- as_tibble(mcols(boundaries)) %>% group_by(com) %>% summarise(n_in_mean = mean(n_boundary_in), p_in_mean = mean(p_boundary_in))
names(community_in_mean)[1] <- ".id"

centrality_insulation <- merge(community_in_mean, net_communities_dt, by=".id", all.x=TRUE)
# boxplot(centrality_insulation$p_in_mean,centrality_insulation$n_in_mean, 
#         log2(centrality_insulation$p_btwn_cent_cnt),log2(centrality_insulation$n_btwn_cent_cnt))


i <- centrality_insulation[c("n_in_mean", "p_in_mean")]
i$ID <- seq(1,nrow(i))
i$color <- ifelse(i$n_in_mean < i$p_in_mean, "up", "down")

i <- reshape2::melt(i, id.var=c("ID", "color"))

# Estimation plots
i_group_paired <- 
  i %>%
  dabest(variable, value, 
         idx = c("n_in_mean", "p_in_mean"), 
         paired = TRUE, id.col = ID)

i_group_paired %>% 
  mean_diff() %>% 
  plot(color.column = color, slopegraph.params = list(alpha=0.5))
ggsave("4_network_TADs/insulation_tadb_subnetworks_paired.pdf")

# Estimation plots
i_group_paired <- 
  i %>%
  dabest(variable, value, 
         idx = c("n_in_mean", "p_in_mean"), 
         paired = FALSE, id.col = ID)

i_group_paired %>% 
  mean_diff() %>% 
  plot(color.column = color, slopegraph.params = list(alpha=0.5))
ggsave("4_network_TADs/insulation_tadb_subnetworks.pdf")


c <- centrality_insulation[c("n_btwn_cent_cnt", "p_btwn_cent_cnt")]
c$p_btwn_cent_cnt <- log2(c$p_btwn_cent_cnt)
c$n_btwn_cent_cnt <- log2(c$n_btwn_cent_cnt)
c$ID <- seq(1,nrow(c))
c$color <- ifelse(c$n_btwn_cent_cnt < c$p_btwn_cent_cnt, "up", "down")

c <- reshape2::melt(c, id.var=c("ID", "color"))

c_group_paired <- 
  c %>%
  dabest(variable, value, 
         idx = c("n_btwn_cent_cnt", "p_btwn_cent_cnt"), 
         paired = TRUE, id.col = ID)

c_group_paired %>% 
  mean_diff() %>% 
  plot(color.column = color, slopegraph.params = list(alpha=0.5))


c_group_unpaired <- 
  c %>%
  dabest(variable, value, 
         idx = c("n_btwn_cent_cnt", "p_btwn_cent_cnt"), 
         paired = FALSE, id.col = ID)

c_group_unpaired %>% 
  mean_diff() %>% 
  plot(color.column = color)
ggsave("4_network_TADs/centrality_tadb_subnetworks.pdf")

all_com <- reshape2::melt(net_communities_dt[,c(".id","p_btwn_cent_cnt", "n_btwn_cent_cnt")], id.var=c(".id"))
all_com$value <- log2(all_com$value)

all_c_group_paired <- 
  all_com %>%
  dabest(variable, value, 
         idx = c("n_btwn_cent_cnt", "p_btwn_cent_cnt"), 
         paired = TRUE, id.col = .id)

all_c_group_paired %>% 
  mean_diff() %>% 
  plot(slopegraph.params = list(alpha=0.5))
ggsave("4_network_TADs/centrality_tadb_subnetworks_paired.pdf")





# 
# 
# p <- ggviolin(i, x = "variable", y = "value",
#               color = "variable", palette = "jco", 
#               line.color = "gray", line.size = 0.4, short.panel.labs = FALSE)
# # Use only p.format as label. Remove method name.
# p + stat_compare_means(label = "p.format", paired = TRUE)
# 
# p <- ggviolin(c, x = "variable", y = "value",
#               color = "variable", palette = "jco", 
#               line.color = "gray", line.size = 0.4, short.panel.labs = FALSE)
# # Use only p.format as label. Remove method name.
# p + stat_compare_means(label = "p.format", paired = TRUE)
# 
# 
# my_comparisons <- list( c("n_in_mean", "p_in_mean"), c("n_btwn_cent_cnt", "p_btwn_cent_cnt"))
# ggboxplot(ci, x = "variable", y = "value",
#           color = "variable", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons, paired=TRUE)+ 
#   stat_compare_means(label.y = 10) 
# 


p <- qplot(centrality_insulation$p_in_mean-centrality_insulation$n_in_mean, 
      centrality_insulation$p_btwn_cent_cnt-centrality_insulation$n_btwn_cent_cnt) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  theme_bw()
ggsave("4_network_TADs/centrality_vs_insulation.pdf", plot = p)
 

# Modularity of each subnetwork ------------------------------------------------

naive_net <- isolate_cell_network(net_all_s, "naive")
primed_net <- isolate_cell_network(net_all_s, "primed")

# get modularity score of each community

keep_nodes <- unlist(strsplit(as.character(centrality_insulation$nodes[1]), ";"))

g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
g <- add_edges(g, c(1,6, 1,11, 6, 11))
wtc <- cluster_walktrap(g)
modularity(wtc)
modularity(g, membership(wtc))

#' Calculate modulatiry of a subnetwork
#' 
#' @param net whole network
#' @param nodes names of nodes to subset from net
#' 
subnet_modulatiry <- function(net, nodes){
  query_nodes <- match(nodes, V(net)$name)
  # Remove NA
  query_nodes <- query_nodes[!is.na(query_nodes)]
  sub_net <- induced.subgraph(graph=net,vids=query_nodes)
  
  # Calculate modularity
  mod <- modularity(sub_net, seq(1, length(V(sub_net))))
  
  return(mod)
}

all_scores <- list()
for(i in seq(nrow(centrality_insulation))){
  keep_nodes <- unlist(strsplit(as.character(centrality_insulation$nodes[i]), ";"))
  all_scores[[i]] <- data.frame(naive=subnet_modulatiry(naive_net, keep_nodes),
                                primed=subnet_modulatiry(primed_net, keep_nodes))

}

all_scores_df <- plyr::ldply(all_scores, cbind)


t_group_unpaired <- 
  test %>%
  dabest(variable, value, 
         idx = c("naive", "primed"), 
         paired = FALSE)

t_group_unpaired %>% 
  mean_diff() %>% 
  plot()
