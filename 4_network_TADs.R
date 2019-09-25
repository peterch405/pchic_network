library(readr)
library(rtracklayer)
library(regioneR)
library(data.table)
library(ggpubr)

source("network_tads_functions.R")
source("network_dynamics_functions.R")

#requires "BSgenome.Hsapiens.UCSC.hg38" to be installed

#needed to edit plot.permTestResultsList (submitted pull request)
source("plot.permTestResults.R")
source("plot.permTestResultsList.R")

#TAD and community correlation -------------------------------------------------

naive_tads <- import("4_network_TADs/naive_hESC-1_exact_TADs_5000_25000.washu_sorted.bed")
primed_tads <- import("4_network_TADs/primed_hESC-1_exact_TADs_5000_25000.washu_sorted.bed")


net_communities_dt <- readRDS("4_network_TADs/net_communities_dt.rds")
hindiii_lookup_dt <- readRDS("4_network_TADs/hindiii_coord.rds")

#Make a gr with community data
com_gr_list <- apply(net_communities_dt, 1, function(x){
  com_gr <- nodes2granges(x["nodes"], hindiii_lookup_dt)
  com_gr$com <- x[".id"]
  return(com_gr)
})

com_gr_all <- do.call("c", com_gr_list)


set.seed(3)
naive_tad_pt <- permTest(A=naive_tads, ntimes=500, randomize.function=randomizeRegions, 
                         genome="hg38", allow.overlaps=TRUE, per.chromosome=TRUE,
                         evaluate.function = tco_median, hindiii_lookup_dt=hindiii_lookup_dt, 
                         B=com_gr_all, verbose=TRUE, mc.cores=7, mc.set.seed=FALSE)

primed_tad_pt <- permTest(A=primed_tads, ntimes=500, randomize.function=randomizeRegions, 
                          genome="hg38", allow.overlaps=TRUE, per.chromosome=TRUE,
                          evaluate.function = tco_median, hindiii_lookup_dt=hindiii_lookup_dt, 
                          B=com_gr_all, verbose=TRUE, mc.cores=7, mc.set.seed=FALSE)

#make plots
pdf("4_network_TADs/primed_tad_permtest.pdf", width = 5, height = 5)
plot(primed_tad_pt, xlim=c(55, 85))
new.grid()
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
naive_tads_random <- regioneR::randomizeRegions(naive_tads, genome="hg38", allow.overlaps=TRUE, per.chromosome=TRUE)
primed_tads_random <- regioneR::randomizeRegions(primed_tads, genome="hg38", allow.overlaps=TRUE, per.chromosome=TRUE)

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