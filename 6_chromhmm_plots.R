

library(readr)
library(lattice)
library(igraph)
library(Gmisc)
library(pheatmap)
library(ggpubr)
library(GenomicRanges)

source("network_chromhmm_functions.R")

#Plot emmission state heatmap (didn't output on cluster) -----------------------

plot_emmisions("6_chromhmm_plots/emissions_16.txt",
               out_path="6_chromhmm_plots/emmisions_16.pdf")



#Make transition plots ---------------------------------------------------------

links_nodes_cat_col_coord_deb2b <- readRDS("2_network_make/links_nodes_cat_col_coord_deb2b_20190829.rds")
nodes_all <- links_nodes_cat_col_coord_deb2b$nodes
links_all <- links_nodes_cat_col_coord_deb2b$links

#all interacting nodes
interacting_nodes <- unique(nodes_all$ID)

naive_intera <- links_all[links_all$origin == "naive",]
primed_intera <- links_all[links_all$origin == "primed",]

naive_interacting_nodes <- unique(c(naive_intera$b_ID, naive_intera$oe_ID))
primed_interacting_nodes <- unique(c(primed_intera$b_ID, primed_intera$oe_ID))


sankey_act_nodes <- nodes_all[nodes_all$ID %in% interacting_nodes,] #same thing, not needed!

shared <- sankey_act_nodes$ID %in% naive_interacting_nodes & sankey_act_nodes$ID %in% primed_interacting_nodes
naive_int <- sankey_act_nodes$ID %in% naive_interacting_nodes & !(sankey_act_nodes$ID %in% primed_interacting_nodes)
primed_int <- !(sankey_act_nodes$ID %in% naive_interacting_nodes) & sankey_act_nodes$ID %in% primed_interacting_nodes



sankey_act_nodes$chromhmm_naive[naive_int] <- paste(sankey_act_nodes$chromhmm_naive[naive_int], "naive", sep = "_")
sankey_act_nodes$chromhmm_naive[primed_int] <- paste(sankey_act_nodes$chromhmm_naive[primed_int], "primed", sep = "_")
sankey_act_nodes$chromhmm_naive[shared] <- paste(sankey_act_nodes$chromhmm_naive[shared], "shared", sep = "_")

sankey_act_nodes$chromhmm_primed[primed_int] <- paste(sankey_act_nodes$chromhmm_primed[primed_int], "primed", sep = "_")
sankey_act_nodes$chromhmm_primed[naive_int] <- paste(sankey_act_nodes$chromhmm_primed[naive_int], "naive", sep = "_")
sankey_act_nodes$chromhmm_primed[shared] <- paste(sankey_act_nodes$chromhmm_primed[shared], "shared", sep = "_")

sankey_act_nodes_full <- nodes_all[nodes_all$ID %in% interacting_nodes,]

output_perc <- function(txt, n) sprintf("%s\n[%.0f%%]", txt, n)


sankey_act_nodes$chromhmm_primed <- gsub(" ", "_", sankey_act_nodes$chromhmm_primed)
sankey_act_nodes$chromhmm_naive <- gsub(" ", "_", sankey_act_nodes$chromhmm_naive)

#Not all used in publication ---------------------------------------------------

state_transition_plot(sankey_act_nodes[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_active_to_primed_interaction_20190114.pdf",
                      "of_subset", state2plot = c("Active_naive", "Active_primed", "Active_shared"), 
                      lvls = c("Active_naive", "Active_primed", "Active_shared",
                               "Background_naive", "Background_primed", "Background_shared",
                               "Bivalent_naive", "Bivalent_primed","Bivalent_shared",
                               "H3K4me1_naive", "H3K4me1_primed", "H3K4me1_shared",
                               "Heterochromatin_Repressed_naive", 
                               "Heterochromatin_Repressed_primed", "Heterochromatin_Repressed_shared", 
                               "Polycomb_Repressed_naive", "Polycomb_Repressed_primed", 
                               "Polycomb_Repressed_shared",
                               "Mixed_naive", "Mixed_primed", "Mixed_shared",
                               "Unclassified_naive", "Unclassified_primed", "Unclassified_shared"),
                      state_colours = list(Active_naive="#3aab04", Active_primed="#3aab04", Active_shared="#3aab04",
                                           Background_naive="#c2c2c2", Background_primed="#c2c2c2", Background_shared="#c2c2c2",
                                           Bivalent_naive="#ff9a01", Bivalent_primed="#ff9a01",Bivalent_shared="#ff9a01",
                                           H3K4me1_naive="#00800B", H3K4me1_primed="#00800B", H3K4me1_shared="#00800B",
                                           Heterochromatin_Repressed_naive="#b340d5", 
                                           Heterochromatin_Repressed_primed="#b340d5", Heterochromatin_Repressed_shared="#b340d5", 
                                           Polycomb_Repressed_naive="#d00101", Polycomb_Repressed_primed="#d00101", 
                                           Polycomb_Repressed_shared="#d00101",
                                           Mixed_naive="#eedd9a", Mixed_primed="#eedd9a", Mixed_shared="#eedd9a",
                                           Unclassified_naive="#6e6e6e", Unclassified_primed="#6e6e6e", Unclassified_shared="#6e6e6e"),
                      table_size = 1)

state_transition_plot(sankey_act_nodes[c("chromhmm_naive", "chromhmm_primed")],
                      percentage = "of_subset", state2plot = c("Active_naive", "Active_primed", "Active_shared"), 
                      lvls = c("Active_naive", "Active_primed", "Active_shared",
                               "Background_naive", "Background_primed", "Background_shared",
                               "Bivalent_naive", "Bivalent_primed","Bivalent_shared",
                               "H3K4me1_naive", "H3K4me1_primed", "H3K4me1_shared",
                               "Heterochromatin_Repressed_naive", 
                               "Heterochromatin_Repressed_primed", "Heterochromatin_Repressed_shared", 
                               "Polycomb_Repressed_naive", "Polycomb_Repressed_primed", 
                               "Polycomb_Repressed_shared",
                               "Mixed_naive", "Mixed_primed", "Mixed_shared",
                               "Unclassified_naive", "Unclassified_primed", "Unclassified_shared"),
                      state_colours = list(Active_naive="#3aab04", Active_primed="#3aab04", Active_shared="#3aab04",
                                           Background_naive="#c2c2c2", Background_primed="#c2c2c2", Background_shared="#c2c2c2",
                                           Bivalent_naive="#ff9a01", Bivalent_primed="#ff9a01",Bivalent_shared="#ff9a01",
                                           H3K4me1_naive="#00800B", H3K4me1_primed="#00800B", H3K4me1_shared="#00800B",
                                           `Heterochromatin Repressed_naive`="#b340d5", 
                                           `Heterochromatin Repressed_primed`="#b340d5", `Heterochromatin Repressed_shared`="#b340d5", 
                                           `Polycomb Repressed_naive`="#d00101", `Polycomb Repressed_primed`="#d00101", 
                                           `Polycomb Repressed_shared`="#d00101",
                                           Mixed_naive="#eedd9a", Mixed_primed="#eedd9a", Mixed_shared="#eedd9a",
                                           Unclassified_naive="#6e6e6e", Unclassified_primed="#6e6e6e", Unclassified_shared="#6e6e6e"))

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_all_to_primed_bi_poly_20190119.pdf",
                      "of_subset", state2plot=c("Bivalent","Polycomb Repressed"), col_filter="right")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      percentage = "of_subset", state2plot = "Active", col_filter = "left")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_active_to_primed_20190107.pdf",
                      "of_subset", "Active")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_background_to_primed_20190107.pdf",
                      "of_subset", "Background")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_Bivalent_to_primed_20190107.pdf",
                      "of_subset", "Bivalent")


state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_Heterochromatin_Repressed_to_primed_20190107.pdf",
                      "of_subset", "Heterochromatin Repressed")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_H3K4me1_to_primed_20190107.pdf",
                      "of_subset", "H3K4me1")


state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_Mixed_to_primed_20190107.pdf",
                      "of_subset", "Mixed")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_Unknown_to_primed_20190107.pdf",
                      "of_subset", "Unclassified")

state_transition_plot(sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")],
                      "6_chromhmm_plots/naive_Polycomb_Repressed_to_primed_20190107.pdf",
                      "of_subset", "Polycomb Repressed")


#Genome wide chromhmm numbers --------------------------------------------------

naive_nodes <- read_delim("0_network_data_preparation/naive_nodes_wchip_20181218.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

primed_nodes <- read_delim("0_network_data_preparation/primed_nodes_wchip_20181218.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

#remove mt
pdf("6_chromhmm_plots/chromhmm_states_heatmap_20190107.pdf", height = 6, width = 6)
pheatmap::pheatmap(log2(table(naive_nodes$chromhmm[!grepl("chrMT", naive_nodes$ID)], 
                              primed_nodes$chromhmm[!grepl("chrMT", primed_nodes$ID)])+1),
                   cluster_rows = F,
                   cluster_cols = F,
                   display_numbers = table(naive_nodes$chromhmm[!grepl("chrMT", naive_nodes$ID)], 
                                           primed_nodes$chromhmm[!grepl("chrMT", primed_nodes$ID)]))
dev.off()

#make heatmap of chromhmm fragments that interact

interacting_nodes_primed <- unique(c(links_all$b_ID[links_all$origin == "primed"], 
                                     links_all$oe_ID[links_all$origin == "primed"]))
interacting_nodes_naive <- unique(c(links_all$b_ID[links_all$origin == "naive"], 
                                    links_all$oe_ID[links_all$origin == "naive"]))

all_int_nodes <- unique(c(interacting_nodes_primed, interacting_nodes_naive))

naive_nodes_sub <- naive_nodes[naive_nodes$ID %in% all_int_nodes,]
primed_nodes_sub <- primed_nodes[primed_nodes$ID %in% all_int_nodes,]

#need to have an equal length table so need to NA fragments not present in cell type instead of removing
naive_nodes_sub$chromhmm[!(naive_nodes_sub$ID %in% interacting_nodes_naive)] <- NA
primed_nodes_sub$chromhmm[!(primed_nodes_sub$ID %in% interacting_nodes_primed)] <- NA


pdf("6_chromhmm_plots/chromhmm_states_heatmap_interacting_20190107.pdf", 
    height = 6, width = 6)
pheatmap::pheatmap(log2(table(naive_nodes_sub$chromhmm, 
                              primed_nodes_sub$chromhmm, useNA = "ifany")+1),
                   cluster_rows = F,
                   cluster_cols = F,
                   display_numbers = table(naive_nodes_sub$chromhmm, 
                                           primed_nodes_sub$chromhmm, useNA = "ifany"))
dev.off()



cnts_table <- table(naive_nodes_sub$chromhmm, primed_nodes_sub$chromhmm, useNA = "ifany")
naive_sum <- rowSums(cnts_table)#-diag(cnts_table)
primed_sum <- colSums(cnts_table)#-diag(cnts_table)

off_diag_cnts <- rbind(data.frame(count=unname(naive_sum), state=names(naive_sum), 
                                  type=rep("naive", length(naive_sum))), 
                       data.frame(count=unname(primed_sum), state=names(primed_sum), 
                                  type=rep("primed", length(primed_sum))))

ggplot(off_diag_cnts[!is.na(off_diag_cnts$state),], aes(x=state, y=count, fill=type))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,30000)) + 
  ylab("HindIII frequency") +
  xlab("ChromHMM states of interacting fragments")
ggsave("6_chromhmm_plots/hindiii_frequency_of_interacting_fragments_chromhmm_state_20190107.pdf", device="pdf",
       width = 8, height = 12, units="cm")



#number of interaction of each chromhmm state ----------------------------------
# test <- as_long_data_frame(as.undirected(net_all_s))

int_nodes_primed_cnt <- data.frame(table(c(links_all$b_ID[links_all$origin == "primed"], 
                                           links_all$oe_ID[links_all$origin == "primed"])))
names(int_nodes_primed_cnt) <- c("ID", "freq")

int_nodes_naive_cnt <- data.frame(table(c(links_all$b_ID[links_all$origin == "naive"], 
                                          links_all$oe_ID[links_all$origin == "naive"])))
names(int_nodes_naive_cnt) <- c("ID", "freq")


n_nodes_sub <- naive_nodes[c("ID", "chromhmm")][naive_nodes$ID %in% all_int_nodes,]
p_nodes_sub <- primed_nodes[c("ID", "chromhmm")][primed_nodes$ID %in% all_int_nodes,]

n_cnts <- merge(n_nodes_sub, int_nodes_naive_cnt, 
                by = "ID", all.x = TRUE)
p_cnts <- merge(p_nodes_sub, int_nodes_primed_cnt, 
                by = "ID", all.x = TRUE)

n_cnts_full <- n_cnts[!is.na(n_cnts$freq),]
np_cnts <- data.frame(table(rep(n_cnts_full$chromhmm, n_cnts_full$freq)))

p_cnts_full <- p_cnts[!is.na(p_cnts$freq),]

np_cnts <- merge(np_cnts, table(rep(p_cnts_full$chromhmm, p_cnts_full$freq)), 
                 by = "Var1", suffixes = c("_naive", "_primed"))

np_cnts_melt <- reshape2::melt(np_cnts)
np_cnts_melt$variable <- factor(np_cnts_melt$variable, levels =  c("Freq_primed", "Freq_naive"))


#This will double count bait to bait interactions
ggplot(np_cnts_melt, aes(x = Var1, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0)) + 
  ylab("Interaction frequency") +
  xlab("ChromHMM states of interacting fragments")
ggsave("6_chromhmm_plots/interaction_frequency_of_each_chromhmm_state_20190107.pdf", device="pdf",
       width = 8, height = 12, units="cm")


#Mixed states in naive cells ---------------------------------------------------
hindiii_coord <- readRDS("6_chromhmm_plots/hindiii_coord.rds")
hindiii_gr <- readRDS("6_chromhmm_plots/hindiii_gr.rds")


mixed_hindiii <- hindiii_gr[hindiii_gr$ID %in% sankey_act_nodes_full$ID[sankey_act_nodes_full$chromhmm_naive == "Mixed" & 
                                                                          sankey_act_nodes_full$chromhmm_primed %in% c("Bivalent", "Polycomb Repressed")]]

naive_16_seg_other <- rtracklayer::import("6_chromhmm_plots/naive_16_segments.bed",
                                          format = "bed")

naive_16_seg_chr <- diffloop::addchr(naive_16_seg_other)

mixed_seg <- findOverlaps(mixed_hindiii, naive_16_seg_chr)


collapse_states <- c('E1' = 'Active', 'E2'= 'Active',
                     'E3'= 'Active', 'E5'= 'Active', 'E6'= 'Active', 'E7'= 'Bivalent',
                     'E4'= 'Active',
                     'E8'= 'Polycomb Repressed', 'E9'= 'Bivalent',
                     'E10'= 'H3K4me1',
                     'E11'= 'H3K4me1',
                     'E12'= 'Unclassified', 'E13'= 'Heterochromatin Repressed',
                     'E14'= 'Background', 'E15'= 'Heterochromatin Repressed', 'E16'= 'Background')

collapse_states <- c('E1' = 'Active', 'E2'= 'Active',
                     'E3'= 'Active', 'E5'= 'Active', 'E6'= 'Active', 'E7'= 'Bivalent',
                     'E4'= 'Active',
                     'E8'= 'K27', 'E9'= 'Bivalent',
                     'E10'= 'Active',
                     'E11'= 'Active',
                     'E12'= 'Unclassified', 'E13'= 'Heterochromatin',
                     'E14'= '', 'E15'= 'Heterochromatin', 'E16'= '')

mcols(naive_16_seg_chr)$state <- sapply(naive_16_seg_chr$name, function(x) unname(collapse_states[x]))

mixed_naive_16_seg_chr <- naive_16_seg_chr[subjectHits(mixed_seg)]


mixed_plot_df <- data.frame(rep("Mixed", length(mixed_naive_16_seg_chr)), 
                            mixed_naive_16_seg_chr$state)


mixed_plot_df <- data.frame(hindiii=queryHits(mixed_seg), state=mixed_naive_16_seg_chr$state)
mixed_plot_df <- mixed_plot_df[!mixed_plot_df$state == "",]
mixed_agg <- aggregate(mixed_plot_df$state, by=list(mixed_plot_df$hindiii), function(x) paste(sort(unique(x)), collapse="_"))

mixed_plot_df <- data.frame(rep("Mixed", nrow(mixed_agg )), 
                            mixed_agg$x)



state_transition_plot(mixed_plot_df, lvls = c("Mixed",names(table(mixed_plot_df$mixed_agg.x))),
                      out_path = "6_chromhmm_plots/naive_mixed_decomposed_20190924.pdf",
                      state2plot = "Mixed", percentage = "of_total", col_filter = "other")








