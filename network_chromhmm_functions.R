

#' Plot the Chromhmm emisions heatmap
#' 
#' @param emissions_path path to emission.txt
#' @param out_path optional PDF output path
#' @param col_names column names of the input datasets
#' 
#' 
plot_emmisions <- function(emissions_path, out_path=NA, col_names=NA){
  
  emissions <- read_delim(emissions_path, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  #plot emmission heatmap
  m_h16a <- t(as.matrix(emissions))
  if(is.na(col_names[1])){
    slice_by <- c("state (Emission order)", "H3K9me3", "H3K27me3", "H3K4me3", 
                  "H3K27ac", "H3K4me1", "T2016_input", "T2014_input", 
                  "Ji_input", "Ac_input", "Ar_input")
  }else{
    slice_by <- col_names
  }
  
  m_h16a <- m_h16a[slice_by,,drop=FALSE]
  
  col_ramp <- colorRampPalette(c('white', 'blue'))(16)
  
  if(!is.na(out_path)){
    pdf(out_path, width = 5, height = 7)
  }
  
  # return(m_h16a)
  print(lattice::levelplot(m_h16a[2:nrow(m_h16a),1:ncol(m_h16a)], col.regions=col_ramp, xlab=NULL, 
                  ylab="State (emission order)",
                  scales=list(x=list(rot=90, cex=1), y=list(at=seq(0,ncol(m_h16a),1))),
                  main="Histone with 16 states", tick.number=ncol(m_h16a)))
  
  if(!is.na(out_path)){
    dev.off()
  }
  
}


#' Move an element in a vector to the first position
#' 
#' @param in_vect A vector for which element order is to be changed
#' @param move_which Which element inside the vector to move to front 
move_to_first <- function(in_vect, move_which){
  if(move_which %in% in_vect){
    rest <- setdiff(in_vect, move_which)
    
    out_vect <- c(move_which, rest)
  }else{
    stop("Element not present in vector provided")
  }
  return(out_vect)
}


#' Make a transitionPlot from Naive chromhmm state changing into Primed chromhmm state
#' 
#' @param two_col_data A data.frame with Naive states as first column and Primed states as second column
#' @param out_path Output pdf path
#' @param percentage of_total is the percentage of all fragments before subseting
#' of_subset is the percentage within the subset dataset
#' @param state2plot which states from the first column to plot
#' @return Saves a pdf to out_path
state_transition_plot <- function(two_col_data, out_path=NULL, percentage=c("of_total","of_subset"),
                                  state2plot=c("Active", "Background", "Bivalent", 
                                               "H3K4me1", "Heterochromatin Repressed", 
                                               "Polycomb Repressed", "Mixed", "Unclassified"),
                                  lvls=NULL, state_colours=NULL, col_filter=c("left", "right", "other"),
                                  table_size=7){
  
  #keep only naive active/selected state and see what they go to
  #by first column
  if(!col_filter %in% c("left", "right", "other")){
    stop("col_filter 'left' or 'right' or 'other'")
  }
  if(col_filter == "right"){
    tp_data <- sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")][sankey_act_nodes_full[c("chromhmm_naive", "chromhmm_primed")][,2] %in% state2plot,]
  }
  if(col_filter %in% c("left", "other")){
    tp_data <- two_col_data[two_col_data[,1] %in% state2plot,]
  }
  
  
  #rename unknown to mixed and mixed to unknown
  #the state2plot needs to come first
  if(is.null(lvls)){
    lvls <- c("Active", "Background", "Bivalent", "H3K4me1", "Heterochromatin Repressed", 
              "Polycomb Repressed", "Mixed", "Unclassified")
    #move state2plot first (moveMe a bit overkill)
    lvls <- move_to_first(lvls, state2plot)
  }
  
  print(lvls)
  # lvls <- moveMe(lvls, paste(state2plot, "first", sep = " "))
  
  tp_data[,1] <- factor(tp_data[,1], levels = lvls)
  tp_data[,2] <- factor(tp_data[,2], levels = lvls)
  
  if(is.null(state_colours)){
    state_colours <- list(Active="#3aab04", Background="#c2c2c2", Bivalent="#ff9a01", 
                          H3K4me1="#00800B", `Heterochromatin Repressed`="#b340d5", 
                          `Polycomb Repressed`="#d00101", Mixed="#eedd9a", Unclassified="#6e6e6e")
    
  }
  
  
  if(percentage == "of_total"){
    box_txt <- cbind(mapply(output_perc,
                            txt = lvls,
                            n = prop.table(table(tp_data[,1]))*100),
                     mapply(output_perc,
                            txt = lvls,
                            n = prop.table(table(tp_data[,2]))*100))
  }else{
    #make the % out of the total interacting fragments to give perspective
    box_txt <- cbind(mapply(output_perc, 
                            txt = lvls, 
                            n = prop.table(table(tp_data[,1]))*100),
                     mapply(output_perc, 
                            txt = lvls, 
                            n = (table(factor(tp_data[,2], levels = lvls))/
                                   table(factor(two_col_data[,2], levels = lvls)))*100))
  }
  
  # print(box_txt_naive)
  # print(table(tp_data[,1], tp_data[,2]))
  
  #need to remove the 0 column so the min_lwd will be dividing by 0
  count_tbl <- table(tp_data[,1], tp_data[,2])
  
  print(count_tbl)
  #which column is 0?
  if(col_filter %in% c("right", "other")){
    count_plot <- count_tbl+0.1
    remove_cols_rows <-NA
  }
  if(col_filter == "left"){
    remove_cols_rows <- names(which(colSums(count_tbl) == 0))
    count_plot <- count_tbl[!(lvls %in% remove_cols_rows),!(lvls %in% remove_cols_rows)]
  }
  
  # return(count_plot)
  
  
  print(remove_cols_rows)
  
  # grid::grid.newpage()
  if(!is.null(out_path)){
    pdf(out_path, height = 5, width = 5)
  }
  
  if(col_filter == "right"){
    Gmisc::transitionPlot(count_plot, #with(tp_data, table(naive, primed))
                          type_of_arrow = "simple",
                          box_label = c("Naive", "Primed"),
                          box_txt = box_txt[!(lvls %in% remove_cols_rows),], #also need to remove 0 count row
                          min_lwd = ggplot2::unit(0, "mm"),
                          max_lwd = ggplot2::unit(15, "mm"),
                          overlap_add_width = ggplot2::unit(1, "mm"), cex = 0.5,
                          fill_start_box = unname(unlist(state_colours[lvls[!(lvls %in% remove_cols_rows)]])), 
                          fill_end_box = unname(unlist(state_colours[state2plot])),
                          new_page = TRUE)
    if(!is.null(out_path)){
      #also save table with %
      grid::grid.newpage()
      gridExtra::grid.table(box_txt[!(lvls %in% remove_cols_rows),], theme=gridExtra::ttheme_default(base_size=table_size))
      #also add table with freq
      grid::grid.newpage()
      gridExtra::grid.table(count_plot, theme=gridExtra::ttheme_default(base_size=table_size))
    }
  }
  if(col_filter %in% c("left", "other")){
    Gmisc::transitionPlot(count_plot, #with(tp_data, table(naive, primed))
                          type_of_arrow = "simple",
                          box_label = c("Naive", "Primed"),
                          box_txt = box_txt[!(lvls %in% remove_cols_rows),], #also need to remove 0 count row
                          min_lwd = ggplot2::unit(0, "mm"),
                          max_lwd = ggplot2::unit(15, "mm"),
                          overlap_add_width = ggplot2::unit(1, "mm"), cex = 0.5,
                          fill_start_box = unname(unlist(state_colours[state2plot])), 
                          fill_end_box = unname(unlist(state_colours[lvls[!(lvls %in% remove_cols_rows)]])),
                          new_page = TRUE)
    
    if(!is.null(out_path)){
      #also save table with %
      grid::grid.newpage()
      gridExtra::grid.table(box_txt[!(lvls %in% remove_cols_rows),], theme=gridExtra::ttheme_default(base_size=table_size))
      #also add table with freq
      grid::grid.newpage()
      gridExtra::grid.table(data.frame(Count=count_plot[state2plot,]), theme=gridExtra::ttheme_default(base_size=table_size))
    }
  }
  if(!is.null(out_path)){
    dev.off()
  }
}
