
#Make S4 class for interactions
#TODO: write the documentation for this http://r-pkgs.had.co.nz/man.html

setClass(Class = "Interactions",
         slots = c(
           node = "character",
           genes = "character",
           naive_enhancers = "data.table",
           primed_enhancers = "data.table",
           naive_states = "table",
           primed_states = "table",
           interactions = "data.table"),
         prototype = list(
           node = NA_character_,
           genes = NA_character_,
           naive_enhancers = data.table(NULL),
           primed_enhancers = data.table(NULL),
           naive_states = table(NULL),
           primed_states = table(NULL),
           interactions = data.table(NULL)
         )
)



setGeneric("node", function(object) standardGeneric("node"))
setMethod("node", "Interactions", function(object) object@node)

setGeneric("genes", function(object) standardGeneric("genes"))
setMethod("genes", "Interactions", function(object) object@genes)

setGeneric("genes<-", function(object, value) standardGeneric("genes<-"))
setReplaceMethod("genes", signature(object="Interactions", value="character"),
                  function(object, value) {
                    object@genes <- value
                    return(object)
                  })


# setClassUnion("missingOrNULL", c("missing", "NULL"))
setGeneric("interactions", function(object, type=NULL) standardGeneric("interactions"), 
           signature = c("object"))

setMethod("interactions", "Interactions",
          function(object, type=NULL) {
            #if type not specified return full table
            if(is.null(type)){
              return(object@interactions)
            }else if(type == "naive"){
              out <- object@interactions[object@interactions$new_origin %in% 
                                           c("naive", "naive_primed"),] 
              return(out)
            }else if(type == "primed"){
              out <- object@interactions[object@interactions$new_origin %in% 
                                           c("primed", "naive_primed"),] 
              return(out)
            }else{
              stop("Type can be either 'naive' or 'primed'")
            }
            
            
          })


setGeneric("primed_interactions", function(object) standardGeneric("primed_interactions"))
setMethod("primed_interactions", "Interactions", function(object) {
  object[object@new_origin %in% c("primed", "naive_primed"),]
})


setGeneric("enhancers", function(object, type) standardGeneric("enhancers"))
setMethod("enhancers", signature(object="Interactions", type="character"),
          function(object, type) {
            if(type == "naive"){
              return(object@naive_enhancers)
            }else if(type == "primed"){
              return(object@primed_enhancers)
            }else{
              stop("Type needs to be either 'naive' or 'primed'")
            }
          })


setGeneric("enhsummary", function(object, type, gene_dt=NULL) standardGeneric("enhsummary"),
           signature = c("object", "type"))
setMethod("enhsummary", signature(object="Interactions", type="character"),
          function(object, type, gene_dt=NULL) {
            node_id <- paste(object@node, collapse = ";") #multiple baits for same gene
            #if not gene lookup dt provided
            if(is.null(gene_dt)){
              print("Gene won't be added to counted enhancers, node used instead")
              
              #aggregate by node
              enh_dt <- attr(object, paste(type, "enhancers", sep = "_"))
              #if empty table
              if(nrow(enh_dt) == 0){
                out <- rbind(enh_dt, list(int_node = node_id), fill = TRUE)
              }else{
                enh_dt$int_node <- node_id
                #produces data.frame
                out <- aggregate(enh_dt[,!"int_node", with=FALSE], 
                                 by = list(int_node=enh_dt$int_node), 
                                 FUN = function(x) paste(x, collapse = ";"))
                setDT(out)
              }
              
            #check the lookup table has the correct columns
            #can have multiple hindiii for a single gene - include count of fragments
            }else if(all(names(gene_dt) %in% c("ID", "prot_genes"))){
              #make sure the key is set to ID
              if(key(gene_dt) != "ID"){
                setkey(gene_dt, key="ID")
              }
              #aggregate by gene
              enh_dt <- attr(object, paste(type, "enhancers", sep = "_"))
              
              found_genes <- gene_dt[object@node]
              node_genes <- paste(unique(found_genes$prot_genes), collapse = ";")
              #if empty table
              if(nrow(enh_dt) == 0){
                out <- rbind(enh_dt, list(prot_genes = node_genes), 
                             fill = TRUE)
              }else{
                
                enh_dt$prot_genes <- node_genes
  
                out <- aggregate(enh_dt[,!"prot_genes", with=FALSE],
                                 by = list(prot_genes=enh_dt$prot_genes),
                                 FUN = function(x) paste(x, collapse = ";"))
                setDT(out)
              }
            #if lookup table has an incorrect format
            }else{
              stop(cat("
                  Gene lookup data.table needs to be in the following format:\n
                  ID      prot_genes
                  chr1_1  genex
                  chr1_2  geney
                  ...   ...
                  "))
            }
            
            out[, paste("hindiii_count", type, sep = "_") := list(length(object@node))]
            setnames(out, old = "ID", new = paste("oe", type, sep = "_"))
            return(out)
          })


setGeneric("statecounts", function(object, type, gene_dt=NULL) standardGeneric("statecounts"),
           signature = c("object", "type"))
setMethod("statecounts", signature(object="Interactions", type="character"),
          function(object, type, gene_dt=NULL) {
            
            state_counts <- attr(object, paste(type, "states", sep = "_"))
            count_col <- paste(state_counts, collapse = ";")
            
            #if not gene lookup dt provided
            if(is.null(gene_dt)){
              print("Gene won't be added to state counts, node used instead")
              
              # out <- data.frame(rbind(state_counts), row.names = NULL)
              
              out <- data.frame(int_node=paste(object@node, collapse = ";"), 
                                tmp=count_col)
              names(out)[names(out) == "tmp"] <- paste("APBH1MBU", type, sep = "_")
             
              
              #check the lookup table has the correct columns
              #can have multiple hindiii for a single gene - include count of fragments
            }else if(all(names(gene_dt) %in% c("ID", "prot_genes"))){
              #make sure the key is set to ID
              if(key(gene_dt) != "ID"){
                setkey(gene_dt, key="ID")
              }
              
              found_genes <- gene_dt[object@node]
              
              node_genes <- paste(unique(found_genes$prot_genes), collapse = ";")
              
              out <- data.frame(prot_genes=node_genes, tmp=count_col)
              names(out)[names(out) == "tmp"] <- paste("APBH1MBU", type, sep = "_")
              # out$prot_genes <- gene_dt[object@node]$prot_genes
              
              #if lookup table has an incorrect format
            }else{
              stop(cat("
                  Gene lookup data.table needs to be in the following format:\n
                  ID      prot_genes
                  chr1_1  genex
                  chr1_2  geney
                  ...   ...
                  "))
            }
            return(out)
          })


setGeneric("enhcounts", function(object, type, gene_dt=NULL) standardGeneric("enhcounts"),
           signature = c("object", "type"))
setMethod("enhcounts", signature(object="Interactions", type="character"),
          function(object, type, gene_dt=NULL) {
            node_id <- paste(object@node, collapse = ";") #multiple baits for same gene
            
            #get enhancer table
            enh_dt <- attr(object, paste(type, "enhancers", sep = "_"))
            
            if(nrow(enh_dt) != 0){
              #count enhancer, count superenhancers, count OSN
              #will produce warning as converting character "NA" to integer
              enhs <- unlist(strsplit(enh_dt[[paste("type", type, sep = "_")]], ";"))
              enh_count <- table(factor(enhs, levels = c("E", "SE")))
              osns <- unlist(strsplit(enh_dt[[paste("OSN", type, sep = "_")]], ";"))
              osn_count <- sum(table(osns), na.rm = TRUE) #numbers - 1 naive 2 shared 3 primed
            }
             
            
            #if not gene lookup dt provided
            if(is.null(gene_dt)){
              print("Gene won't be added to counted enhancers, node used instead")
              
              #if empty table
              if(nrow(enh_dt) == 0){
                out <- data.table(int_node=paste(object@node, collapse = ";"),
                                  SE=0,
                                  E=0, 
                                  OSN_count=0)
              }else{
                out <- data.table(int_node=paste(object@node, collapse = ";"),
                                  SE=unname(enh_count["SE"]),
                                  E=unname(enh_count["E"]), 
                                  OSN_count=osn_count)
              }
              
              #check the lookup table has the correct columns
              #can have multiple hindiii for a single gene - include count of fragments
            }else if(all(names(gene_dt) %in% c("ID", "prot_genes"))){
              #make sure the key is set to ID
              if(key(gene_dt) != "ID"){
                setkey(gene_dt, key="ID")
              }
             
              found_genes <- gene_dt[object@node]
              node_genes <- paste(unique(found_genes$prot_genes), collapse = ";")
              #if empty table
              if(nrow(enh_dt) == 0){
                out <- data.table(prot_genes=node_genes,
                                  SE=0,
                                  E=0, 
                                  OSN_count=0)
              }else{
                out <- data.table(prot_genes=node_genes,
                                  SE=unname(enh_count["SE"]),
                                  E=unname(enh_count["E"]), 
                                  OSN_count=osn_count)
              }
              #if lookup table has an incorrect format
            }else{
              stop(cat("
                  Gene lookup data.table needs to be in the following format:\n
                  ID      prot_genes
                  chr1_1  genex
                  chr1_2  geney
                  ...   ...
                  "))
            }
            #rename columns with type in name
            setnames(out, old = "SE", new = paste("SE", type, sep = "_"))
            setnames(out, old = "E", new = paste("E", type, sep = "_"))
            setnames(out, old = "OSN_count", new = paste("OSN_count", type, sep = "_"))
            return(out)
          })




