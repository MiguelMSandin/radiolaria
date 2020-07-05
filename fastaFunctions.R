#---- 
#---- loading packages ----

library(ape)
library(phangorn)
library(seqinr)

library(dplyr)
library(tidyr)
library(data.table)

library(ggplot2)

#---- 
#---- Loading functions ----

# All files in fasta format (".fasta")

# SplitFasta: cuts a concatenated file vertically (e.g., 18S+28S) into two different files (e.g. 18S and 28S)
SplitFasta <- function(input, output_1, output_2, position, clean=F){
  
  if(grepl("\\.fasta$", input)){input<-gsub("\\.fasta$","",input)}
  
  file <- seqinr::read.fasta(paste0(input, ".fasta"), seqtype= "AA", as.string = T)
  file <- data.frame(name=paste(getAnnot(file)), sequence=paste(file))
  file$name <- gsub(">", "", file$name)
  
  file$seq1 <- substr(as.character(file$sequence), 1, as.numeric(position)-1)
  file$seq2 <- substr(as.character(file$sequence), as.numeric(position), nchar(as.character(file$sequence[1])))
  
  s1 <- lapply(file$seq1, function(x) grepl(".[^-]", as.character(x)))
  s2 <- lapply(file$seq2, function(x) grepl(".[^-]", as.character(x)))
  
  cat("From ", length(file$name), " sequences in input file ", input, ":", "\n", sep="")
  
  if(clean == F){
    cat(sum(s1 == F), " sequences can be removed from ", output_1, "\n", sep="")
    seqinr::write.fasta(sequences=as.list(file$seq1), names=file$name, nbchar=80, file.out=paste(output_1, ".fasta", sep=""))
    cat(sum(s2 == F), " sequences can be removed from ", output_2, "\n", sep="")
    seqinr::write.fasta(sequences=as.list(file$seq2), names=file$name, nbchar=80, file.out=paste(output_2, ".fasta", sep=""))
  } else if (clean == T){
    ss1 <- file[unlist(s1),]
    seqinr::write.fasta(sequences=as.list(ss1$seq1), names=ss1$name, nbchar=80, file.out=paste(output_1, "_clean.fasta", sep=""))
    cat(sum(s1 == F), " sequences were removed from ", output_1, "\n", sep="")
    ss2 <- file[unlist(s2),]
    seqinr::write.fasta(sequences=as.list(ss2$seq2), names=ss2$name, nbchar=80, file.out=paste(output_2, "_clean.fasta", sep=""))
    cat(sum(s2 == F), " sequences were removed from ", output_2, "\n", sep="")
    
  }
  
}   

# concatenation: Concatenate two different files (e.g. 18S and 28S) into a single file (e.g., 18S+28S). Sequences must have an identifier followed by "_"
concatenation <- function(file1, file2, output_file, name="full", order=TRUE, export=TRUE){
  if(grepl("\\.fasta$", file1)){ file1 <- gsub("\\.fasta$","",file1)}
  if(grepl("\\.fasta$", file2)){ file2 <- gsub("\\.fasta$","",file2)}
  if(grepl("\\.fasta$", output_file)){ output_file <- gsub("\\.fasta$","",output_file)}
  
  fasta18 <- seqinr::read.fasta(paste(file1, ".fasta", sep=""), seqtype= "AA", as.string = T)
  fasta28 <- seqinr::read.fasta(paste(file2, ".fasta", sep=""), seqtype= "AA", as.string = T)
  
  fasta18df <- data.frame(name=paste(getAnnot(fasta18)), sequence=paste(fasta18))
  fasta18df$name <- gsub(">", "", fasta18df$name)
  fasta18df$name <- gsub(" ", "_", fasta18df$name)
  fasta18df$ID <- gsub("_.*", "", fasta18df$name) %>% gsub("\\|.*", "", .)
  length(fasta18df$name)==length(unique(fasta18df$name))
  length(fasta18df$ID)==length(unique(fasta18df$ID))
  
  fasta28df <- data.frame(name=paste(getAnnot(fasta28)), sequence=paste(fasta28))
  fasta28df$name <- gsub(">", "", fasta28df$name)
  fasta28df$name <- gsub(" ", "_", fasta28df$name)
  fasta28df$ID <- gsub("_.*", "", fasta28df$name) %>% gsub("\\|.*", "", .)
  length(fasta28df$name)==length(unique(fasta28df$name))
  length(fasta28df$ID)==length(unique(fasta28df$ID))
  
  if(name=="full"){
    fastadf <- merge(fasta18df, fasta28df, by.x = "name", by.y = "name", all = T, sort = F)
    colnames(fastadf) <- c("name", "seq18S", "ID_18S", "seq28S", "ID_28S")        
  }else if(name=="id") {
    if(order==TRUE){
      fastadf <- merge(fasta18df, fasta28df, by.x = "ID", by.y = "ID", all = T, sort = F)
      colnames(fastadf) <- c("ID", "name18S", "seq18S", "name28S", "seq28S")
    } else {
      fastadf <- plyr::join(fasta18df, fasta28df, by=c("ID", "ID"))
      colnames(fastadf) <- c("name18S", "seq18S", "ID", "name28S", "seq28S")
    }
  }else{
    stop("'name' must be either 'full' or 'id' to be able to find the same sequence identifiers in both fasta files.")
  }
  
  
  
  length(fastadf$ID)==length(unique(fastadf$ID))
  
  nchar(as.character(fasta18df$sequence[min(which(!is.na(fastadf$seq18S)))]))
  fastadf$seq18S <- ifelse(is.na(fastadf$seq18S), 
                           paste(rep("-", nchar(as.character(fastadf$seq18S[min(which(!is.na(fastadf$seq18S)))]))), sep="", collapse=""),
                           paste(fastadf$seq18S))
  
  
  nchar(as.character(fasta28df$sequence[min(which(!is.na(fastadf$seq28S)))]))
  fastadf$seq28S <- ifelse(is.na(fastadf$seq28S), 
                           paste(rep("-", nchar(as.character(fastadf$seq28S[min(which(!is.na(fastadf$seq28S)))]))), sep="", collapse=""),
                           paste(fastadf$seq28S))
  
  fastadf$seq <- paste(as.character(fastadf$seq18S), as.character(fastadf$seq28S))
  
  if(name=="full"){
    cat("Used full name for concatenating sequences", "\n", sep="")
    fastadf$sum <- ifelse(!is.na(fastadf$ID_18S)&!is.na(fastadf$ID_28S), "18S + 28S", 
                          ifelse(!is.na(fastadf$ID_18S)&is.na(fastadf$ID_28S), "18S______",
                                 ifelse(is.na(fastadf$ID_18S)&!is.na(fastadf$ID_28S), "______28S",
                                        ifelse(is.na(fastadf$ID_18S)&is.na(fastadf$ID_28S),"OUPS...","OUPS... x2"))))
  } else {
    fastadf$name <- ifelse(is.na(fastadf$name18S), fastadf$name28S, fastadf$name18S)
    fastadf$sum <- ifelse(!is.na(fastadf$name18S)&!is.na(fastadf$name28S), "18S + 28S", 
                          ifelse(!is.na(fastadf$name18S)&is.na(fastadf$name28S), "18S______",
                                 ifelse(is.na(fastadf$name18S)&!is.na(fastadf$name28S), "______28S",
                                        ifelse(is.na(fastadf$name18S)&is.na(fastadf$name28S),"OUPS...","OUPS... x2"))))
  }
  
  if(order==TRUE){
    fastadf <- fastadf[order(fastadf$ID),]
  } else {
  }
  if(export==TRUE){
    seqinr::write.fasta(sequences=as.list(fastadf$seq), names=fastadf$name, nbchar=80, file.out=paste0(output_file, ".fasta"))        
  }else if(export==FALSE){
    return(fastadf)
  }
  rm()    
  
  cat("The concatenation has a total length of ",
      nchar(as.character(fastadf$seq[1])), " bp, and ",
      length(fastadf$seq), " sequences:", "\n", sep="")
  cat("  The 18S has: ", nchar(as.character(fastadf$seq18S[1]))," bp and ", length(grep("18S", fastadf$sum)), " sequences.", "\n", sep="")
  cat("  The 28S has: ", nchar(as.character(fastadf$seq28S[1]))," bp and ", length(grep("28S", fastadf$sum)), " sequences.", "\n", sep="")
}
#################################################################################-
# Fix bug for full name: Error in order(fastadf$ID) : argument 1 is not a vector
#################################################################################-

# Split, trim and concatenate again a fasta file. Requires SplitFasta & concatenation
trim2genes <- function(file, position, gt, st, concatenateBy="full"){
  if(grepl("\\.fasta$", file)) {file<-gsub("\\.fasta$","",file)}
  cat("  Splitting ", file, ".fasta", "\n", sep="")
  SplitFasta(input=file, 
             output_1=paste0("tmp_", file, "_18S"), 
             output_2=paste0("tmp_", file, "_28S"), 
             position, clean=TRUE)
  
  cat("  Trimming files: ", "tmp_", file, "_18S_clean.fasta & ", "tmp_", file, "_28S_clean.fasta", "\n", sep="")
  
  if(is.numeric(st)){
    system(paste0("trimal -in ", "tmp_", file, "_18S_clean.fasta -out ", "tmp_", file, "_18S_clean_trim", gt, "_st.fasta -gt ", gt/100, " -st ", st))
  }else{
    system(paste0("trimal -in ", "tmp_", file, "_18S_clean.fasta -out ", "tmp_", file, "_18S_clean_trim", gt, ".fasta -gt ", gt/100))
  }
  tmp <-  seqinr::read.fasta(paste0("tmp_", file, "_18S_clean_trim", gt, ifelse(is.numeric(st),"_st",""), ".fasta"), seqtype= "AA", as.string = TRUE)
  tmp <- data.frame(name=paste(getAnnot(tmp)), sequence=paste0(tmp))
  tmp$name <- gsub(">", "", tmp$name)
  tmp$name <- gsub(" \\d+ bp$", "", tmp$name)
  tmp$sequence <- gsub(" ", "", tmp$sequence)
  seqinr::write.fasta(sequences=as.list(tmp$sequence), names=tmp$name, 
                      nbchar=80, file.out=paste0("tmp_", file, "_18S_clean_trim", gt, ifelse(is.numeric(st),"_st",""), ".fasta"))
  cat("  ", "tmp_", file, "_18S_clean_trim", gt, ".fasta done", "\n", sep="")
  
  if(is.numeric(st)){
    system(paste0("trimal -in ", "tmp_", file, "_28S_clean.fasta -out ", "tmp_", file, "_28S_clean_trim", gt, "_st.fasta -gt ", gt/100, " -st ", st))
  }else{
    system(paste0("trimal -in ", "tmp_", file, "_28S_clean.fasta -out ", "tmp_", file, "_28S_clean_trim", gt, ".fasta -gt ", gt/100))
  }
  tmp <-  seqinr::read.fasta(paste0("tmp_", file, "_28S_clean_trim", gt, ifelse(is.numeric(st),"_st",""), ".fasta"), seqtype= "AA", as.string = TRUE)
  tmp <- data.frame(name=paste(getAnnot(tmp)), sequence=paste0(tmp))
  tmp$name <- gsub(">", "", tmp$name)
  tmp$name <- gsub(" \\d+ bp$", "", tmp$name)
  tmp$sequence <- gsub(" ", "", tmp$sequence)
  seqinr::write.fasta(sequences=as.list(tmp$sequence), names=tmp$name, 
                      nbchar=80, file.out=paste0("tmp_", file, "_28S_clean_trim", gt, ifelse(is.numeric(st),"_st",""), ".fasta"))
  cat("  ", "tmp_", file, "_28S_clean_trim", gt, ".fasta done", "\n", sep="")
  
  cat("  Concatenating files ", "\n", sep="")
  
  concatenation(file1= paste0("tmp_", file, "_18S_clean_trim", gt, ifelse(is.numeric(st),"_st","")), 
                file2= paste0("tmp_", file, "_28S_clean_trim", gt, ifelse(is.numeric(st),"_st","")), 
                output_file= paste0(file, "_trim", gt, ifelse(is.numeric(st),"_st","")), 
                name=concatenateBy, order=TRUE, export=TRUE)
  
  file.remove(c(paste0("tmp_", file, "_18S_clean.fasta"),
                paste0("tmp_", file, "_28S_clean.fasta"),
                paste0("tmp_", file, "_18S_clean_trim", gt, ifelse(is.numeric(st),"_st",""), ".fasta"),
                paste0("tmp_", file, "_28S_clean_trim", gt, ifelse(is.numeric(st),"_st",""), ".fasta")))
}

# MergeFasta: Add sequences from two or more different files into a single file
MergeFasta <- function(files, output_file, addFileName=TRUE, sep="_", verbose=FALSE){
  file <- c()
  for (i in 1:length(files)){
    fasta <- files[i]
    fasta <- seqinr::read.fasta(fasta, seqtype= "AA", as.string = T)
    fasta <- data.frame(name=paste(getAnnot(fasta)), sequence=paste(fasta))
    fasta$name <- gsub(">", "", fasta$name)
    if(addFileName){
      beg <- files[i]
      beg <- gsub("\\..*", "", beg)
      fasta$name <- paste0(beg, sep, fasta$name)
    }
    file <- rbind(file, fasta)
    if(verbose){
      message(round((i/length(files)*100),2), "%. File '", files[i], "' merged.", sep="")
    }
  }
  seqinr::write.fasta(sequences=as.list(file$seq), names=file$name, nbchar=80, 
                      file.out=output_file)
  cat("Merged file exported to '", output_file, "\n", sep="")
}

# reorderFastaTree: Reorder a fasta file based on the tips of a tree by either the "full" label or only the "id"
reorderFastaTree <- function(fasta, tree, by="full", output){
  if(!missing(output)){
    if(output==fasta){
      stop("'Output' file has the same name as 'input' file. Please rename 'output' file.")
    }
  }
  if(grepl("\\.fasta$", fasta)){ fasta <- gsub("\\.fasta$","",fasta)}
  if(grepl("\\.tre$", tree)){ tree <- gsub("\\.tre$","",tree)}
  
  treeTips <- as.data.frame(ape::read.nexus(paste(tree, ".tre", sep=""))$tip.label)
  colnames(treeTips) <- c("label")
  treeTips$label <- gsub("'", "", treeTips$label)
  
  fastaFile <- seqinr::read.fasta(paste(fasta, ".fasta", sep=""), seqtype= "AA", as.string = T)
  fastaFile <- data.frame(name=paste(seqinr::getAnnot(fastaFile)), sequence=paste(fastaFile))
  fastaFile$name2 <- gsub(">", "", fastaFile$name)
  
  if(by=="full"){
    out <- merge(treeTips, fastaFile, by.x="label", by.y="name2", sort=FALSE, all=TRUE)   
  }else if(by=="id"){
    treeTips$ID <- gsub("_.*|\\|.*|\\..*","", treeTips$label)
    fastaFile$ID <- gsub("_.*|\\|.*|\\..*","", fastaFile$name2)
    out <- merge(treeTips, fastaFile, by="ID", sort=FALSE, all=TRUE)   
  }else{
    stop("Please especify either by='full' or by='id', in order to link the name of the fasta file sequences with the name of the tree tips.")
  }
  if(!exists("out")){
    stop("Fasta file couldn't be merged with tree file.")
  }else{
    if(missing(output)){
      seqinr::write.fasta(sequences=as.list(out$sequence), names=out$label, nbchar=80, file.out=paste0(fasta, "_treeOrder.fasta"))
      cat("Fasta file exported in: '", getwd(), "/", fasta, "_treeOrder.fasta' \n", sep="")
    }else{
      seqinr::write.fasta(sequences=as.list(out$sequence), names=out$label, nbchar=80, file.out=output)
      cat("Fasta file exported in: '", getwd(), "/", output, "' \n", sep="")
    }
    
    if(any(is.na(out$label))){message("Some sequences were named 'NA'. Please check that all sequences in fasta file are in the selected tree.")}
    if(any(is.na(out$sequence))){message("Some sequences were named 'NA'. Please check that all tree tips are in the selected fasta file.")}
  }
}

#
#----
setwd("~/")
rm(list=ls()[!ls() %in% c("MergeFasta", "SplitFasta", "concatenation", "reorderFastaTree", "trim2genes")])
#---- 
#---- Merge ----

# MergeFasta: Add sequences from two or more different files into a single file

(files <- c(grep(".fasta$", dir(), value=T)))
output <- "fasta"

MergeFasta(files=files, 
           output_file=output,
           addFileName=FALSE, sep="_")
#_______________________________________________________________    


#---- Split ----

# SplitFasta: cuts a concatenated file vertically (e.g., 18S+28S) into two different files (e.g. 18S and 28S)

input_file <- "fasta"        ;if(grepl("\\.fasta$", input_file)){ input_file <- gsub("\\.fasta$", "", input_file)}
add <- ""
position <- "1782" # Select the first position of the second file

SplitFasta(input=input_file, 
           output_1=paste0(input_file, "_18S"), 
           output_2=paste0(input_file, "_28S"), 
           position, clean=T)


#---- Concatenate ----

# concatenation: Concatenate two different files (e.g. 18S and 28S) into a single file (e.g., 18S+28S). Sequences must have an identifier followed by "_"
#_______________________________________________________________

concatenation(file1= "fasta_18S", 
              file2= "fasta_28S", 
              output_file= "fasta", 
              name="id", order=F, export=T)


#----
#---- Splitting, Trimming and concatenating    ----

#######################################################-
#### Be careful with the pathways of the files!!!! ####-
#######################################################-

# Trim a file separatedly in two files
trim2genes(file="fasta", 
           position=1920,                           # Position to cut the alignment
           gt=30,                                   # %: -gt -gapthreshold <n>    1 - (fraction of sequences with a gap allowed).
           st=FALSE,                                #    -st -simthreshold <n>    Minimum average similarity allowed.
           concatenateBy="id")                      # Concatenate by "full" name or "id".


#----      
#---- Reordering fasta (Tree Order)          ----


  reorderFastaTree(fasta="fasta",
                 tree="tree",
                 by="id", # by="full" or "id"
                 output="fasta_reOrder")


#


#---- 