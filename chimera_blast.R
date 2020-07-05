#---- 
#---- loading packages ----

library(seqinr)

library(data.table)

#---- 
setwd("~/")
#---- 
file <- "new.fasta"                  
length <- 300
threshold <- 20
#---- 
#---- Opening data and exporting dataset for BLAST ----

if(grepl("\\.fasta$", file)){ file <- gsub("\\.fasta$","",file)}

fasta <- seqinr::read.fasta(paste(file, ".fasta", sep=""), seqtype= "AA", as.string = T)
fasta <- data.frame(name=gsub(">", "", paste(seqinr::getAnnot(fasta))), sequence=paste(fasta))


fasta$id <- gsub("\\|.*","", fasta$name) %>% gsub("\\..*", "", .) %>% gsub(" .*", "", .) 


# removing duplicates
fasta$id[duplicated(fasta$id)]


# Creating the fasta file that is going to be blasted
out <- data.frame()
for(i in unique(fasta$id)){
    ss <- subset(fasta, id==i)
    tmp <- data.frame(name=c(paste0(i, "_begining"), paste0(i, "_end")),
                      sequence=c(substr(ss$sequence, 1, length), 
                                 substr(ss$sequence, nchar(as.character(ss$sequence))-length+1, nchar(as.character(ss$sequence)))))
    out <- rbind(out, tmp)
    if(nchar(as.character(ss$sequence))<(length*2)){
        warning("Sequence '", i, "' has a length lower than ", length*2, 
                " (", nchar(as.character(ss$sequence))," bp). Begining and end are overlapping.")
        }
    if(grep(i, fasta$id) %% round(nrow(fasta)*0.1) == 0){
        message(round(grep(i, fasta$id)/nrow(fasta),1)*100, "%")
        }
}; rm(ss, tmp, i)


# Exporting
dir.create("chimeras")
seqinr::write.fasta(sequences=as.list(out$sequence), names=out$name, nbchar=80, 
                    file.out=paste0("chimeras/", file, "_toBlast.fasta"))



#---- 
#---- BLAST ----

# Submitting the BLAST as a job:

# To automatically create the shell script to submit the job:
# But remember to check the options and paths !!!!
sh <- data.frame(c("#!/bin/bash", 
                   "#$ -S /bin/bash",
                   "#$ -V",
                   "#$ -cwd",
                   "#$ -q long.q",
                   "#$ -M XXXXXXXXXXXXXXXX",
                   "#$ -m beas",
                   "#$ -N blastChim",
                   "",
                   paste0('FASTA="',file, '_toBlast.fasta"'),
                   "",
                   'BLAST_TSV="${FASTA/_toBlast.fasta/_blast.tsv}"',
                   'OUT_FMT="6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore"',
                   "",
                   '/usr/local/genome2/conda3/envs/blast-2.9.0/bin/blastn -num_threads 8 -max_target_seqs 100 -evalue 1.00e-10 -query $FASTA -out $BLAST_TSV -db /db/nt/nt_2020-2-26/flat/nt -outfmt "$OUT_FMT"'))

write.table(sh, paste0("chimeras/", file, "_toBlast.sh"), quote=FALSE, sep='', row.names=FALSE, col.names=FALSE, eol="\n")

# Send this files to cluster
# Again remember to check the paths !!!!
system(paste0("rsync -r", 
              " chimeras/", file, "_toBlast.fasta",
              " chimeras/", file, "_toBlast.sh",
              " XXXX_cluster_Address_XXXX"))

# Now importing the output back in the PC 
# Again remember to check the paths !!!!
system(paste0("rsync ", 
              "XXXX_cluster_Address_XXXX", file, "_blast.tsv ",
              "XXXX_working_Address_XXXX"))

#

#---- 
#---- Opening blast results and exploring matches ----

df <- fread(paste0("chimeras/", file, "_blast.tsv"))
colnames(df) <- c("qseqid", "sseqid", "sacc", "stitle", "sscinames", "staxids", "sskingdoms", "sblastnames", "pident", "slen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  
head(df)


df$id <- gsub("_.*", "", df$qseqid)
df$pos <- gsub(".*_", "", df$qseqid)

chim <- data.frame(sequence=c(), chimera=c(), matchBoth=c(), comment=c())
chimeras <- c()
for(i in unique(df$id)){
    sb <- subset(df, id==i & pos=="begining")
    se <- subset(df, id==i & pos=="end")
    
    tmp <- data.frame()
    if(sum(sb$sacc %in% se$sacc)<=threshold | sum(se$sacc %in% sb$sacc)<=threshold){
        cat("Sequence '", i, "' has ", threshold, " or less blast matches in common.", sep="", "\n")
        tmp <- data.frame(sequence=i, chimera="maybe", 
                          matchBoth=sum(min(sb$sacc %in% se$sacc, se$sacc %in% sb$sacc)), 
                          comment="")
        if(sum(sb$sacc %in% se$sacc)==0 | sum(se$sacc %in% sb$sacc)==0){
            cat("    '", i, "' has no similar sequences in all matches. Detected as chimera.", sep="", "\n")
            tmp <- data.frame(sequence=i, chimera="yes", 
                              matchBoth="0", 
                              comment="")
            chimeras <- c(chimeras, i)
        }
        else if((sum(sb$sacc %in% se$sacc)==1 & sb$sacc[sb$sacc %in% se$sacc]==i) | (sum(se$sacc %in% sb$sacc)==1 & se$sacc[se$sacc %in% sb$sacc]==i)){
            cat("    '", i, "' has one similar match that corresponds with itself. Detected as chimera.", sep="", "\n")
            tmp <- data.frame(sequence=i, chimera="yes", 
                              matchBoth="1", 
                              comment="Self-match")
            chimeras <- c(chimeras, i)
        }
        else if(all(sb$sacc[sb$sacc %in% se$sacc] %in% chimeras) | 
                all(se$sacc[se$sacc %in% sb$sacc] %in% chimeras)){
            cat("        '", i, "' matches with previously detected chimeras only!!!!", sep="", "\n")
            tmp <- data.frame(sequence=i, chimera="yes", 
                              matchBoth=sum(min(sb$sacc %in% se$sacc, se$sacc %in% sb$sacc)), 
                              comment=paste("Match_only_with_other_chimeras!: ", 
                                            paste(chimeras[chimeras %in% sb$sacc[sb$sacc %in% se$sacc]], collapse=","), sep=""))
            chimeras <- c(chimeras, i)
        }
        else if(any(sb$sacc[sb$sacc %in% se$sacc] %in% chimeras) & !all(sb$sacc[sb$sacc %in% se$sacc] %in% chimeras) |
                any(se$sacc[se$sacc %in% sb$sacc] %in% chimeras) & !all(se$sacc[se$sacc %in% sb$sacc] %in% chimeras)){
            cat("        '", i, "' some matches with previously detected chimeras!", sep="", "\n")
            tmp <- data.frame(sequence=i, chimera="maybe", 
                              matchBoth=sum(sb$sacc %in% se$sacc), 
                              comment=paste("Match_with_", length(chimeras[chimeras %in% sb$sacc[sb$sacc %in% se$sacc]]) ,"_other_chimeras: ", 
                                        paste(chimeras[chimeras %in% sb$sacc[sb$sacc %in% se$sacc]], collapse=","), sep=""))
        }
    }
    chim <- rbind(chim, tmp)

    if(length(unique(df$id))<=10){
        message("Reading")
        } else {
            if (grep(i, unique(df$id)) %% round(length(unique(df$id))*0.1) == 0){
                message(round(grep(i, unique(df$id))/length(unique(df$id)),1)*100, "%")
            }
        }
}; rm(i, sb, se, tmp)



subset(chim, chimera=="yes")
cat("Identified ", nrow(chim), " (", round(nrow(chim)/length(unique(df$id))*100, 2), "%) sequences as possible chimeras,", "\n", 
    "  of which ", sum(chim$chimera=="yes"), " are most likely chimeras. \n", sep="")

write.table(chim, paste0("chimeras/", file, "_chimeras.tsv"), quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)


# ----



