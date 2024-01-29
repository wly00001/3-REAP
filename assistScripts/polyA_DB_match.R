
library(reshape2)
library(tidyverse)
library(ggrepel)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomicAlignments)
require(org.Hs.eg.db)
library(plyr)
library(ComplexHeatmap)


args<- commandArgs(TRUE)
args_v<-NA
if(length(args)>1){
   args_v <- args[c(FALSE, TRUE)]
   names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
   print(args_v)
}
args<-NULL

bed_dir <- ifelse(is.na(args_v["bedLAP"]), "_bedLAP_", args_v["bedLAP"])
out_dir <- ifelse(is.na(args_v["out"]), "_out_", args_v["out"])

anno <- read.table("./assistScripts/hg19_polyA_DB_v4.1s.anno", header = TRUE, sep = "\t") #### input polyA_DB annotation

sample <- gsub(".bed.sorted","",tail(strsplit(bed_dir,"/")[[1]],1))
print(paste0("sample name, ", sample))

bed <- read.table(bed_dir, header = F, sep = "\t", quote = "", stringsAsFactors = FALSE)
names(bed) <- c("chr", "start", "end", "readID", "MAPQ", "strand", "CIGAR")

bed$from <- ""
bed$to <- ""
bed$from[bed$strand == "-"] <- bed$start[bed$strand == "-"]
bed$to[bed$strand == "-"] <- bed$end[bed$strand == "-"]
bed$from[bed$strand == "+"] <- 125-bed$end[bed$strand == "+"]
bed$to[bed$strand == "+"] <- 125-bed$start[bed$strand == "+"]
bed$from <- as.numeric(bed$from)
bed$to <- as.numeric(bed$to)





print(paste0("mappable read #, ", length(unique(bed$readID))))

bed$strand_rev <- ""
bed$strand_rev[bed$strand == "-"] <- "+"
bed$strand_rev[bed$strand == "+"] <- "-"

df <- bed

df <- df %>% separate(chr, sep = "::", into = c("PasID","miniREFid"), remove = T)
df <- df %>% separate(PasID, sep = ":", into = c("Chromosome","Strand","Position"), remove = F)

df2 <- subset(df, Strand==strand_rev)
df2$d_LAP2PAS <- abs(df2$to-100)

df3_0 <- df2 %>% 
  group_by(readID) %>% 
  dplyr::slice(which.min(d_LAP2PAS))   ### deal with multimappers, useless for unique reads


df3_0$CIGAR_head_symbol <- sapply(explodeCigarOps(df3_0$CIGAR), head, 1)
df3_0$CIGAR_tail_symbol <- sapply(explodeCigarOps(df3_0$CIGAR), tail, 1)

df3_0$CIGAR_head_value <- sapply(explodeCigarOpLengths(df3_0$CIGAR), head, 1)
df3_0$CIGAR_tail_value <- sapply(explodeCigarOpLengths(df3_0$CIGAR), tail, 1)


startTime <- Sys.time()

df3_0$read_from <-""
df3_0$read_from[df3_0$Strand == "+" & df3_0$CIGAR_tail_symbol=="S"] <- df3_0$CIGAR_tail_value[df3_0$Strand == "+" & df3_0$CIGAR_tail_symbol=="S"]+1
df3_0$read_from[df3_0$Strand == "+" & df3_0$CIGAR_tail_symbol!="S"] <- 1
df3_0$read_from[df3_0$Strand == "-" & df3_0$CIGAR_head_symbol=="S"] <- df3_0$CIGAR_head_value[df3_0$Strand == "-" & df3_0$CIGAR_head_symbol=="S"]+1
df3_0$read_from[df3_0$Strand == "-" & df3_0$CIGAR_head_symbol!="S"] <- 1

df3_0$read_from <- as.numeric(df3_0$read_from)

endTime <- Sys.time()
print(endTime - startTime)


df2D2 <- df3_0 %>% 
  group_by(to,read_from) %>% 
  dplyr::summarise(n = n())

df2D3 <- dcast(df2D2, read_from~to, value.var = "n")
matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$read_from)
matrix2D <- as.matrix(matrix2D)

pdf(paste0(out_dir,"_heatmap.pdf"), width=10,height=10)
Heatmap(matrix2D, name = sample,
        row_order = rownames(matrix2D),
        column_order = colnames(matrix2D))
dev.off()
############

df3_1 <- subset(df3_0, d_LAP2PAS<=24)
print(paste0("PAS matched (PASS) # (<=24), ", nrow(df3_1)))

df3 <- subset(df3_1, read_from<=6)
print(paste0("PAS matched (PASS) # (<=24 & <=6), ", nrow(df3)))

df3$Position <- as.numeric(df3$Position)


PASS_bw_bed <- df3[c("Chromosome", "Position", "start", "end", "readID", "MAPQ", "strand")]

PASS_bw_bed$start_GENOME <- ""
PASS_bw_bed$end_GENOME <- ""
PASS_bw_bed$start_GENOME <- as.numeric(PASS_bw_bed$start_GENOME)
PASS_bw_bed$end_GENOME <- as.numeric(PASS_bw_bed$end_GENOME)
PASS_bw_bed$end_GENOME[PASS_bw_bed$strand == "-"] <- PASS_bw_bed$end[PASS_bw_bed$strand == "-"]-100+PASS_bw_bed$Position[PASS_bw_bed$strand == "-"]
PASS_bw_bed$start_GENOME[PASS_bw_bed$strand == "-"] <- PASS_bw_bed$end_GENOME[PASS_bw_bed$strand == "-"]-1
PASS_bw_bed$start_GENOME[PASS_bw_bed$strand == "+"] <- PASS_bw_bed$start[PASS_bw_bed$strand == "+"]-25+PASS_bw_bed$Position[PASS_bw_bed$strand == "+"]
PASS_bw_bed$end_GENOME[PASS_bw_bed$strand == "+"] <- PASS_bw_bed$start_GENOME[PASS_bw_bed$strand == "+"]+1

PASS_bw_bed <- PASS_bw_bed[c("Chromosome", "start_GENOME", "end_GENOME", "readID", "MAPQ", "strand")]

options(scipen = 999)
write.table(PASS_bw_bed, paste0(out_dir,"_PASS_bw.bed"), col.names = F, row.names = F, sep = "\t", quote = F)






pdf(paste0(out_dir,"_distribution.pdf"), width=10,height=5); layout(matrix(1:5,nrow=5))
ggplot(df3_0, aes(x=to)) + 
  geom_histogram()

ggplot(df3_0, aes(x=read_from)) + 
  geom_histogram()

ggplot(df3_1, aes(x=read_from)) + 
  geom_histogram()

dev.off()






write.table(df3[c("PasID","readID")], paste(out_dir, "RoSpm.tbl", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)

cluster.all.reads <- data.frame(table(df3$PasID))
print(paste0("PASS# in csv, ", sum(cluster.all.reads$Freq)))

names(cluster.all.reads) <- c("PasID",sample)

write.csv(cluster.all.reads, paste0(out_dir,"_cluster.all.reads.csv"), col.names = T, row.names = F, quote = F)





