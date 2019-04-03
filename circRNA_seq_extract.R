library("Biostrings")
library("stringi")
library("seqRFLP")

args <- commandArgs(trailingOnly = TRUE)
circSeq_1<-function (fasta_file, e_count, e_sizes, e_offsets, strand, circ_type, out_filename){
  C=NULL
  for (i in 1:nrow(fasta_file)){
    a<-fasta_file[i,2]
    a<-as.character(a)
    if(strand[i]=="-" & circ_type[i]!="intergenic" ) a<-chartr("acgtACGT", "tgcaTGCA", a)
    B=NULL
    pp<-do.call("rbind", strsplit(as.character(e_sizes[i]),","))
    ppp<-data.frame(apply(pp, 2, as.numeric))
    pp1<-do.call("rbind", strsplit(as.character(e_offsets[i]),","))
    ppp1<-data.frame(apply(pp1, 2, as.numeric))
    for (j in 1:e_count[i]){
      B[j]<-substr(a, ppp1[j,1]+1, ppp1[j,1]+ppp[j,1])
    }
    C[i]<-paste(B,collapse="")
    if(strand[i]=="-" & circ_type[i]!="intergenic") C[i]<-stri_reverse(C[i])
  }
  csv<-data.frame(fasta_file[,1],C)
  colnames(csv)<-c("name", "seq")
  df.fasta = dataframe2fas(csv, out_filename)
}

fastatodframe <- function(filename){
dnaString<-readDNAStringSet(filename)
seq_name = names(dnaString)
sequence = paste(dnaString)
df <- data.frame(seq_name, sequence)
return(df)
}


fasta_file<-fastatodframe(args[1]) # args[1] : circRNA_genomic_sequence.fasta

txt<-read.table(args[2], header=T )# genpmic interval(start, end) and args[2] : circRNA_class.txt

circSeq_1(fasta_file, txt$e_count, txt$e_sizes, txt$e_offsets, txt$circ_strand, txt$circ_type, file.path(args[3], 'circRNA_sequence.fasta'))#extract circRNA sequence and args[3] : output directory



