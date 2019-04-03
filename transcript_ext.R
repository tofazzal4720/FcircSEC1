################################################################################################

args <- commandArgs(trailingOnly = TRUE)

getFileNameExtension <- function (fn) {
# remove a path
splitted    <- strsplit(x=fn, split='/')[[1]]   
# or use .Platform$file.sep in stead of '/'
fn          <- splitted [length(splitted)]
ext         <- ''
splitted    <- strsplit(x=fn, split='\\.')[[1]]
l           <-length (splitted)
if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l] 
# the extention must be the suffix of a non-empty name    
ext
}

ext<-getFileNameExtension(args[1])

if(ext=="gff" | ext=="gff3"){
annot_file<- read.delim(args[1], header = FALSE, sep = '\t', skip = 8)
} else if(ext=="gtf") {
annot_file<-read.table(args[1], header = FALSE, sep = '\t')
}else {
print("Error: please input the annotation file in gff or gtf format as the first argument")
}



if(args[2]=="ncbi"){
  gff_exon<-annot_file[which(annot_file$V3=="exon"),]

  trans_id<-sub(".*transcript_id=", "", gff_exon$V9)
  if (length(which(grepl("ID=",trans_id)==T))==length(trans_id) | length(which(grepl("gene_id ",trans_id)==T))==length(trans_id)){
  print("Error: your provided annoation file is not an ncbi annotation file")
  }else{
  gene_id<-sub('.*gene=', '', gff_exon$V9)
  gene_id_final<-sub('\\;.*', '', gene_id)
  gff_exon_trans_id<-data.frame(gff_exon,trans_id, gene_id_final)

  wh<-which(grepl("ID=",trans_id)==T)

  gff_exon_trans_id_final<-gff_exon_trans_id[-wh,]


  trans_chr_id<-paste(gff_exon_trans_id_final$trans_id, gff_exon_trans_id_final$V1, sep=':')

  gff_exon_trans_chr_id_final<-data.frame(gff_exon_trans_id_final,trans_chr_id)

  gff_exon_final<-gff_exon_trans_chr_id_final

  gff_exon_final<-gff_exon_final[-9]

  gff_exon_final$V4<-gff_exon_final$V4-1

  colnames(gff_exon_final)<-c("V1","V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11")

  trans_chr_id_1<-gff_exon_final$V11
  dup<-duplicated(trans_chr_id_1)
  trans_chr_id_uni<-trans_chr_id_1[-which(dup==T)]

  transcript_id=NULL
  chr=NULL
  strand=NULL
  trans_start=NULL
  trans_end=NULL
  exon_count=NULL
  exon_starts=NULL
  exon_ends=NULL
  gene=NULL

  for(i in 1:length(trans_chr_id_uni)){
  total_id<-gff_exon_final[which(gff_exon_final$V11==trans_chr_id_uni[i]),]
  exon_id<-data.frame(total_id$V4,total_id$V5)
  exon_id_sort<-exon_id[order(exon_id[,1]),]
  transcript_id[i]=as.character(total_id$V9[1])
  chr[i]=as.character(total_id$V1[1])
  strand[i]=as.character(total_id$V7[1])
  trans_start[i]=exon_id_sort[1,1]
  trans_end[i]=exon_id_sort[dim(exon_id_sort)[1],2]
  exon_count[i]=dim(exon_id_sort)[1]
  exon_starts[i]=paste(exon_id_sort[,1],collapse=",")
  exon_ends[i]=paste(exon_id_sort[,2],collapse=",")
  gene[i]=as.character(total_id$V10[1])
  }
  trans_data<-data.frame(transcript_id, chr, strand, trans_start, trans_end, exon_count, exon_starts, exon_ends,gene)

  write.table(trans_data,file.path(args[3], 'transcriptdata.txt'),sep="\t",quote=F,row.names=F)
  }
}else if (args[2]=="ucsc"){
  gtf_exon<-annot_file[which(annot_file$V3=="exon"),]

  trans_id<-sub(".*transcript_id ", "", gtf_exon$V9)
  if(length(which(grepl("ID=",trans_id)==T))==length(trans_id) | length(which(grepl("gene_id ",trans_id)==T))==length(trans_id)){
  print("Error: your provided annoation file is not an ucsc annotation file")
  }else{
  trans_id_final<-sub('\\;.*', '', trans_id)

  gene_id<-sub('.*gene_name ', '', gtf_exon$V9)
  gene_id_final<-sub('\\;.*', '', gene_id)

  gtf_exon_trans_id<-data.frame(gtf_exon,trans_id_final, gene_id_final)

  wh<-which(grepl("gene_id",trans_id_final)==T)

  if (length(wh)>0)gtf_exon_trans_id_final<-gtf_exon_trans_id[-wh,] else gtf_exon_trans_id_final<-gtf_exon_trans_id

  trans_chr_id<-paste(gtf_exon_trans_id_final$trans_id_final, gtf_exon_trans_id_final$V1, sep=':')

  gtf_exon_trans_chr_id_final<-data.frame(gtf_exon_trans_id_final,trans_chr_id)

  gtf_exon_final<-gtf_exon_trans_chr_id_final

  gtf_exon_final<-gtf_exon_final[-9]

  gtf_exon_final$V4<-gtf_exon_final$V4-1

  colnames(gtf_exon_final)<-c("V1","V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11")

  trans_chr_id_1<-gtf_exon_final$V11
  length(trans_chr_id_1)
  dup<-duplicated(trans_chr_id_1)
  trans_chr_id_uni<-trans_chr_id_1[-which(dup==T)]


  transcript_id=NULL
  chr=NULL
  strand=NULL
  trans_start=NULL
  trans_end=NULL
  exon_count=NULL
  exon_starts=NULL
  exon_ends=NULL
  gene=NULL

  for(i in 1:length(trans_chr_id_uni)){
  total_id<-gtf_exon_final[which(gtf_exon_final$V11==trans_chr_id_uni[i]),]
  exon_id<-data.frame(total_id$V4,total_id$V5)
  exon_id_sort<-exon_id[order(exon_id[,1]),]
  transcript_id[i]=as.character(total_id$V9[1])
  chr[i]=as.character(total_id$V1[1])
  strand[i]=as.character(total_id$V7[1])
  trans_start[i]=exon_id_sort[1,1]
  trans_end[i]=exon_id_sort[dim(exon_id_sort)[1],2]
  exon_count[i]=dim(exon_id_sort)[1]
  exon_starts[i]=paste(exon_id_sort[,1],collapse=",")
  exon_ends[i]=paste(exon_id_sort[,2],collapse=",")
  gene[i]=as.character(total_id$V10[1])
  }
  trans_data<-data.frame(transcript_id, chr, strand, trans_start, trans_end, exon_count, exon_starts, exon_ends,gene)

  write.table(trans_data,file.path(args[3], 'transcriptdata.txt'),sep="\t",quote=F,row.names=F)
  }
}else if (args[2]=="other"){
  gtf_exon<-annot_file[which(annot_file$V3=="exon"),]

  trans_id<-sub(".*transcript_id ", "", gtf_exon$V9)
  if(length(which(grepl("ID=",trans_id)==T))==length(trans_id) | length(which(grepl("gene_id ",trans_id)==T))==length(trans_id)){
  print("Error: your provided annoation file is not supported by our method, please use another annotation file")
  }else{
  trans_id_final<-sub('\\;.*', '', trans_id)

  gene_id<-sub('.*gene_name ', '', gtf_exon$V9)
  gene_id_final<-sub('\\;.*', '', gene_id)

  gtf_exon_trans_id<-data.frame(gtf_exon,trans_id_final, gene_id_final)

  wh<-which(grepl("gene_id",trans_id_final)==T)

  if (length(wh)>0)gtf_exon_trans_id_final<-gtf_exon_trans_id[-wh,] else gtf_exon_trans_id_final<-gtf_exon_trans_id

  trans_chr_id<-paste(gtf_exon_trans_id_final$trans_id_final, gtf_exon_trans_id_final$V1, sep=':')

  gtf_exon_trans_chr_id_final<-data.frame(gtf_exon_trans_id_final,trans_chr_id)

  gtf_exon_final<-gtf_exon_trans_chr_id_final

  gtf_exon_final<-gtf_exon_final[-9]

  gtf_exon_final$V4<-gtf_exon_final$V4-1

  colnames(gtf_exon_final)<-c("V1","V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11")

  trans_chr_id_1<-gtf_exon_final$V11
  length(trans_chr_id_1)
  dup<-duplicated(trans_chr_id_1)
  trans_chr_id_uni<-trans_chr_id_1[-which(dup==T)]


  transcript_id=NULL
  chr=NULL
  strand=NULL
  trans_start=NULL
  trans_end=NULL
  exon_count=NULL
  exon_starts=NULL
  exon_ends=NULL
  gene=NULL

  for(i in 1:length(trans_chr_id_uni)){
  total_id<-gtf_exon_final[which(gtf_exon_final$V11==trans_chr_id_uni[i]),]
  exon_id<-data.frame(total_id$V4,total_id$V5)
  exon_id_sort<-exon_id[order(exon_id[,1]),]
  transcript_id[i]=as.character(total_id$V9[1])
  chr[i]=as.character(total_id$V1[1])
  strand[i]=as.character(total_id$V7[1])
  trans_start[i]=exon_id_sort[1,1]
  trans_end[i]=exon_id_sort[dim(exon_id_sort)[1],2]
  exon_count[i]=dim(exon_id_sort)[1]
  exon_starts[i]=paste(exon_id_sort[,1],collapse=",")
  exon_ends[i]=paste(exon_id_sort[,2],collapse=",")
  gene[i]=as.character(total_id$V10[1])
  }
  trans_data<-data.frame(transcript_id, chr, strand, trans_start, trans_end, exon_count, exon_starts, exon_ends,gene)

  write.table(trans_data,file.path(args[3], 'transcriptdata.txt'),sep="\t",quote=F,row.names=F)
  }
}else{
print("Error: please input 'ncbi' or 'ucsc' or 'other' as the second argument")
} 