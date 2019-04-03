
args <- commandArgs(trailingOnly = TRUE)

transcript<-read.table(args[1],header=T) # args[1]: transcript data file

output_circ_tool<-read.table(args[2]) # args[2]: bed file from the output of circRNA prediction tool#
Splice_length=NULL
CircRNA_type=NULL
b_transcript=NULL
b_strand=NULL
b_trans_start=NULL
b_trans_end=NULL
b_gene=NULL
e_count=NULL
e_sizes=NULL
e_offsets=NULL
for(i in 1:dim(output_circ_tool)[1]){
s=as.numeric(as.character(output_circ_tool$V2[i]))
e=as.numeric(as.character(output_circ_tool$V3[i]))
wh<-which(as.character(transcript$chr)==as.character(output_circ_tool$V1[i]) & transcript$trans_start<e & transcript$trans_end>s)
if (length(wh)==0){
circRNA_type="intergenic"
splice_length=e-s
B_transcript="NA"
B_strand="NA"
B_trans_start="NA"
B_trans_end="NA"
B_gene="NA"
e_count_1=1
e_sizes_1=e-s
e_offsets_1=0
}else{
trans<-transcript[wh,]
splice_L=NULL
dm=NULL
convse_1=NULL
t_length=NULL
for(j in 1: dim(trans)[1]){
a<-do.call("rbind", strsplit(as.character(trans$exon_starts[j]),","))
a1<-data.frame(apply(a, 2, as.numeric))
b<-do.call("rbind", strsplit(as.character(trans$exon_ends[j]),","))
b1<-data.frame(apply(b, 2, as.numeric))
colnames(a1)<-"start"
colnames(b1)<-"end"
exon_start_end<-cbind(a1,b1)

wh1<-exon_start_end[which(exon_start_end$start<e & exon_start_end$end>s),]
if (dim(wh1)[1]==0){
convse_1[j]<-1000000
splice_L[j]<-e-s
}else if(wh1$start[1]==s & wh1$end[dim(wh1)[1]]==e){
convse_1[j]<-0
splice_L[j]<-sum(wh1$end-wh1$start)
}else if (wh1$start[1]==s | wh1$end[dim(wh1)[1]]==e){
convse_1[j]<-0.5
wh1$start[1]=s
wh1$end[dim(wh1)[1]]=e
splice_L[j]<-sum(wh1$end-wh1$start)
}else{
convse_1[j]<-(abs(wh1$start[1]-s)+abs(wh1$end[dim(wh1)[1]]-e))/2
wh1$start[1]=s
wh1$end[dim(wh1)[1]]=e
splice_L[j]<-sum(wh1$end-wh1$start)
}
dm[j]=dim(wh1)[1]
t_length[j]=trans$trans_end[j]-trans$trans_start[j]
}

d_l<-data.frame(dm,splice_L,convse_1,t_length)

if (min(d_l$convse_1)==0){
d_2<-d_l[which(d_l$convse_1==min(d_l$convse_1)),]
d_3<-d_2[which(d_2$splice_L==max(d_2$splice_L)),]
d_4<-d_3[which(d_3$t_length==max(d_3$t_length)),]
}else {

d_2<-d_l[which(d_l$dm==max(d_l$dm)),]
d_3<-d_2[which(d_2$t_length==max(d_2$t_length)),]
d_4=d_3[which(d_3$splice_L==max(d_3$splice_L)),]
}

tr<-as.numeric(as.character(row.names(d_4)))

trans_1<-trans[tr,]

trans_3=trans_1

if(dim(trans_3)[1]>1) trans_3=trans_3[order(trans_3$transcript_id),]
trans_3=trans_3[1,]

B_transcript=as.character(trans_3$transcript_id)
B_strand=as.character(trans_3$strand)
B_trans_start=as.character(trans_3$trans_start)
B_trans_end=as.character(trans_3$trans_end)
B_gene=as.character(trans_3$gene)
pp<-do.call("rbind", strsplit(as.character(trans_3$exon_starts),","))
ppp<-data.frame(apply(pp, 2, as.numeric))
pp1<-do.call("rbind", strsplit(as.character(trans_3$exon_ends),","))
ppp1<-data.frame(apply(pp1, 2, as.numeric))
colnames(ppp)<-"start"
colnames(ppp1)<-"end"
exon_start_end<-cbind(ppp,ppp1)

wh1<-exon_start_end[which(exon_start_end$start<e & exon_start_end$end>s),]
if (as.character(B_strand)!=as.character(output_circ_tool$V4[i])){
circRNA_type="antisense"
if (dim(wh1)[1]>0){
wh1$start[1]=s
wh1$end[dim(wh1)[1]]=e
splice_length<-sum(wh1$end-wh1$start)
e_count_1=dim(wh1)[1]
e_sizes_1=wh1$end-wh1$start
e_offsets_1=wh1$start-s
}else{
splice_length<-e-s
e_count_1=1
e_sizes_1=e-s
e_offsets_1=0
}
}else if (as.character(B_strand)==as.character(output_circ_tool$V4[i]) & dim(wh1)[1]==0){
circRNA_type="intronic"
splice_length=e-s
e_count_1=1
e_sizes_1=e-s
e_offsets_1=0
}else if (as.character(B_strand)==as.character(output_circ_tool$V4[i]) & dim(wh1)[1]>0 & wh1$start[1]==s & wh1$end[dim(wh1)[1]]==e){
circRNA_type="exonic"
splice_length<-sum(wh1$end-wh1$start)
e_count_1=dim(wh1)[1]
e_sizes_1=wh1$end-wh1$start
e_offsets_1=wh1$start-s
}else{
circRNA_type="sense_overlapping"
wh1$start[1]=s
wh1$end[dim(wh1)[1]]=e
splice_length<-sum(wh1$end-wh1$start)
e_count_1=dim(wh1)[1]
e_sizes_1=wh1$end-wh1$start
e_offsets_1=wh1$start-s
}
}
Splice_length[i]=splice_length
CircRNA_type[i]=circRNA_type
b_transcript[i]=B_transcript
b_strand[i]=B_strand
b_trans_start[i]=B_trans_start
b_trans_end[i]=B_trans_end
b_gene[i]=B_gene
e_count[i]=e_count_1
e_sizes[i]=paste(e_sizes_1, collapse=",")
e_offsets[i]=paste(e_offsets_1, collapse=",")
}
cRNA_id<-paste(output_circ_tool$V1,output_circ_tool$V2, sep=':')
cRNA_id_final<-paste(cRNA_id,output_circ_tool$V3, sep='-')

circ_class<-cbind(cRNA_id_final, output_circ_tool[,1:4],Splice_length, CircRNA_type, e_count, e_sizes, e_offsets, b_transcript, b_strand, b_trans_start, b_trans_end, b_gene)
colnames(circ_class)<-c("ID", "chr", "circ_start", "circ_end", "circ_strand", "splice_L", "circ_type", "e_count", "e_sizes", "e_offsets", "b_transcript", "b_strand", "b_trans_start", "b_trans_end", "b_gene")
write.table(circ_class,file.path(args[3], 'circRNA_class.txt'),sep="\t",quote=F,row.names=F) # args[3]: output directory#

circ_class_bed<-circ_class[,2:4]
write.table(circ_class_bed, file.path(args[3], 'circRNA_class.bed'), sep="\t",quote=F,row.names=F,col.names=F) # args[2]: output directory#




