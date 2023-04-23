# library(ShortRead)
.libPaths('/home/nishuai/R/x86_64-pc-linux-gnu-library/4.2/')

invisible(capture.output(suppressPackageStartupMessages(library(Biostrings))))
invisible(capture.output(suppressPackageStartupMessages(library(msa))))
invisible(capture.output(suppressPackageStartupMessages(library(tidyr))))
invisible(capture.output(suppressPackageStartupMessages(library(writexl))))
invisible(capture.output(suppressPackageStartupMessages(library(dplyr))))
invisible(capture.output(suppressPackageStartupMessages(library(getopt))))
invisible(capture.output(suppressPackageStartupMessages(library(readxl))))
invisible(capture.output(suppressPackageStartupMessages(library(igraph))))

# spec <- matrix(
#   c("outdir",  "o", 1, "character", "outdir",
#     "indir", "i", 1, "character", "indir",
#     "tagtablefile", "t", 1, "character",  "tagtablefile",
#     "isovertimefile", "l", 1, "character","isovertimefile",
#     "percent", "p", 2, "numeric", "percent ratio of topIS"),
#   byrow=TRUE, ncol=5)
spec <- matrix(
  c("config",  "c", 1, "character", "config",
    "samfile","s",1,"character","samfile"),
  byrow=TRUE, ncol=5)

opt=getopt(spec=spec)
config=opt$config
samfile=opt$samfile
print(c(config,samfile))

config_data=read.delim(config,sep="=",header=F,as.is=T)

tagtablefile=config_data[match("tagtablefile",config_data$V1),'V2'] #*** tagtablefile文件路径
isovertimefile=config_data[match("top_is_file",config_data$V1),'V2']  #*** top_is_file文件路径,
percent=as.numeric(config_data[match("select_is_percent",config_data$V1),'V2'])  #*** 选择top_is_percent比例信息
indir=config_data[match("indir",config_data$V1),'V2']   #***输入路径 
outdir=config_data[match("outdir",config_data$V1),'V2'] #***输出路径
genome_fa = config_data[match("genome_fa", config_data$V1), "V2"]
print(c(tagtablefile,isovertimefile,percent,indir,outdir,genome_fa))


#*******************************************判断文件是否存在

#*****************************************提取genome序列信息
# genome_fa="D:/项目/工作/项目/MiS026/zea/GCF_000005005.2_B73_RefGen_v4_genomic.changechromome.fa"
# genome_fa="/mount/labdatastorage001/input/konghl/genome/zea_GCF_000005005.2/ncbi_dataset/data/GCF_000005005.2/GCF_000005005.2_B73_RefGen_v4_genomic.changechromome.fa"

genome = readDNAStringSet(genome_fa, "fasta")

#******** 提取其他的fa序列信息, 根据tagtable表进行读取，并且进行is和is.r文件的合并
#*然后再进行提取序列信息等再做序列比对信息，会更方便一些
# tagtablefile="/mount/labdatastorage001/input/konghl/data/plant/MiS026/zea/A001_1/tagTabMiS026.txt"
# indir="/mount/labdatastorage001/input/konghl/data/plant/MiS026/zea/A001_1/results_A001_1"
# outdir=indir
# isot_path="ISOT_results_MiS026_A001_1"
# isovertimefile=file.path(indir,isot_path,'isovertime_results_by_replicate_MiS026_A001_1.xlsx')
# # num=20
# percent=0.1

tagtable=read.delim(tagtablefile,sep="\t",as.is=T,header=F,check.names=F)
topis=data.frame(read_excel(file.path(isovertimefile)),check.names = F,stringsAsFactors = F)

#samples=sapply(tagtable$V4,function(t){unlist(strsplit(t,split="_"))[5]})
samples=sapply(strsplit(colnames(topis)[7],split="x"),function(t)t[2])

tagtable=tagtable[grep(paste0(samples,"x"),tagtable$V4),]

is_r1_list=paste0(tagtable$V4,'_._25.is.r1')
is_list=paste0(tagtable$V4,'_._25.is')
is=c()
is_r1=c()

for(i in 1:nrow(tagtable)){
  tmp_isfile=file.path(indir,'is_file',is_list[i])
  tmp_isr1file=file.path(indir,'is_file',is_r1_list[i])
  print(c(tmp_isfile,tmp_isr1file))
  if(file.exists(tmp_isfile)){
    # print(tmp_isfile)  
    tmp_is=read.delim(tmp_isfile,sep=" ",as.is=T,header=F)
    if(nrow(tmp_is)>=1){
      tmp_is$sample=samples[1]
    }
    is=rbind(is,tmp_is)  
  }
  if(file.exists(tmp_isr1file)){
    # print(tmp_isr1file)
    tmp_isr1=read.delim(tmp_isr1file,sep=" ",as.is=T,header=F)
    if(nrow(tmp_isr1)>=1){
      tmp_isr1$sample=samples[1]
    }
    is_r1=rbind(is_r1,tmp_isr1)
  }
}
if(!is.null(is)){
  is$type="is"
}
if(!is.null(is_r1)){
  is_r1$type="is_r1"
}
total_is=rbind(is,is_r1)

if(nrow(topis)>=4){
  topis=topis[4:nrow(topis),,drop=FALSE]
}
which_id=grep('MegaPCR',colnames(topis))
for(i in 7:(which_id-1)){
  topis[,i]=as.numeric(unlist(topis[,i]))
}
mean_value=rowMeans(topis[,7:(which_id-1),drop=FALSE])
topis=topis[order(topis[,7],decreasing = T),] ###排序，按照value值排序
num=length(which(mean_value>=percent))
num=min(num,nrow(topis))

topis=as.data.frame(topis[1:num,])

#*************************************读取sam文件获取比对上的序列信息
# set the path to your SAM file

sam_list=read.delim(samfile,sep="\t",as.is=T,header=T,check.names=F,stringsAsFactors = F)
sam_name_split=sapply(strsplit(sam_list$sample,split="_"),function(t)t[5])
sam_list$sample=sam_name_split

#*************************************添加检查，发现sam_list和total_is是否有错误
print("#######检查sam_list 和total_is 之间的交集信息")
print("sam_list\t total_is\t length(intersect(sam_list[,1],total_is[,1]))\n")
print(c(dim(sam_list),dim(total_is),length(intersect(sam_list[,1],total_is[,1]))))
print("########检查sam_list 和total_is 之间的交集结束")

substitutionMatrix=matrix(c(1,-0.1,-0.1,-0.1,
                            -0.1,1,-0.1,-0.1,
                            -0.1,-0.1,1,-0.1,
                            -0.1,-0.1,-0.1,1),ncol=4,byrow = F)
colnames(substitutionMatrix)=c('A','T','C','G')
rownames(substitutionMatrix)=c('A','T','C','G')

#strand=topis[1,"Sense"]
align_seq <- function(ischr, isloc, total_is, sam_list, genome_fa) {

  samtop1=(sam_list%>%filter(V3==ischr,V4==isloc)%>%arrange(desc(match_length)))
  if(nrow(samtop1)==0){
    samtop1=(sam_list%>%filter(V3==ischr,reference_end==isloc)%>%arrange(desc(match_length)))
  }
  
  samtop2=samtop1[1,,drop=FALSE]
  type = gsub("[0-9]+", "", samtop2$V6)
  number = as.numeric(unlist(strsplit(samtop2$V6, split = "S|M|I|D")))
  type_value = unlist(strsplit(type, split = ""))
  
  print(c(type, number))
  
  queryseq = switch(
    type,
    "MS" = substr(samtop2$V10, 1, number[1]),
    "SM" = substr(samtop2$V10, number[1] + 1, sum(number)),
    "SMS" = substr(samtop2$V10, number[1] + 1, sum(number[1:2])),
    "MIMS" = {
      queryseq1 = substr(samtop2$V10, 1, number[1])
      queryseq2 = substr(samtop2$V10, sum(number[1:2]) + 1, sum(number[1:3]))
      queryseq = paste0(queryseq1, queryseq2)
    },
    "SMIM" = {
      queryseq1 = substr(samtop2$V10, number[1] + 1, sum(number[1:2]))
      queryseq2 = substr(samtop2$V10, sum(number[1:3]) + 1, sum(number))
      queryseq = paste0(queryseq1, queryseq2)
    },
    "M"=as.character(samtop2$V10),
    "MDM"=as.character(samtop2$V10),
    "MDMS"={
      queryseq=substr(samtop2$V10,1,sum(number[-length(number)]))
    },
    "MIM"=as.character(samtop2$V10)
  )
  
  align1 = pairwiseAlignment(
    pattern = queryseq,
    subject = genome_fa,
    type = "global-local",
    substitutionMatrix = substitutionMatrix,
    gapOpening = 2,
    gapExtension = 1
  )
  
  align2 = pairwiseAlignment(
    pattern = as.character(reverseComplement(DNAString(queryseq))),
    subject = genome_fa,
    type = "global-local",
    substitutionMatrix = substitutionMatrix,
    gapOpening = 2,
    gapExtension = 1
  )
  
  score_list = c(score(align1), score(align2))
  
  score_id = which.max(score_list)
  if(score_id==1){
    sub_aln=subject(align1)
  }else{
    sub_aln=subject(align2)
  }
  
  score_ratio = round(score_list[score_id] / nchar(queryseq), 3)
  seq_all = c(
    as.character(sub_aln),
    start(sub_aln),
    end(sub_aln),
    ifelse(as.numeric(samtop2$V2) == 16, "-", "+"),
    round(score_list[score_id], 3),
    score_ratio,
    nchar(queryseq),
    type
  )
  
  return(list(seq_all, c(samtop2, queryseq)))
}

#**************************先根据位点获取信息
# chr=paste0("chr",topis[num,'Chrom'])
# loc=as.numeric(topis[num,"Loc"])
# strand=topis[num,"Sense"]

seq_all_list = list()
combined_id = list()
for (j in seq_len(num)) {
  outs1 = list()
  chr = paste0("chr", topis[j, 'Chrom'])
  loc = as.numeric(topis[j, "Loc"])
  strand = topis[j, "Sense"]
  genome_fa = as.character(Views(genome[names(genome) == as.character(chr)][[1]], start = loc -
                                   300, end = loc + 300))
  ####检查genomefa中是否有N
  genome_fa=gsub('N','',genome_fa)
  for (i in seq_len(num)) {
    top2ischr = paste0("chr", topis[i, 'Chrom'])
    top2isloc = as.numeric(topis[i, 'Loc'])
    #ischr=top2ischr;isloc=top2isloc;total_is=total_is;sam_list=sam_list;genome_fa=genome_fa
    seq_all_data = align_seq(
      ischr = top2ischr,
      isloc = top2isloc,
      total_is = total_is,
      sam_list = sam_list,
      genome_fa = genome_fa
    ) #ischr,isloc,total_is,sam_list,genome_fa
    name_ = paste0(top2ischr, ":", top2isloc)
    if (!name_ %in% names(seq_all_list)) {
      seq_all_list[[name_]] = seq_all_data[[2]]
    }
    outs1[[i]]= if (i == j) {
      rep(1,8)
    } else{
      seq_all_data[[1]][1:8]
    }
  }
  combined_id[[j]] = outs1
}

score_list=data.frame(matrix(rep(0,length=num*num),ncol=num))
colnames(score_list)=paste('topis',1:num,sep = "")
rownames(score_list)=paste('topis',1:num,sep = "")
for(s in 1:num){
  t1=do.call(rbind,combined_id[[s]])
  score_list[s,]=as.numeric(t1[,6])
}

#******************长矩阵的信息,存储取出来的reads信息等
score_info=c()
for(i in 1:num){
  t1=data.frame(do.call(rbind,combined_id[[i]]))
  t1$from=paste('topis',1:num,sep="")
  rownames(t1)=NULL
  t1$to=paste0('topis',i)
  score_info=rbind(score_info,t1)
}

#***********************判断位点的合并信息*********
new_score_list=score_list
for(i in 1:num){
  for(j in 1:num){
    new_score_list[i,j]=max(score_list[i,j],score_list[j,i])
  }
}
######################保存score值得分
results_score_list=new_score_list
names_id=paste0(paste(paste0("chr",topis[,1],sep=""),topis[,3],sep=":"),":","topis",1:nrow(topis),sep = "")
rownames(results_score_list)=names_id
colnames(results_score_list)=names_id
saveRDS(results_score_list,file.path(outdir,paste0(samples,'_new_score_list.RDS')))
write.table(results_score_list,file.path(outdir,paste0(samples,'_new_score_list.txt')),sep="\t",quote=F,row.names=T,col.names=T)

#*
#*************************判断位点的距离近似***********************

# Add a column with the topic names called "from"
# df <- as.data.frame(new_score_list)
# 
# df <- df %>% gather(from, value, 1:ncol(new_score_list))
# 
# df$to=rep(rownames(new_score_list),num)
# loc_list=as.list(apply(topis[1:num,c(1,3)],1,function(t)paste0('chr',t[1],":",t[2])))
# names(loc_list)=paste0("topis",1:num,sep="")
# 
# df$from=unlist(loc_list[df$from])
# df$to=unlist(loc_list[df$to])

# Add a column with the topic names called "to"
# df <- df %>% separate(from, into = c("tmp", "to"), sep = "_") %>% select(-tmp)

# View the resulting data frame
# 当只有一个TOPIS点的时候，不过滤IS位点
# 当有多个TOPIS点的时候，不满足条件的时候，过滤后可能为空
# 当有多个TOPIS点的时候，有存在满足条件的点，过滤后不为空,即对不为空的IS位点做聚类，剩下的位点单独为一类
# check_filter_id=sum((df$value<=0.8)&(df$value!=1))
loc_list=as.list(apply(topis[1:num,c(1,3)],1,function(t)paste0('chr',t[1],":",t[2])))
names(loc_list)=paste0("topis",1:num,sep="")
# 
colnames(new_score_list)=loc_list[colnames(new_score_list)]
rownames(new_score_list)=loc_list[rownames(new_score_list)]

new_score_list=as.matrix(new_score_list)
g <- graph_from_adjacency_matrix(new_score_list, mode = "undirected", weighted = TRUE, diag = FALSE)

g_filtered <- subgraph.edges(g, E(g)[E(g)$weight > 0.8], delete.vertices = FALSE)
clusters <- components(g_filtered)
# clusters$membership,聚类的位点信息统计号
unique_loc=split(clusters$membership,f=clusters$membership)
unique_loc_names=sapply(unique_loc,function(t)names(t)[1])
unique_loc=sapply(unique_loc,function(t)names(t))
names(unique_loc)=unique_loc_names
unique_is=unique_loc

if(length(unique_is)>=1){
  print(paste0("there is ",length(unique_is),"clusters in the loc!\n"))
}
topis$new=paste0('chr',topis$Chrom,":",topis$Loc)
topis$id=1:nrow(topis)

results=c()
j=1
for(i in 1:nrow(topis)){
  if(topis$new[i]%in%names(unique_is)){
    cluster_tmp=paste0('cluster',j)
    results=rbind(results,data.frame(loc=unique_is[[topis$new[i]]],cluster=cluster_tmp))
    j=j+1
  }
}

results$Chrom=do.call(rbind,strsplit(results$loc,split=":"))[,1]
results$LOCATION=do.call(rbind,strsplit(results$loc,split=":"))[,2]
# clusterid=sapply(strsplit(results$cluster,split="cluster"),function(t)as.numeric(t[2]))
# results$clusterid=clusterid
# results=results%>%group_by(clusterid)%>%arrange(desc())

topis_add_results=merge(results,topis,by.x=c('loc'),by.y=c('new'))
topis_add_results=topis_add_results[order(topis_add_results$id,decreasing = F),]
topis_add_results$Chrom.x=NULL
topis_add_results$LOCATION=NULL
colnames(topis_add_results)[grep('Chrom',colnames(topis_add_results))]="Chrom"
#******************************cluster 信息和sam信息存储
sam_select=data.frame(do.call(rbind,seq_all_list))
sam_select$loc_info=names(seq_all_list)
colnames(sam_select)[16]='queryseq'
rownames(sam_select)=sam_select$loc_info
sam_select$clusterinfo=topis_add_results[match(sam_select$loc_info,topis_add_results$loc),"cluster"]
sam_select=unique(sam_select[topis_add_results$loc,])
sam_select$value=topis_add_results[match(sam_select$loc_info,topis_add_results$loc),9]
sam_select$clusteid=sapply(strsplit(sam_select$clusterinfo,split="cluster"),function(t)as.numeric(t[2]))
sam_select=sam_select[order(sam_select$clusteid,decreasing = F),] ###按照cluster排序
sam_select$clusteid=NULL

saveRDS(sam_select,file.path(outdir,paste0(samples,'_sam_select.RDS')))

write_xlsx(unique(topis_add_results),file.path(outdir,paste0(samples,"_topis_combin_loc.xlsx")))
#*******************************************************

