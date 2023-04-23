# tagtable=read.delim(tagtablefile,sep="\t",as.is=T,header=F,check.names=F)
# samples=sapply(tagtable$V4,function(t){unlist(strsplit(t,split="_"))[5]})
#'tagtable为tagtable文件tagTabMiS026.txt
#'indir为输入文件路径/mount/labdatastorage001/input/konghl/data/plant/MiS026/zea/A001_2
#'samples为样本名称

get_is_sam_data <- function(tagtable,indir,samples){
  
  is_r1_list=paste0(tagtable$V4,'_._20.is.r1')
  is_list=paste0(tagtable$V4,'_._20.is')
  is=c()
  is_r1=c()
  sam_list=c()
  
  for(i in 1:nrow(tagtable)){
    tmp_isfile=file.path(indir,'is_file',paste0(tagtable$V4[i],"_._20.is"))
    if(!file.exists(tmp_isfile)){
      tmp_isfile=file.path(indir,paste0(tagtable$V4[i],"_._20.is"))
    }
    if(file.exists(tmp_isfile)&file.size(tmp_isfile)!=0){
      tmp_is=read.delim(tmp_isfile,sep=" ",as.is=T,header=F)
      if(nrow(tmp_is)>=1){
        tmp_is$sample=samples[i]
      }
      is=rbind(is,tmp_is)  
    }
    tmp_isr1file=file.path(indir,'is_file',paste0(tagtable$V4[i],"_._20.is.r1"))
    if(!file.exists(tmp_isr1file)){
      tmp_isr1file=file.path(indir,paste0(tagtable$V4[i],"_._20.is.r1"))
    }
    if(file.exists(tmp_isr1file)&file.size(tmp_isr1file)!=0){
      tmp_isr1=read.delim(tmp_isr1file,sep=" ",as.is=T,header=F)
      if(nrow(tmp_isr1)>=1){
        tmp_isr1$sample=samples[i]
      }
      is_r1=rbind(is_r1,tmp_isr1)
    }
  }
  
  is$type="is"
  is_r1$type="is_r1"
  total_is=rbind(is,is_r1)
  
  for(i in 1:nrow(tagtable)){
    samfile=file.path(indir,"sam_file",paste0(tagtable$V4[i],'_.fastq.trimmed2.sam'))
    if(!file.exists(samfile)){
      samfile=file.path(indir,paste0(tagtable$V4[i],'_.fastq.trimmed2.sam'))
    }
    if(file.exists(samfile)&file.size(samfile)!=0){
      tmp_sam=read.delim(samfile,sep="\t",as.is=T,header=F)
      if(nrow(tmp_sam)>=1){
        tmp_sam$sample=samples[i]
      }
      sam_list=rbind(sam_list,tmp_sam)
    }
  }
  return(list(total_is))  
}