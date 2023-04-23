.libPaths('/home/nishuai/R/x86_64-pc-linux-gnu-library/4.2/')
invisible(capture.output(suppressPackageStartupMessages(library(getopt))))
invisible(capture.output(suppressPackageStartupMessages(library(readxl))))
invisible(capture.output(suppressPackageStartupMessages(library(data.table))))
invisible(capture.output(suppressPackageStartupMessages(library(dplyr))))
invisible(capture.output(suppressPackageStartupMessages(library(Biostrings))))
invisible(capture.output(suppressPackageStartupMessages(library(stringr))))
invisible(capture.output(suppressPackageStartupMessages(library(magrittr))))
invisible(capture.output(suppressPackageStartupMessages(library(officer))))
invisible(capture.output(suppressPackageStartupMessages(library(ShortRead))))
invisible(capture.output(suppressPackageStartupMessages(library(flextable))))
# source('/mount/labdatastorage001/input/konghl/script/UISRIS_combine/clusterIS_20230329.R')
clusterIS_script="/mount/labdatastorage001/input/konghl/script/UISRIS_combine/clusterISbySample_20230423.R"
source('/mount/labdatastorage001/input/konghl/script/plant/get_is_sam_data_20230330.R')
seqtk_path="/usr/bin/seqtk"

spec <- matrix(
  c("config",  "c", 1, "character", "config"),
  byrow=TRUE, ncol=5)

opt=getopt(spec=spec)
config=opt$config
config_data=read.delim(config,sep="=",header=F,as.is=T)

#****************************************测试数据信息
# tagtablefile="/mount/labdatastorage001/input/konghl/data/plant/MiS026/zea/A001_1/tagTabMiS026.txt"
# tagtable=read.delim(tagtablefile,sep="\t",as.is=T,header=F)
# top_is_file="/mount/labdatastorage001/input/konghl/data/plant/MiS026/zea/A001_1/results_A001_1/ISOT_results_MiS026_A001_1/isovertime_results_by_replicate_MiS026_A001_1.xlsx"
# fq1="/mount/labdatastorage001/input/fastq/2023/MiS026fastq/MiS026_S1_L001_R1_001.fastq.gz"
# select_is_percent=0.1
# readslength=301
# indir="/mount/labdatastorage001/input/konghl/data/plant/MiS026/zea/A001_1/results_A001_1"
# outdir=file.path(getwd())
#***********************read argument from config_data
tagtablefile=config_data[match("tagtablefile",config_data$V1),'V2'] #*** tagtablefile文件路径
tagtable=read.delim(tagtablefile,sep="\t",as.is=T,header=F)
top_is_file=config_data[match("top_is_file",config_data$V1),'V2']  #*** top_is_file文件路径，即isovertime文件信息
fq1=config_data[match("fq1",config_data$V1),'V2']  #*** fq1 文件路径
select_is_percent=as.numeric(config_data[match("select_is_percent",config_data$V1),'V2'])  #*** 选择top_is_percent比例信息
readslength=as.numeric(config_data[match("readslength",config_data$V1),'V2'])  #*** read_length信息
indir=config_data[match("indir",config_data$V1),'V2']   #***输入路径 
outdir=config_data[match("outdir",config_data$V1),'V2'] #***输出路径
print(c(tagtablefile,top_is_file,fq1,select_is_percent,readslength,indir,outdir))

t1=Sys.time()

topis=read_excel(top_is_file)
topis=topis[!is.na(topis$Loc),]
megapcr_id=which('MegaPCR'==colnames(topis))
samples=sapply(strsplit(colnames(topis)[7],split="x"),function(t)t[2])
sample_name=samples[1]

tagtable=tagtable[grep(paste0(sample_name,"x"),tagtable$V4),]

is_sam=get_is_sam_data(tagtable=tagtable,indir=outdir,samples=samples)
is_data=is_sam[[1]]
sam_data=is_sam[[2]]

value=data.frame(do.call(cbind,lapply(data.frame(topis[,7:(megapcr_id-1)],check.names = F,stringsAsFactors = F),as.numeric)))
select_is=which(rowMeans(value)>=select_is_percent)
topis_select=data.frame(topis[select_is,],check.names=F)

# sam_data=read.delim(sam_file,sep="\t",as.is=T,header=F)

########对第一个位点进行信息提取，选择match长度最长的reads信息
# get_sam_info= as.data.frame(t(sapply(get_sam_info, unlist)), stringsAsFactors = FALSE)

get_matchseq_info<-function(align_type,get_sam_info){
  if (substr(align_type,1,1)=='M' & substr(align_type,nchar(align_type),nchar(align_type))=="S") {
    if(align_type=="MIMS"|align_type=="MDMS"){
      value = as.numeric(unlist(strsplit(get_sam_info[1, 6], split = "S|M|I|H|D")))
      clip_loc=which(unlist(strsplit(align_type,split=""))=="S")
      matchseq1 = substr(get_sam_info[1, 10], 1,value[1])
      matchseq2 = substr(get_sam_info[1, 10], sum(value[1:2]),sum(value[1:2]))
      matchseq3 = substr(get_sam_info[1, 10], sum(value[1:2])+1,sum(value[1:3]))
      matchseq=paste0(matchseq1,matchseq2,matchseq3)
      matchcseq_list=c(align_type,matchseq1,matchseq2,matchseq3)
    }else{
      value = as.numeric(unlist(strsplit(get_sam_info[1,6], split = "S|M|I|H")))
      clip_loc=which(unlist(strsplit(align_type,split=""))=="S")
      matchseq = substr(get_sam_info[1, 10],1, value[1])
      matchcseq_list=c(align_type,matchseq)
    }
  }
  if (align_type == "SMS") {
    value = as.numeric(unlist(strsplit(get_sam_info[1,6], split = "S|M|I|H|D")))
    matchseq = substr(get_sam_info[1, 10], value[1] + 1, sum(value[1:2]))
    matchcseq_list=c(align_type,matchseq)
  }
  if (substr(align_type,1,1)=='S' & substr(align_type,nchar(align_type),nchar(align_type))=="M") {
    if(align_type=="SMDM"|align_type=="SMIM"|align_type=="SDM"){
      value = as.numeric(unlist(strsplit(get_sam_info[1, 6], split = "S|M|I|H|D")))
      clip_loc=which(unlist(strsplit(align_type,split=""))=="S")
      match_search=which(unlist(strsplit(align_type,split=""))=="M")
      match_id=which.max(value[which(unlist(strsplit(align_type,split=""))=="M")])
      match_max=value[match_search[match_id]]
      matchseq_id=which(value==match_max)
      matchseq=substr(get_sam_info[1,10],sum(value[matchseq_id-1])+1,sum(value[1:matchseq_id]))
      matchcseq_list=c(align_type,matchseq)
    }else{
      value = as.numeric(unlist(strsplit(get_sam_info[1,6], split = "S|M|I|H")))
      matchseq = substr(get_sam_info[1, 10], value[1] + 1, sum(value[1:2]))
      matchcseq_list=c(align_type,matchseq)
    }
  }
  if(align_type=="M"|align_type=="MDM"){
    matchseq=get_sam_info[1,10]
    matchcseq_list=c('M',matchseq)
  }
  return(list(matchseq,matchcseq_list))
}

megaprimer_list=unlist(sapply(tagtable$V4,function(t)unlist(strsplit(t,split="_"))[18]))
lkbc_list=unlist(sapply(tagtable$V4,function(t)unlist(strsplit(t,split="_"))[23]))
names(megaprimer_list)=paste0(samples,'xR',1:length(megaprimer_list))
names(lkbc_list)=paste0(samples,'xR',1:length(megaprimer_list))
######## 生成sam文件信息
samfile=file.path(outdir,paste0(sample_name,'.test1.sam'))
cmd=paste0("Rscript ",clusterIS_script," -c ",config," -s ",samfile)
print(cmd)
system(cmd)

if(file.exists(file.path(outdir,paste0(sample_name,'_sam_select.RDS')))){
  sam_select=readRDS(file.path(outdir,paste0(sample_name,'_sam_select.RDS')))
  sam_select_filter=sam_select[!duplicated(sam_select$clusterinfo),]
}

readid=gsub('HWI','M0',as.character(sam_select_filter[,1]))
readidfile=file.path(outdir,paste0(sample_name,'readid.txt'))
write.table(readid,readidfile,sep="\t",quote=F,row.names=F,col.names=F)
readidfastqfile=file.path(outdir,paste0(sample_name,'read.fastq'))
fq1_cmd = paste0("zcat ", fq1, "| ",seqtk_path," subseq - ", readidfile, " > ", readidfastqfile)
print(fq1_cmd)
system(fq1_cmd)
read_fastq<-function(readidfastqfile){
  data=read.delim(readidfastqfile,sep="\t",as.is=T,header=F)
  readidlist=split(data,f=rep(1:(nrow(data)/4),each=4))
  names(readidlist)=sapply(strsplit(gsub('@M0','HWI',data[seq(1,nrow(data),4),1]),split=" "),function(t)t[1])
  return(readidlist)
}

fastq_list=read_fastq(readidfastqfile)

loc_list = list()

# for (i in 1:nrow(topis_select)) {
for(i in 1:nrow(sam_select_filter)){
  chr = unlist(strsplit(sam_select_filter[i, 'loc_info'],split=":"))[1]
  loc = unlist(strsplit(sam_select_filter[i, 'loc_info'],split=":"))[2]
  strand=ifelse(sam_select_filter[i,'V2']==0,"+","-") ##不知道这里是否会有小bug，先这样子做吧
  # strand = topis_select[i, 2]
  unknown_1="" ## 比对到基因组前面的未知序列
  unknown_2="" ##比对到基因组后面的未知序列
  unknown_3="" ##lkbc后面的未知序列
  # is_r1_filter = data.frame(is_data %>% filter(V2 == chr, V3 == loc)) #根据染色体和loc筛选reads信息
  is_r1_filter=sam_select_filter[i, ,drop=FALSE]
  #****************计算megaprimer和lkbc的信息
  
  print(nrow(is_r1_filter))
  if(nrow(is_r1_filter)>=1){
    
    readid=is_r1_filter[,1]
    
    align_type = gsub('[0-9]', '', is_r1_filter[, 6])#获取比对信息
    align_type=align_type
    print(c(i,align_type))
    
    # is_r1_filter=as.data.frame(t(sapply(is_r1_filter, unlist)), stringsAsFactors = FALSE)
    
    matchseqlist=list()
    matchseqlist[[1]]=get_matchseq_info(align_type,is_r1_filter)
    matchseq=do.call(rbind,lapply(matchseqlist,function(t)t[[1]]))
    matchcseq_list=do.call(rbind,lapply(matchseqlist,function(t)t[[2]]))
    print(c(matchseq,matchcseq_list))
    matchseq = sapply(1:length(matchseq),function(t)ifelse(is_r1_filter$V4[t]==16,as.character(reverseComplement(DNAString(matchseq[t]))),matchseq[t]))
    
    
    if(unlist(readid)%in%names(fastq_list)){
      fq1_data=data.frame(fastq_list[[unlist(readid)]])
    }else{
      outfile = file.path(outdir,paste0(chr, "_", loc, '_', 'top', i))
      outreadidfile=file.path(outdir,paste0(chr,"_",loc,"_readsid"))
      write.table(gsub("HWI","M0",readid),outreadidfile,sep="\t",quote=F,row.names=F,col.names=F)
      
      fq1_cmd = paste0("zcat ", fq1, "| ",seqtk_path," subseq - ", outreadidfile, " > ", outfile)
      
      if (!file.exists(outfile) | file.size(outfile) == 0) {
        system(fq1_cmd)
      }else{
        file.remove(outfile)
        system(fq1_cmd)
      }
      
      fq1_data = read.delim(outfile,
                            sep = '\t',
                            as.is = T,
                            header = F)
      if(file.exists(outreadidfile)){
        file.remove(outreadidfile)
      }
    }
    lkbc=lkbc_list[as.character(is_r1_filter[1,"sample"])]
    megaprimer=megaprimer_list[as.character(is_r1_filter[1,"sample"])]
    
    readslength=nchar(fq1_data[2,1])
    mismatch=max(3,round(ceiling(as.numeric(is_r1_filter[1,'queryseq'])/10)))
    
    if(is_r1_filter[1,2]==16){
      ref_align_seq=as.character(reverseComplement(DNAString(as.character(is_r1_filter[1,'V18']))))
    }else{
      ref_align_seq=as.character(is_r1_filter[1,'V18']) ###序列信息V18列
    }
    hits = matchPattern(pattern =ref_align_seq,fq1_data[2, 1],max.mismatch = mismatch,with.indels = T)
    
    ref_seq_id=unlist(as.data.frame(ranges(hits)))[1:2]
    
    megaprimer_id = unlist(str_locate_all(fq1_data[2, 1], pattern = megaprimer))
    if(length(megaprimer_id)==0){
      hits = matchPattern(pattern = megaprimer,fq1_data[2, 1],max.mismatch = 2)
      megaprimer_id=unlist(data.frame(ranges(hits))[1,1:2])
    }
    lkbc_id=unlist(str_locate_all(fq1_data[2, 1], pattern = lkbc))
    if(length(lkbc_id)==0){
      hits = countPattern(lkbc,fq1_data[2, 1],max.mismatch = 2)
      if(hits==1){
        hits = matchPattern(pattern = lkbc,fq1_data[2, 1],max.mismatch = 2)
        lkbc_id=unlist(data.frame(ranges(hits))[1,1:2])
      }else{
        lkbc_id=unlist(str_locate_all(fq1_data[2, 1], pattern = 'CCTAAC'))
        if(length(lkbc_id)==0){
          if(readslength-ref_seq_id[2]<=10){
            lkbc_id = unlist(str_locate_all(fq1_data[2, 1], pattern = substr(lkbc, 1, ifelse(readslength-ref_seq_id[2]>=5,4,readslength-ref_seq_id[2]))))
            if(length(lkbc_id)>=10){
              lkbc_id=c(ref_seq_id[2],nchar(fq1_data[2,1]))
            }
          }else{
            lkbc_id=c(ref_seq_id[2]+1,nchar(fq1_data[2,1]))
          }
        }
      }
    }
    
    barcodes = substr(fq1_data[2, 1], 1, megaprimer_id[1] - 1)
    
    if (megaprimer_id[2] + 1 !=ref_seq_id[1]) {
      unknown_1 = substr(fq1_data[2, 1], megaprimer_id[2] + 1, ref_seq_id[1] -1)
    }
    if(length(lkbc_id)==0){
      unknown_2=substr(fq1_data[2, 1], ref_seq_id[2] + 1, nchar(fq1_data[2, 1]))
      unknown_3=""
      tail_seq=""
    }else{
      if ((ref_seq_id[2] + 1) != lkbc_id[1]) {
        unknown_2 = substr(fq1_data[2, 1], ref_seq_id[2] + 1, lkbc_id[1] - 1)
      }
      if(lkbc_id[2]<readslength){
        unknown_3=substr(fq1_data[2,1],lkbc_id[2]+1,nchar(fq1_data[2,1]))
      }
      tail_seq = substr(fq1_data[2, 1], lkbc_id[1], lkbc_id[2])
    }
    
    loc_list[[paste(chr, loc, sep = "_")]] = list(c(barcodes,
                                                    as.character(megaprimer),
                                                    as.character(unknown_1),
                                                    as.character(matchseq),
                                                    as.character(unknown_2),
                                                    as.character(tail_seq),as.character(unknown_3)),
                                                  is_r1_filter,list(unlist(matchcseq_list[1,])))
    ### 删除 outreadidfile 文件

  }
}
saveRDS(loc_list,file.path(outdir,paste0(sample_name,"_loc_list.RDS")))


copy_numbers=nrow(sam_select_filter)
#步骤四：写入内容


fp_text_properties <- function(color="black",bold=TRUE,font.size=14,shading.color=NA){
  properties <- fp_text(color = color, #normal
                        bold = bold,
                        font.size = font.size,shading.color=shading.color)
  return(properties)
}
properties5=fp_text_properties() #black normal
properties1=fp_text_properties(color = "#00F5FF") #blue megaprimer
properties2=fp_text_properties(color = "#778899") #grey unknown
properties3=fp_text_properties(color = "#FF1493") #green genome
properties4=fp_text_properties(color = "#CDCD00") #red bclinker
properties6=fp_text_properties(color = "red",bold =TRUE) #red bclinker


colorlist_total=c("#00CED1","#EE82EE","#BE2A3E", "#EC754A",
                  "#0780cf", "#765005", "#fa6d1d", "#0e2c82", "#b6b51f", "#da1f18",
                  "#701866", "#f47a75", "#009db2", "#024b51","#0780cf", "#765005",
                  "#63b2ee", "#76da91", "#f8cb7f", "#f89588", "#7cd6cf", "#9192ab", "#7898e1" ,"#efa666",
                  "#eddd86", "#9987ce", "#63b2ee","#76da91")
color_list=c('#F26659','#596235','orange',colorlist_total)
# refA_color="#EE859A"
# tdna_color="#65A9DD"
# insert_lb_color="#E5C4E4"
# insert_rb_color="#C8DE88"
# refB_color="#FFB90F"
std_border = fp_border(color="black", width = 2)
std_border_inner=fp_border(color="black", width = 1)


word_path = file.path(outdir, paste0('LTA_zea_',sample_name,'_topISfastq.docx'))

my_doc <- read_docx()

#固定格式模板信息
title1=paste0("样本整合位点分析结果及序列信息展示")
info1="检测位点信息说明:"
# info2="实验结果：以下列出在任意复孔中克隆占比大于1%的t-DNA整合位点，在克隆占比大于5%时，本方法可以精确定量。在克隆占比在1%-5%时，本方法可以稳定检出。"
# info3="玉米基因组中存在大量重复序列，其中包含一个156 bp的着丝粒重复序列，一个9349 bp的45S rDNA重复序列和一个341 bp的5S rDNA重复序列。此外，玉米还包含两种丰富的knob重复序列，分别为主要存在的knob重复序列（~180 bp）和TR-1重复序列（~360 bp）；总体重复序列占比大约为玉米基因组的85%。由于本次实验的技术条件限制，如果t-DNA恰好整合在玉米基因组的重复区域内，我们无法确定t-DNA整合在基因组上的具体位置。因此我们对于该类t-DNA整合事件输出所有可能的整合位点以供验证。"
#********************以下文字的结果很不智能
info4="本次样本中共检测到"
info5=copy_numbers
info6=paste0("个整合事件，其中列表top","1-",nrow(topis_select),"位点列出了该整合事件在基因组上基于1%克隆占比检测线检出的最有可能的整合位点,属于同一个整合事件的整合位点用相同的颜色表示。")

#。附带文件“002样本整合事件blat比对结果.bst”为该整合事件所有可能的整合位点以供验证。"

my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0(title1),fp_text_properties(font.size=18)),fp_p=fp_par(text.align="center")))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))
# my_doc <- my_doc %>% body_add_fpar(fpar(ftext("\n",properties5)))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0(info1),properties5)))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))
# my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0(info2),properties5)))
# my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))
# my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0(info3),properties5)))
# my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0(info4),properties5),
                                        ftext(paste0(info5),properties6),
                                        ftext(paste0(info6,"\n"),properties5)))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))

info1=ftext('Fastq颜色说明: ',properties5)

info2=ftext(' 未知序列',properties2)
info3=ftext(' megaprimer',properties1)
info4=ftext(' 比对到基因组序列',properties3)
info5=ftext(' bclinker\n',properties4)
info6=ftext('\n',properties4)

my_doc <- my_doc %>% body_add_fpar(fpar(info1,info2,info3,info4,info5,info6))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))

my_doc <- my_doc %>% body_add_fpar(fpar(ftext("topISinfo:",properties5)))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext("\n",properties5)))

check_id=which(colnames(topis_select)=='MegaPCR')
tabledata=topis_select[,c(1:(check_id-1)),drop=FALSE]
rownames(tabledata)=paste0('chr',tabledata[,1],":",tabledata[,3])
###判断样本数目，展示6列样本信息

tabledata=data.frame(cluster=sam_select[rownames(tabledata),'clusterinfo'],tabledata)
tabledata$newid=as.numeric(gsub("cluster","",tabledata$cluster))
tabledata=tabledata[order(tabledata$newid,decreasing=F),]
tabledata$newid=NULL

###对tabledata从第8列至最后一列数据，字符串转数字类型
for(m in 8:ncol(tabledata)){
  tabledata[,m]=as.numeric(tabledata[,m])
}
####判断样本数目，如果样本数目大于6个，则需要分列进行展示
tablelist=list()
if(ncol(tabledata)>13){ #7+6
  split_col_id=split(c(8:ncol(tabledata)),f=rep(1:6,each=6))
  split_col=sapply(split_col_id,length)
  split_col_id_select=split_col_id[split_col>0]
  #print(split_col_id_select)
  tablelist=lapply(split_col_id_select,function(t)tabledata[,c(1:7,unlist(t))])
}else{
  tablelist[[1]]=tabledata
}
####选择不同的table输出格式信息


####
for(m in 1:length(tablelist)){
  table_split=tablelist[[m]]
  ft_2 <- flextable(table_split)
  ft_2 <- border_remove(x = ft_2)
  ft_2 <- hline_top(ft_2, part="header", border = std_border) ##添加标题行线
  ft_2 <- hline_bottom(ft_2, part="body", border = std_border) ##添加最后面的行
  ft_2 <- hline(ft_2, part="header", border = std_border) ##添加标题行线
  ft_2=border_inner_h(ft_2, border = std_border_inner, part = "body") ###添加内部的行线
  # ft_2=autofit(ft_2,  part = c("body", "header")) ##自动调整行 add_w = 0.1, add_h = 0.1,
  ft_2=width(ft_2, j = NULL,width=0.8) ##设置宽度, width=0.5,改为默认
  ft_2=height(ft_2, i = NULL, height=0.8) ##设置宽度，height=0.5，改为默认
  
  
  ft_2 <- color(ft_2, color = "black", part = "header") ##标题字体颜色
  ft_2 <- bold(ft_2, part = "all",bold = TRUE) ###所有字体加粗
  ft_2 <- fontsize(ft_2, part = "all",size = 9) ###所有字体加粗
  
  
  for(i in 1:length(unique(table_split$cluster))){
    id=which(tabledata$cluster==unique(table_split$cluster)[i])
    ft_2 <- color(ft_2, i = id, color = color_list[i])  ###设置不同cluster颜色，即不同拷贝的颜色
  }
  
  ft_2 <- set_formatter_type(ft_2)
  my_doc <- body_add_flextable(my_doc, ft_2)
  my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))
  my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0("\n"),properties5)))
  
}
# print(my_doc, target = word_path)

# my_doc <- body_add_table(my_doc, topis_select[,c(1:9,14),drop=FALSE], style = "table_template")

my_doc <- my_doc %>% body_add_fpar(fpar(ftext("\n",properties5)))
#### 输出cluster 按照顺序来输出

for(s in 1:nrow(sam_select_filter)) {
  # chr=paste0('chr',sam_select_filter[s,1])
  # loc=sam_select_filter[s,3]
  # name=paste(chr, loc, sep = "_")
  name=gsub(':','_',sam_select_filter$loc_info[s])
  chr=unlist(strsplit(name,split="_"))[1]
  loc=unlist(strsplit(name,split="_"))[2]
  strand=ifelse(sam_select_filter$V2[s]==16,"-","+")
  
  if(name%in%names(loc_list)){
    fq_data=loc_list[[name]]
    clusterinfo=tabledata[gsub('_',":",name),'cluster'] ##是提取cluster位点信息
    my_doc <- my_doc %>% body_add_fpar(fpar(ftext(
      paste0(clusterinfo, ": ", chr, "\t", loc,"\t", strand), properties5
    )))
    my_doc <- my_doc %>% body_add_fpar(fpar(ftext("\n", properties5)))
    ##reads id 信息
    my_doc <- my_doc %>% body_add_fpar(fpar(ftext(paste0(fq_data[[2]][1,1], '\n'), properties5)))
    
    ftext1 <- ftext(fq_data[[1]][2], properties1)#megaprimer
    ftext2 <- ftext(fq_data[[1]][3], properties2)#unknown_1
    if(fq_data[[3]][[1]][1]=="MIMS"){
      ftext3_1=ftext(fq_data[[3]][[1]][2], properties3)
      ftext3_2=ftext(fq_data[[3]][[1]][3], properties5)
      ftext3_3=ftext(fq_data[[3]][[1]][4], properties3)
    }else{
      ftext3 <-ftext(fq_data[[1]][4], properties3)#matchseq
    }
    ftext4 <- ftext(fq_data[[1]][6], properties4)#tail seq
    ftext5 <- ftext(fq_data[[1]][1], properties5)#barcodes
    ftext6 <- ftext(fq_data[[1]][5], properties2)#unknown_2
    ftext7 <- ftext(fq_data[[1]][7], properties2)#unknown_2
    
    # print(c(fq_data[[1]][1],fq_data[[1]][2],fq_data[[1]][3],fq_data[[1]][4],fq_data[[1]][5],fq_data[[1]][6],fq_data[[1]][7]))
    
    if(fq_data[[3]][[1]][1]!="MIMS"){
      paragraph <- fpar(ftext5, ftext1, ftext2, ftext3, ftext6, ftext4,ftext7)
    }else{
      paragraph <- fpar(ftext5, ftext1, ftext2, ftext3_1,ftext3_2,ftext3_3, ftext6, ftext4,ftext7)
    }
    
    my_doc <- my_doc %>% body_add_fpar(paragraph)
    my_doc <-
      my_doc %>% body_add_fpar(fpar(ftext(paste0('+\n'), properties5)))
    my_doc <-
      my_doc %>% body_add_fpar(fpar(ftext(paste0(fq1_data[4, 1], '\n'), properties5)))
    my_doc <- my_doc %>% body_add_fpar(fpar(ftext("\n", properties5)))
  }
}

print(my_doc, target = word_path)

t2=Sys.time()
print(paste0('consume: ',round(t2-t1,4),' s\n'))
