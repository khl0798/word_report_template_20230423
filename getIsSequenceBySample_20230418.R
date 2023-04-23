####按照每个样本出报告，
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
invisible(capture.output(suppressPackageStartupMessages(library(writexl))))
invisible(capture.output(suppressPackageStartupMessages(library(parallel))))

# source('/mount/labdatastorage001/input/konghl/script/UISRIS_combine/clusterIS_20230329.R')
clusterIS_script="/mount/labdatastorage001/input/konghl/script/UISRIS_combine/clusterISbySample_20230423.R"
getISsequence_script="/mount/labdatastorage001/input/konghl/script/UISRIS_combine/getIStosequencebySample_genomefa_20230423.R"
source('/mount/labdatastorage001/input/konghl/script/plant/get_is_sam_data_20230330.R')
extractFromSam_script="/mount/labdatastorage001/input/konghl/data/plant/MiS027/results_os_new/test/extractFromSamFile.py"
python_script="/home/adminKHL/anaconda3/envs/python3.8/bin/python"

spec <- matrix(
  c("config",  "c", 1, "character", "config",
    "help", "h", 0,"loical", "帮助文档"),
  byrow=TRUE, ncol=5)

opt=getopt(spec=spec)

if (is.null(args$config)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

config=opt$config
config_data=read.delim(config,sep="=",header=F,as.is=T)
###读取输入文件
tagtablefile=config_data[match("tagtablefile",config_data$V1),'V2'] #*** tagtablefile文件路径
tagtable=read.delim(tagtablefile,sep="\t",as.is=T,header=F)
top_is_file=config_data[match("top_is_file",config_data$V1),'V2']  #*** top_is_file文件路径
fq1=config_data[match("fq1",config_data$V1),'V2']  #*** fq1 文件路径
select_is_percent=as.numeric(config_data[match("select_is_percent",config_data$V1),'V2'])  #*** 选择top_is_percent比例信息
readslength=as.numeric(config_data[match("readslength",config_data$V1),'V2'])  #*** read_length信息
indir=config_data[match("indir",config_data$V1),'V2']   #***输入路径 
outdir=config_data[match("outdir",config_data$V1),'V2'] #***输出路径
genome_fa=config_data[match("genome_fa",config_data$V1),'V2']
print(c(tagtablefile,top_is_file,fq1,select_is_percent,readslength,indir,outdir,genome_fa))

####拆分isovertime_sample 文件
topis=data.frame(read_excel(top_is_file))
megapcr_id=which('MegaPCR'==colnames(topis))
total_counts=topis[2,]

samples_counts=topis[2,7:(megapcr_id-1)] ##获取样本信息
samples_IsNum=topis[1,7:(megapcr_id-1)]
samples_uniqueIsNum=topis[3,7:(megapcr_id-1)]

# i=1
cmd_list=c()

for(i in seq_len(length(7:(megapcr_id-1)))){
  print(i)
  #### 每次只取一个样本的信息
  samples1=topis[,c(1:6,6+i,megapcr_id:ncol(topis))]
  
  samples1_new=samples1[4:nrow(samples1),] ##提取每个样本的IS位点信息
  samples1_new=samples1_new[samples1_new[,7]>=select_is_percent,] ##考虑是大于等于percent百分比的信息
  
  samples1_new[,5]=round(as.numeric(samples1_new[,7])*as.numeric(unlist(samples_counts[i]))/100) ## 重新计算每个样本的值
  #samples1_new[,6]=1  ###只有一个样本，待定，是否为1
  samples1_new=rbind(samples1[1:3,],samples1_new)
  each_sample_configfile=file.path(outdir,paste0(names(samples_counts)[i],'.xlsx'))
  write_xlsx(samples1_new,each_sample_configfile)
  ######重新写入配置文件中，每个样本一个配置文件信息
  #new_config=read.delim('config.txt',sep="=")
  ######开始配备每个样本的信息
  new_config_data=config_data
  new_config_data[match('top_is_file',config_data$V1),'V2']=each_sample_configfile
  each_sample_configpath=file.path(outdir,paste0(names(samples_counts)[i],'_config.txt'))
  write.table(new_config_data,each_sample_configpath,sep="=",quote=F,row.names=F,col.names=F)
  ######修正sam文件，添加from 和end 信息
  cmd=paste(python_script,extractFromSam_script,"-t",tagtablefile,"-o",outdir,"-i",file.path(indir,"sam_file"),"-s",unlist(strsplit(names(samples_counts)[i],split="x"))[2])
  print(cmd)
  system(cmd)
  ######重新生成每个样本的配置文件
  cmd=paste("Rscript",getISsequence_script,"-c",each_sample_configpath)
  print(cmd)
  cmd_list=c(cmd_list,cmd)
  system(cmd)
}
####并行运行计算
cl <- makeCluster(length(cmd_list))
run_cmd<-function(cmd){system(cmd)}
result <- parLapply(cl, cmd_list, run_cmd)
stopCluster(cl)

######输出word信息
#####获取序列位点前后的2000bp的序列信息
###从sample_loc获取信息
genome = readDNAStringSet(genome_fa, "fasta")
sample_loc=readRDS(paste0(unlist(strsplit(names(samples_counts)[i],split="x"))[2],'_sam_select.RDS'))
library(Biostrings)
sample_loc=sample_loc[!duplicated(sample_loc$clusterinfo),]
sequence_list=list()
j=1
#genome[names(genome)==as.character(id1[1,3])][[1]]
loc_info=t(sapply(strsplit(sample_loc$loc_info,split=":"),function(t)t))
x1=as.character(Views(genome[names(genome)==as.character(loc_info[j,1])][[1]], start=as.numeric(loc_info[j,2])-2000, 
      end=as.numeric(loc_info[j,2])+2000))
if(sample_loc$V2[j]==16){
  x1=as.character(reverseComplement(DNAStringSet(x1)))
}
sequence_list[[sample_loc$loc_info[j]]]=c(sample_loc$loc_info[j],
                                          ifelse(sample_loc$V2[j]==16,'-','+'), x1)
word_path = file.path(outdir, paste0('LTA_zea_',sample_name,'_topISfastq.docx'))

my_doc <- read_docx()
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
my_doc <- my_doc %>% body_add_fpar(fpar(ftext("topISinfo:",properties5)))
my_doc <- my_doc %>% body_add_fpar(fpar(ftext("\n",properties5)))


###在运行脚本前
# library(lintr)
# lint('D://项目/工作/脚本/getIsSequenceBySample_20230418.R')
# samplesID=c('0011','0012','0021','0022','0031','0032')
# path="/mount/labdatastorage001/input/konghl/data/plant/MiS027"
# uis_list=c()
# ris_list=c()
# for(s in samplesID){
#   tmp_sam=paste0("results_os_new_",s)
#   uis=file.path(path,tmp_sam,"MiS027_os.25-0.9-0.9-150-1000.ResultsClusteredAnnotated.csv")
#   ris=file.path(path,tmp_sam,"MiS027_os.25-0.9-0.9-150-1000.repeats.ResultsClusteredAnnotated.csv")
#   uis_data=read.csv(uis,as.is=T,header=F)
#   ris_data=read.csv(ris,as.is=T,header=F)
#   if(s!=samplesID[1]){
#      uis_data=uis_data[-1,]
#      ris_data=ris_data[-1,]
#   }
#   uis_list=rbind(uis_list,uis_data)
#   ris_list=rbind(ris_list,ris_data)
# }
# write.table(uis_list,file.path(path,'MiS027_os.25-0.9-0.9-150-1000.ResultsClusteredAnnotated.csv'),sep=",",
#             quote=F,row.names=F,col.names=F)
# write.table(ris_list,file.path(path,'MiS027_os.25-0.9-0.9-150-1000.repeats.ResultsClusteredAnnotated.csv'),sep=",",
#             quote=F,row.names=F,col.names=F)

