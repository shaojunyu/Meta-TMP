#################################################################
# Function: Multivariate statistical analysis based on distance matrix
# Call: Rscript PM_Bdiversity.R -m map_file -d dist_file -o output
# R packages used: reshape,ggplot2,pheatmap,pROC,combinat,plyr,vegan,optparse
# Last update: 2015-07-06, Shi Huang, Xiaoquan Su, Gongchao jing
#################################################################

## install necessary libraries
p <- c("reshape","ggplot2","pheatmap","pROC","combinat","plyr","vegan","optparse")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
sourcedir <- Sys.getenv("ParallelMETA")
source(sprintf('%s/Rscript/util.R',sourcedir))
#sourcedir <- Sys.getenv("HOME")
#source(sprintf('%s/util.R',sourcedir))
# make option list and parse command line
option_list <- list(
    make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix table [Required]"),
    make_option(c("-m", "--meta_data"), type="character", help="Input metadata mapping file [Required]"),
    make_option(c("-o", "--outdir"), type="character", default='Beta_diversity', help="Output directory [default %default]"),
    make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]"),
    make_option(c("-n", "--dist_name"), type="character", default='Default', help="The distance metrics name such as Meta-Storms, Jensen-Shannon, Euclidean et al. [Optional, default %default]")    
    )
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$meta_data)) stop('Please supply a meta-data file.')
if(is.null(opts$dist_file)) stop('Please supply an distance matrix table.')

# create output directory if needed
#if(opts$outdir != ".") 
dir.create(outpath1<-paste(opts$outdir,"/",sep=""),showWarnings=FALSE, recursive=TRUE)

filename<-opts$dist_file                       
metadata.filename<-opts$meta_data               
dm_name<-opts$dist_name
prefix_name<-opts$prefix                         
#all_group<-strsplit(opts$category,",")[[1]]  #c("Habitat","Status","Timepoint","HostID","Gender")

con <- file(paste(opts$outdir,'/',prefix_name,'.','Beta_diversity_Values.txt',sep=''))

sink(con, append=TRUE)
sink(con, append=TRUE, type='message')

#--------------------------------
dm<-read.table(filename,header=T,row.names=1)
dm<-dm[order(rownames(dm)),order(colnames(dm))]
allmetadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1)
if(length(allmetadata)==1){metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
                           all_group<-colnames(metadata)<-colnames(allmetadata)
                           }else{
                           allmetadata<-allmetadata[order(rownames(allmetadata)),]
                           metadata<-allmetadata[,which(lapply(allmetadata,var)!=0)]
                           metadata<-metadata[sapply(metadata,class)=="factor"][order(rownames(metadata)),]
                           all_group<-colnames(metadata)
                           }
cat("The sample groupings: ",all_group, "\n",sep=" ")
#--------------------------------Data Check
if(any((colnames(dm)==rownames(dm))==FALSE)) 
  {cat("The column names do not exactly match the row names! Please revise!")}
if(any((rownames(metadata)==rownames(dm))==FALSE)) 
  {cat("The row names in Map file do not exactly match the row names in the distance matrix! Please revise!\n")}

#--------------------------------
dm_v<-matrix(NA,ncol=length(all_group))
#--------------------------------
suppressWarnings(
for(group in all_group) {
    #--------------------------------
    dir.create(outpath2<-paste(outpath1,dm_name,".",group,"/",sep=""))
    if(var(allmetadata[,group])!=0){
    d<-DistBoxplot(dm,dm_name=dm_name,group=metadata[,group],group_name=group,outpath=outpath2)
    unlink(outpath2, recursive=TRUE)
    #--------------------------------
    # print(paste(" All_between ",group," VS All_within ",group," P value (T-test)=", d$p_t,sep=""))
    # print(paste(" All_between ",group," VS All_within ",group," P value (Wilcox-test)=", d$p_w,sep=""))
    #--------------------------------
    ano<-anosim(dm,metadata[,group])
    ano.P<-ano$signif
    ano.R<-ano$statistic
    cat(paste(group,": \n",sep=""))
    print(paste("P value (ANOSIM)=", ano.P,sep=""))
    print(paste("R value (ANOSIM)=", ano.R,sep=""))
    #--------------------------------
    ado<-adonis(dm~metadata[,group])
    ado.P<-ado$aov.tab$P[1]
    ado.F<-ado$aov.tab$F.Model[1]
    cat(paste(group,": \n",sep=""))
    ado$aov.tab
    print(paste("P value (ADONIS/PERMANOVA)=", ado.P,sep=""))
    print(paste("F value (ADONIS/PERMANOVA)=", ado.F,sep=""))
    cat("\n")
    #--------------------------------
    if(is.na(dm_v[1])){
                       dm_v<-d$dm_all_values}else{
                       dm_v<-Cbind(dm_v,d$dm_all_values)}
    }
}
)
colnames(dm_v)<-paste(rep(all_group,each=2),colnames(dm_v),sep=".")
#--------------------------------
suppressMessages(dm_v.melt<-melt(data.frame(dm_v))) ; dm_v.melt<-data.frame(dm_v.melt[which(!is.na(dm_v.melt$value)),])

dm_v.melt<-data.frame(cbind(do.call(rbind,strsplit(as.character(dm_v.melt$variable),".",fixed=TRUE)),dm_v.melt))
colnames(dm_v.melt)<-c("Group","DistType","DataName","Distance")
p<-qplot(x=Group,y=Distance,  data=dm_v.melt, geom="boxplot", fill=DistType, position="dodge",main="", ylab=paste(dm_name," Distance",sep=""))+
        coord_flip()+
        theme_bw()
all_groups<-paste(all_group,collapse="-")
suppressMessages(ggsave(filename=paste(outpath1,"/",prefix_name,".Beta_diversity_Boxplot",".pdf",sep=""),plot=p, limitsize=TRUE, width=6, height=length(all_group)))


sink()
sink(type='message')
