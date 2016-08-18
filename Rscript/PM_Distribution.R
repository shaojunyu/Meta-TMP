#################################################################
# Function: Abundance distribution analysis
# Call: Rscript PM_Distribution.R -i abund_file -o outfile
# R packages used: optparse, reshape2, ggplot2,RColorBrewer, grDevices
# Last update: 2016-05-26, Zheng Sun, Xiaoquan Su
#################################################################
# install necessary libraries
p <- c("optparse","reshape2","ggplot2","RColorBrewer")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
## clean R environment
rm(list = ls())
setwd('./')
## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--abund_file"), type="character", help="Input feature table with Relative Abundance (*.Abd) [Required]"),
  make_option(c("-o", "--outfile"), type="character", default='distribution.pdf', help="Output distribution file [default %default]"),
  make_option(c("-v", "--threshold"), type="double", default=0.01, help="Average value threshold [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$abund_file)) stop('Please supply an abundance table.')
# load data
matrixfile <- opts$abund_file
ave_t <- opts$threshold
#------------------------------------------------------------------------------------
disbar <- t(read.table(matrixfile,header = T, row.names = 1,sep="\t"))
disbar <- disbar[names(sort(rowSums(disbar),decreasing = T)),]
disbar <- floor(disbar*1000000)/1000000
data_matrix_other <- disbar[which(apply(disbar,1,mean) <= ave_t),]

invisible(if (sum(data_matrix_other) ==0 ) (data_matrix_big <- disbar))
invisible(if (sum(data_matrix_other) !=0 ) (data_matrix_big <- disbar[-(which(apply(disbar,1,mean) <= ave_t)),]))

widforpdf <- ncol(disbar)
data_matrix_other <- as.matrix(data_matrix_other)
if (dim(data_matrix_other)[2] ==1 ) data_matrix_other <- t(data_matrix_other)

if (is.null(data_matrix_other)==F) {
  disbar <- rbind(data_matrix_big,"Other"=c(colSums(data_matrix_other)),deparse.level = 2)
}

if (mean(colSums(disbar))<0.9999) {                         #Complete to 100%
  if (rownames(disbar)[nrow(disbar)]=="Other") {
    disbar[nrow(disbar),] <- disbar[nrow(disbar),]+(1-colSums(disbar))
  }
  else {
    disbar <- rbind(disbar,"other"=(sapply((1-colSums(disbar)),function(x)max(x,0))),deparse.level = 2)
  }
}
colours <- c(brewer.pal(9, "Set1"),brewer.pal(9, "Pastel1"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),brewer.pal(8, "Accent"),brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel2"),sample(rainbow(length(colnames(t(disbar)))),length(colnames(t(disbar)))))
data_melt <- melt(abs(data.matrix(t(disbar))),varnames=c("Samples","Cutline"),value.name="Relative_Abundance")
#-----------------------------------------------------------------------------------------
pp<-ggplot(data_melt,aes(x=Samples,y=Relative_Abundance,fill=Cutline))+
  geom_bar(stat='identity')+ ylab("Relative Abundance")+ xlab("Samples")+
  scale_x_discrete(limits=c(colnames(disbar)))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","25%","50%","75%","100%"))+
  guides(fill = guide_legend(ncol = (ceiling(nrow(disbar)/35))))+
  scale_fill_manual (values=colours) +
  theme(legend.position="right",axis.text.x=element_text(size=12,colour="black",angle=90,vjust=0.5),
        axis.text.y=element_text(size=12,colour="black"), axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),panel.grid.major=element_line(colour=NA))
suppressWarnings(ggsave(opts$outfile,plot=pp,width=ceiling(16+widforpdf/8),height=10, limitsize=FALSE))
