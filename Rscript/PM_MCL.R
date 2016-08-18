## install necessary libraries
p <- c("networkD3","ProNet" ,"magrittr","optparse")
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
option_list <- list(
  make_option(c("-i", "--inputfile"), type="character", help="Input cytoscape_file [Required]"),
  make_option(c("-o", "--outfile"), type="character", default='MCL_network.html', help="Output Network [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$inputfile)) stop('Please check the input file.')

corr_data <- read.table(opts$inputfile)
local <- data.frame(corr_data$V1,corr_data$V2)
corr <- data.frame(corr_data$V3)
network <- construct_local(input = local, node.attribute = corr)
cluster_out <- cluster(network, method = "MCL", plot = FALSE)
names <- names(cluster_out)
groups <- cluster_out
names(groups) <- NULL

nodes <- data.frame(c(1), names, c(1), c(1))
names(nodes) <- c("id", "name", "group", "size")
#nodes to group
j <- 1
for (node_name in nodes$name) {
  group_id <- which(names == node_name)
  #print(groups[group_id])
  nodes$group[j] <- groups[group_id]
  j <- j+1
}

link <- MisLinks
links <- corr_data
names(links) <- c("source","target","value")
#links$value <- round(links$value * 10)
#links$value <- abs(links$value)
links$source <- as.character(links$source)
links$target <- as.character(links$target)

j <- 1;
for (i in links$source) {
  id <- which(names == i)
  links$source[j] <- id
  #print(j)
  #print(links$source[j])
  j <- j+1
}
j <- 1
for (i in links$target) {
  id <- which(names == i)
  links$target[j] <- id
  j <- j+1
}

links$source <- as.numeric(links$source)
links$target <- as.numeric(links$target)
links$source <- links$source - 1
links$target <- links$target - 1 

nd <- forceNetwork(Links = links, Nodes = nodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 1)

saveNetwork(network = nd,file = opts$outfile)

#