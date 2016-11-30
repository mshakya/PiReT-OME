library(limma)

pairs <- read.table(file = "./Deseq_EdgeRpairs.txt",sep = "\t", stringsAsFactors = FALSE,  header=TRUE)
pair1 <- pairs$group1
pair2 <-pairs$group2

for (n in 1:length(pair1) ) {
if (length(pair1) ==0) next
filename=paste(pair1[n],pair2[n], sep="__")
filename_sig=paste(pair1[n],pair2[n],"sig", sep="__")


fileout_EdgeR=paste('./EdgeR/EdgeR.',filename, ".txt", sep="")
fileout_sig_EdgeR=paste('./EdgeR/EdgeR.',filename_sig, ".txt", sep="")

fileout_Deseq=paste('./Deseq/Deseq.',filename, ".txt", sep="")
fileout_sig_Deseq=paste('./Deseq/Deseq.',filename_sig, ".txt", sep="")

result_EdgeR <- read.delim(file =fileout_sig_EdgeR, stringsAsFactors = FALSE, skip=1)
result_Deseq <- read.delim(file =fileout_sig_Deseq, stringsAsFactors = FALSE, skip=1)

set1<-result_Deseq[,1]
set2<-result_EdgeR[,1]

# What are the possible letters in the universe?
universe <- sort(unique(c(set1, set2)))
#universe <- sort(unique(c(set1, set2, set3)))

# Generate a matrix, with the sets in columns and possible letters on rows
Counts <- matrix(0, nrow=length(universe), ncol=2)
# Populate the said matrix
for (i in 1:length(universe)) {
   Counts[i,1] <- universe[i] %in% set1
   Counts[i,2] <- universe[i] %in% set2
 #  Counts[i,3] <- universe[i] %in% set3
}

# Name the columns with the sample names
#colnames(Counts) <- c("set1","set2","set3")
colnames(Counts) <- c("Deseq","EdgeR")

# Specify the colors for the sets
cols<-c("Red", "Green")

fileout_vennDiagram=paste('./figures/vennDiagram.',filename, ".pdf", sep="")

pdf(fileout_vennDiagram)
vennDiagram(vennCounts(Counts),  circle.col=cols)
title(main = list(filename_sig, cex = 1.5, col = "red", font = 3))
dev.off()

}


