count<-read.table("count.csv",sep = "\t",header=T,check.names=F,row.names = 1)
batch<-batch[order(batch$tissue,decreasing = T),]
data<-count[,batch$sample]
library(DESeq2)
colData <- batch
rownames(colData)<-batch$sample
dds <- DESeqDataSetFromMatrix(data, colData, design= ~tissue+batch)
dds <- DESeq(dds)
res1 = results(dds, contrast=c("tissue","post","pre"))
