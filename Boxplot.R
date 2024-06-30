tpm<-read.table("final_pctpm_Combat_convertID.txt",sep = ",",header=T,check.names=F,row.names = 1)
gene<-c("ASCL1","NEUROD1","POU2F3","YAP1")
tlsmarker1<-tpm[as.vector(gene),as.vector(batch$sample)]
tlsmarker1<-as.data.frame(t(as.matrix(tlsmarker1)))
tlsmarker1$treatment<-batch$tissue[match(rownames(tlsmarker1),batch$sample)]
library(reshape2)
ciber<-melt(tlsmarker1,id.vars = "treatment")
colnames(ciber)[2]<-"cell"
ciber$treatment<-factor(ciber$treatment,levels = c("pre","post"))
table(ciber$treatment)
library(ggpubr)
p=ggboxplot(ciber, x="cell", y="value", color = "treatment",
            ylab="Expression",
            xlab="",
            palette = c("blue","red"))
p+rotate_x_text(45)
