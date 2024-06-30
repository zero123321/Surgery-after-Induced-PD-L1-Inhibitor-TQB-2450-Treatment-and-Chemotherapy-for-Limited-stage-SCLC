library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name,gene_symbol)
head(KEGG_df)
HALL_df = msigdbr(species = "Homo sapiens",category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)
head(HALL_df)
ALL_df<-rbind(KEGG_df,HALL_df)

###response
result1<-read.table("tumor post response.vs.nonresponse pc count with batch.csv",header = T,sep = ",",check.names = F,row.names = 1)
ge = result1$log2FoldChange
names(ge) = result1$gene
ge = sort(ge,decreasing = T)
head(ge)
em <- GSEA(ge, TERM2GENE = ALL_df,pvalueCutoff = 0.05)
