res<-read.table("tumor response post.vs.pre pc count with batch.csv",header = T,sep = ",",check.names = F,row.names = 1)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
xl <- -1
xr <- 1
yp <- -log10(0.05)
dif.data<- na.omit(res) %>% 
  mutate(logP = -log10(padj)) %>% 
  mutate(color = ifelse(log2FoldChange > xr & logP > yp, 
                        yes = "Tumor", no = ifelse(log2FoldChange < xl & logP > yp, yes = "Normal",  no = "none")))
table(dif.data$color)
###extract top differentially expressed genes
top_n <- 10
up <- arrange(subset(dif.data, color == "Tumor"),desc(log2FoldChange)) 
down <- arrange(subset(dif.data, color == "Normal"),log2FoldChange)
top_labelled <- rbind.data.frame(up[1:top_n,], down[1:top_n,])

ggplot(dif.data, aes(x = log2FoldChange, y = logP)) + 
  geom_point(aes(color = factor(color)), size = 1.55, alpha = 0.8, na.rm = TRUE) + # add gene points
  theme_bw(base_size = 16) + # clean up theme
  theme(legend.position = "none") + # remove legend 
  theme(axis.text=element_text(size=14),axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        legend.title=element_text(face="bold",size=10),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  ylab('-log10 (Pvalue)')+xlab('mRNA (Log2 FC, post-treatment vs pre-treatment)')+
  #xlim(-24,12)+
  geom_vline(xintercept = xl, colour = "black", linetype = "dashed") + # add line at 0
  geom_vline(xintercept = xr, colour = "black", linetype = "dashed") + # add line at 0
  geom_hline(yintercept = yp, colour = "black", linetype = "dashed") + # p(0.05) = 1.3
  scale_color_manual(values = c("Tumor" = "#E64B35", 
                                "Normal" = "#3182bd", 
                                "none" = "#636363")) + # change colors
  annotate(geom = "text", 
           label = "pre", 
           x = -10, y = 19, 
           size = 5, colour = "#3182bd") + # add Down text
  annotate(geom = "text", 
           label = "post", 
           x = 8, y = 19, 
           size = 5, colour = "#E64B35") + # add Up text
  scale_y_continuous(trans = "log1p")+  # Scaled Y-axis with log1p function
  geom_label_repel(data = top_labelled, 
                   max.overlaps =100,
                   aes(label = gene), fill = "white",
                   fontface = 'bold',box.padding = 0.15, color = 'black',
                   label.size = 0.1,point.padding = 0.4,segment.color = 'gold')