assemblies<-read.delim("allPan_references.tsv")
assemblies$genomeName <- sub("GCA_013052645.3_Mhudiblu_PPA_v2", "Mhudiblu", assemblies$genomeName)
assemblies$genomeName <- sub("GCF_002880755.1_Clint_PTRv2_genomic", "Clint", assemblies$genomeName)
assemblies$genomeName <- sub("GCF_000001405.40_GRCh38.p14_genomic", "GRCh38", assemblies$genomeName)
assemblies$genomeName <- sub("GCF_009914755.1_T2T-CHM13v2.0_genomic", "CHM13", assemblies$genomeName)
#assemblies$genomeName <- sub("mPanPan1", "Bonobo_T2T", assemblies$genomeName)
#assemblies$genomeName <- sub("mGorGor1", "Gorilla_T2T", assemblies$genomeName)
#assemblies$genomeName <- sub("mPanTro3", "Chimpanzee_T2T", assemblies$genomeName)
assemblies$type <- sub("^$", "primary", assemblies$type)
dataset <- read.delim("Pan_updated_Final.tsv", header=TRUE)
meta_statistics <- merge(assemblies,dataset,by=c("genomeName", "type"), all.x=TRUE)
NGx<-meta_statistics[meta_statistics $stat == 'NGx',]
​
​
​
colors <- c("black", "#9dced9", "#8F8F8F", "#026078", 
            "#FF450040", "#FFA50040", "#FF8C0040", "#FF7F5040", "#FF634740", "#FFD70040", "#FFA07A40", 
            "#FF7F0040", "#FF550040", "#FF993340", "#FF751840", "#FF572140", "#FF6F0040", "#FFB34740", 
            "#FF9F0040", "#FF681F40", "#FFAE4240")
​
​
plotNGx<-ggplot(data=NGx, aes(x = v1, y = v2, colour = NewName, linetype = type)) +
  geom_line() +  xlab("NG(x)") + ylab("Sequence length (bp)") + #scale_y_continuous(trans='log10') +
  scale_linetype_manual(values=c("dotted", "dotdash", "solid")) +
  scale_color_manual(values=colors) + guides(colour = guide_legend(title = c("Sample", "Type"))) +  
  ylim(1e4,NA) +
  theme_classic() 
ggsave("plotNGx.pdf", plotNGx, width = 10, height =4, dpi=300)
​
​
plotNGx_highligted<-ggplot(data=NGx, aes(x = v1, y = v2, colour = NewName, linetype = type)) +
  geom_line() +  xlab("NG(x)") + ylab("Sequence length (bp)") + scale_linetype_manual(values=c("dotted", "dotdash", "solid")) +
  ylim(1e4,NA) + guides(colour = guide_legend(title = c("Sample", "Type"))) +
  gghighlight(genomeName %in%  c("PR00115", "Clint", "Mhudiblu", "CHM13", "GRCh38"), use_direct_label = FALSE) +  
  scale_color_manual(values= c("pink", "green","black", "orange", "#9dced9")) +   
  theme_classic() +
  theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.title = element_blank()
)