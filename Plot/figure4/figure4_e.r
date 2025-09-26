library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hs1)
library(ggplot2)
library(reshape2) 
library(tidyr)
library(do)
library(patchwork)
library(lsa)
library(dplyr)

new_cosmic <- read.csv("COSMIC_v3.4_SBS_GRCh37.txt",sep="\t",row.names=1)

all_mut <- read.csv("all.SNV.bed",sep="\t",header=F)

all_mut$group <- factor(all_mut$V4, levels=c("BP_shared","PP_shared","P_unique"), order=TRUE)

all_sig <- mut.to.sigs.input(mut.ref = all_mut, sample.id = "group", chr = "V1", pos = "V3", ref = "V5", alt = "V6", bsg = BSgenome.Hsapiens.UCSC.hs1)

w=lapply(unique(all_mut$group) , function(i){
  sample_1 = whichSignatures(tumor.ref = all_sig[,], 
                             signatures.ref = as.data.frame(t(new_cosmic)), 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  print(i)
  return(sample_1$tumor)
})

w=do.call(rbind,w)

data <- as.data.frame(t(w))
data$type <- rownames(data)
data$mut <- mid(data$type,3,3)

mutation_order <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mycolor <- c("skyblue","black","red","gray","lightgreen","pink")
# data$mut <- factor(data$mut,levels = c("C>A","C>G","C>T","T>A","T>C","T>G"))
# data$type <- factor(data$type,levels = data$type)
sorted_data <- data %>% mutate(mut = factor(mut, levels = mutation_order), type = factor(type, levels = unique(type[order(mut)]))) %>% arrange(mut, type)
data_long <- pivot_longer(sorted_data, cols = c(BP_shared,PP_shared,P_unique), names_to = "category", values_to = "value")
data_long$category <- factor(data_long$category,levels = c("BP_shared", "PP_shared", "P_unique"))
# labeller = labeller(category = c("BP_shared" = "BP_shared (n=171)", "PP_shared" = "PP_shared (n=379)", "P_unique" = "P_unique (n=2577)"))
p1 <- ggplot(data_long, aes(x = type, y = value, fill=mut)) +
  geom_bar(stat = "identity", width = 0.7) + theme_bw() +
  facet_wrap(~category, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = mycolor)+
  scale_x_discrete(labels = paste0(mid(data$type,1,1),mid(data$type,3,1),mid(data$type,7,1)))+
  labs(y = "Fraction of SNVs",x=NULL) +
  ylim(0,0.09)+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 18, color = "black"),axis.text.y = element_text(size = 25, color = "black"),strip.background = element_rect(fill = "#f5f5f5", linewidth = 0.5, color = "black"),strip.text = element_text(size = 30,color = "black"),panel.grid = element_blank(),text = element_text(size = 30),panel.border = element_rect(color = "black"),legend.position = "top",legend.text = element_text(size = 30,color = "black")) + guides(fill = guide_legend(nrow = 1))

sim1 <- cosine(data_long[data_long$category=="BP_shared",]$value, data_long[data_long$category=="PP_shared",]$value)
sim2 <- cosine(data_long[data_long$category=="BP_shared",]$value, data_long[data_long$category=="P_unique",]$value)
sim3 <- cosine(data_long[data_long$category=="PP_shared",]$value, data_long[data_long$category=="P_unique",]$value)

combined_plot <- p1 + plot_annotation(title = paste0("Cosine Similarity","\n","BP_shared vs. PP_shared: ",round(sim1[1], 3),"\n","BP_shared vs. P_unique: ",round(sim2[1], 3),"\n","PP_shared vs. P_unique: ",round(sim3[1], 3)), theme = theme(plot.title = element_text(size = 30)))
