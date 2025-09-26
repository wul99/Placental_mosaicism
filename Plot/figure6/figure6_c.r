library(ggplot2)
library(RColorBrewer)
library(scales)

data <- read.csv("SV.info.final.txt",header=F,sep="\t")
data$Mechanism <- factor(data$V13, levels = c('alt-EJ','NHEJ','TEI'))  
mycolor <- c("#769fcd","#b9d7ea","#d6e6f2")
# data$length <- log(data$V12,base=10)
breaks <- 10^seq(0, 7, by = 1)
pic <- ggplot(data, aes(x = V12, fill = Mechanism)) +
  geom_histogram(binwidth = 1,breaks = breaks, color = "white", alpha = 0.8) +
  scale_fill_manual(values = mycolor)+ 
  theme_bw()+
  labs(x = "SV length (bp)", y = "Number of SVs") +
  # scale_x_continuous(breaks = seq(0, 8, by = 1), labels = seq(0, 8, by = 1)) +
  scale_x_log10(limits = c(1e1, 1e7),breaks = 10^seq(1, 7, by = 1),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(breaks = seq(0, 25, by = 5))+
  # coord_cartesian(xlim = c(1, 7)) +
  theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30),panel.border = element_rect(color = "black"))+theme(legend.position = c(0.85, 0.89))
