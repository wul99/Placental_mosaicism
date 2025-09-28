library(ggplot2)
library(RColorBrewer)
library(dplyr)

data <- read.csv("cell_division.txt", sep='\t', header=TRUE,stringsAsFactors=FALSE)

mycolor1 <- c(brewer.pal(9,"Blues")[3:7])
mycolor2 <- c(brewer.pal(9,"Greens")[3:4])
mycolor3 <- c(brewer.pal(9,"YlOrBr")[4:5])
mycolor <- c(mycolor1,mycolor2,mycolor3)
data$family <- factor(data$family,levels=c("M1","M2","M3","M4","M5","N1","N2","N3","C1","C2","C3"),ordered = TRUE)

pic <- ggplot(data, aes(x=family, fill = family))+
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, group = interaction(family, node)),position = position_dodge(width = 0.4), width = 0.1, color = "black", linewidth = 0.6) +
    geom_point(aes(y = mean_proportion, group = interaction(family, node), fill = family), position = position_dodge(width = 0.4),shape = 21,size = 8, color = "black", stroke = 1)+
    scale_fill_manual(values=mycolor)+
    scale_x_discrete(expand = expansion(add = 0.4))+
    theme_bw()+
    theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"), legend.position="none",panel.border = element_rect(color = "black")) +
    labs(y="Cell divisions", x = NULL)
