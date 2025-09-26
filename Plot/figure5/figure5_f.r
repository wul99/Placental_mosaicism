library(ggplot2)
library(dplyr)
library(RColorBrewer)

data1 <- read.csv("clonal_nodes_and_mutations.txt", sep='\t', header=T,stringsAsFactors=FALSE)

data1 <- data1 %>%
  group_by(family, lineage) %>%
  arrange(family, lineage, node) %>% 
  mutate(cumulative_mut = cumsum(mut)) %>% 
  ungroup()

data1$family <- factor(data1$family, levels=c("M1","M2","M3","M4","M5","N1","N2","C2","C3"))

mycolor1 <- c(brewer.pal(9,"Blues")[2:6])
mycolor2 <- c(brewer.pal(9,"Greens")[2:3])
mycolor3 <- c(brewer.pal(9,"YlOrBr")[3:4])
mycolor <- c(mycolor1,mycolor2,mycolor3)


p1 <- ggplot(data1, aes(x = node, y = cumulative_mut, color = family, group = interaction(family, lineage))) +
        # annotate("rect", 
        #    xmin = 1, xmax = 4,
        #    ymin = -Inf, ymax = Inf,
        #    alpha = 0.2, fill = "gray80") +
        geom_line(linewidth=1.5) + 
        geom_point(size=4) + 
        scale_color_manual(values = mycolor)+
        labs(x = "Node depth", y = "Cumulative mutations") +
        theme_bw()+
        scale_x_continuous(breaks = 0:6)+
        theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black")) + geom_vline(xintercept = c(1.75, 4.25), color = "gray40", linetype = "dashed", size = 1)
