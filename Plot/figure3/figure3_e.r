library(ggplot2)
library(RColorBrewer)
library(dplyr)

data <- read.csv("difSNVtype-proportion-all.txt", sep='\t', header=TRUE,stringsAsFactors=FALSE)

data_summary <- data %>%
  group_by(mut_type, family) %>%
  summarise(
    mean_proportion = first(all_proportion),  
    sd = sd(proportion),                      
    lower = mean_proportion - sd,
    upper = mean_proportion + sd,
    .groups = "drop" 
  )

mycolor1 <- c(brewer.pal(9,"Blues")[2:6])
mycolor2 <- c(brewer.pal(9,"Greens")[2:4])
mycolor3 <- c(brewer.pal(9,"YlOrBr")[2:4])
mycolor <- c(mycolor1,mycolor2,mycolor3)
data_summary$family <- factor(data_summary$family,levels=c("M1","M2","M3","M4","M5","N1","N2","N3","C1","C2","C3"),ordered = TRUE)

# data$sample <- factor(data$sample,levels=c("P1","P2","P3","P4","P5"),ordered = TRUE)
# shape_mapping <- c(P1=15, P2=17, P3=18, P4=10, P5=20)

###bar_plot
pic <- ggplot(data_summary, aes(x=mut_type,fill=family))+
geom_bar(aes(y=mean_proportion,fill = family), stat='identity', position=position_dodge(width=0.9))+
# geom_point(aes(y = mean_proportion, fill = family), position = position_dodge(width = 0.9),shape = 21,size = 6, color = "black", stroke = 1)+
# geom_point(aes(y=proportion, shape=sample, group=interaction(mut_type, family)), position = position_dodge(width=0.9), size=0.7, show.legend=TRUE, alpha = 0.8)+
# scale_shape_manual(values=shape_mapping)+
scale_fill_manual(values=mycolor)+
geom_errorbar(aes(ymin = lower, ymax = upper),position = position_dodge(width = 0.9), width = 0.4, color = "black",linewidth = 0.25) +
theme_bw()+
theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black")) +
labs(y="SNV type proportions", x = NULL)
