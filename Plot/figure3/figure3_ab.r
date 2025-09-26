library(nlme) 
library(multcomp)
library(ggsignif)
library(ggplot2)
library(RColorBrewer)
library(patchwork)


## figure 3a
snv_data <- read.csv("SNV-num-group-VAF-depth.txt",header=T,sep="\t")
snv_data$group <- factor(snv_data$group, levels = c("mosaic embryos", "euploid embryos", "natural conceptions"))
mod <- lme(number~group+depth, random=~1|family, data=snv_data)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))
pic1 <- ggplot(snv_data,aes(x = group,y = number)) +
    geom_boxplot(aes(fill = group),alpha = 0.7, width = 0.6)+
    scale_fill_manual(values = mycolor)+
    scale_color_manual(values = mycolor)+
    theme_bw()+
    scale_x_discrete(labels = c("mosaic embryos" = "Mosaic\nembryos","euploid embryos" = "Euploid\nembryos","natural conceptions" = "Natural\nconceptions"))+
    labs(y="SNV burden", x = NULL)+
    ylim(0,170)+
    theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.position = "none",panel.border = element_rect(color = "black"))+
    geom_signif(xmin = c(1, 2, 1),xmax = c(2, 3, 3),annotations = annotations, y_position = c(max(snv_data$number)+15,max(snv_data$number)+30,max(snv_data$number)+45),textsize = 10)

indel_data <- read.csv("indel-num-group-VAF-depth.txt",header=T,sep="\t")
indel_data$group <- factor(indel_data$group, levels = c("mosaic embryos", "euploid embryos", "natural conceptions"))
mod <- lme(number~group+depth, random=~1|family, data=indel_data)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))
pic2 <- ggplot(indel_data,aes(x = group,y = number)) +
    geom_boxplot(aes(fill = group),alpha = 0.7, width = 0.6)+
    scale_fill_manual(values = mycolor)+
    scale_color_manual(values = mycolor)+
    theme_bw()+
    scale_x_discrete(labels = c("mosaic embryos" = "Mosaic\nembryos","euploid embryos" = "Euploid\nembryos","natural conceptions" = "Natural\nconceptions"))+
    labs(y="Indel burden", x = NULL)+
    ylim(0,30)+
    theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.position = "none",panel.border = element_rect(color = "black"))+
    geom_signif(xmin = c(1, 2, 1),xmax = c(2, 3, 3),annotations = annotations, y_position = c(max(indel_data$number)+3,max(indel_data$number)+6,max(indel_data$number)+9),textsize = 10)

combined_plot <- pic1 + pic2 + plot_layout(ncol = 2, guides = "collect") + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))


## figure 3b
mut_data <- read.csv("mut-num-group-VAF-depth.txt",header=T,sep="\t")
mut_data$group <- factor(mut_data$group, levels = c("mosaic embryos", "euploid embryos", "natural conceptions"))
mod <- lme(median_VAF~group+depth, random=~1|family, data=mut_data)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))
pic3 <- ggplot(mut_data,aes(x = group,y = median_VAF)) +
    geom_boxplot(aes(fill = group),alpha = 0.7, width = 0.6)+
    scale_fill_manual(values = mycolor)+
    scale_color_manual(values = mycolor)+
    theme_bw()+
    scale_x_discrete(labels = c("mosaic embryos" = "Mosaic\nembryos","euploid embryos" = "Euploid\nembryos","natural conceptions" = "Natural\nconceptions"))+
    labs(y="Median VAF of mutations", x = NULL)+
    ylim(0,0.40)+
    theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.position = "none",panel.border = element_rect(color = "black"))+
    geom_signif(xmin = c(1, 2, 1),xmax = c(2, 3, 3),annotations = annotations, y_position = c(max(mut_data$median_VAF)+0.03,max(mut_data$median_VAF)+0.07,max(mut_data$median_VAF)+0.11),textsize = 10)
