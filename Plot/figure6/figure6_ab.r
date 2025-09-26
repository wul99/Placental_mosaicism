library(nlme) 
library(multcomp)
library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(dplyr)


## figure 6a
data <- read.csv("SV-num-depth.txt",header=T,sep="\t")
data$group <- factor(data$group, levels = c("mosaic embryos", "euploid embryos", "natural conceptions"))
mod <- lme(number~group+depth, random=~1|family, data=data)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))
pic <- ggplot(data,aes(x = group,y = number)) +
    geom_boxplot(aes(fill = group),alpha = 0.7, width = 0.6)+
    scale_fill_manual(values = mycolor)+
    scale_color_manual(values = mycolor)+
    theme_bw()+
    scale_x_discrete(labels = c("mosaic embryos" = "Mosaic\nembryos","euploid embryos" = "Euploid\nembryos","natural conceptions" = "Natural\nconceptions"))+
    labs(y="SV burden",x=NULL)+
    theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.position = "none",panel.border = element_rect(color = "black"))+
    geom_signif(xmin = c(1, 2, 1),xmax = c(2, 3, 3),annotations = annotations, y_position = c(max(data$number)+0.9,max(data$number)+1.8,max(data$number)+2.7),textsize = 10)+
    ylim(0,9)


## figure 6b
data <- read.csv("SV-VAF-depth.txt",header=T,sep="\t", row.names = NULL)
data$group <- factor(data$group, levels = c("mosaic embryos", "euploid embryos", "natural conceptions"))

result <- data %>%
    group_by(family, sample) %>%
    summarise(
        median_vaf = median(VAF, na.rm = TRUE),
        group = first(group),
        depth = first(depth),
        .groups = "drop"
)

mod <- lme(median_vaf~group+depth, random=~1|family, data=result)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))
pic2 <- ggplot(result,aes(x = group,y = median_vaf)) +
    geom_boxplot(aes(fill = group),alpha = 0.7, width = 0.6)+
    scale_fill_manual(values = mycolor)+
    scale_color_manual(values = mycolor)+
    theme_bw()+
    scale_x_discrete(labels = c("mosaic embryos" = "Mosaic\nembryos","euploid embryos" = "Euploid\nembryos","natural conceptions" = "Natural\nconceptions"))+
    labs(y="Median VAF of SVs", x = NULL)+
    ylim(0,0.5)+
    theme(panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.position = "none",panel.border = element_rect(color = "black"))+
    geom_signif(xmin = c(1, 2, 1),xmax = c(2, 3, 3),annotations = annotations, y_position = c(max(result$median_vaf)+0.04,max(result$median_vaf)+0.08,max(result$median_vaf)+0.12),textsize = 10)
