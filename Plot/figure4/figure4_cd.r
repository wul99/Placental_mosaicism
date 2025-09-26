library(ggplot2)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(tidyr)
library(dplyr)
library(rstatix)
library(ggsignif)

args <- commandArgs(trailingOnly = TRUE)

data1 <- read.csv(paste0("data/BP_shared.",args[1],".txt"),sep="\t",header = F)
data2 <- read.csv(paste0("data/PP_shared.",args[1],".txt"),sep="\t",header = F)
data3 <- read.csv(paste0("data/P_unique.",args[1],".txt"),sep="\t",header = F)
data4 <- read.csv(paste0("data/BP_random.",args[1],".txt"),sep="\t",header = F)
data5 <- read.csv(paste0("data/PP_random.",args[1],".txt"),sep="\t",header = F)
data6 <- read.csv(paste0("data/P_random.",args[1],".txt"),sep="\t",header = F)

data1$group <- "BP_shared"
data2$group <- "PP_shared"
data3$group <- "P_unique"
data4$group <- "BP_random"
data5$group <- "PP_random"
data6$group <- "P_random"

all_data <- bind_rows(data1, data2, data3, data4, data5, data6)
all_data$group <- factor(all_data$group, levels=c("BP_shared","PP_shared","P_unique","BP_random","PP_random","P_random"))

all_data$alpha_group <- ifelse(all_data$group %in% c("BP_shared", "PP_shared", "P_unique"), 0.7, 0.4)
mycolor <- c("#D74B34", "#F3AAA4", "#4BB0C8", "#D74B34", "#F3AAA4", "#4BB0C8")

all_data <- all_data[all_data$V3>50,]

if (args[1]=="phyloP30way"){
        my_comparisons <- list(c("BP_shared", "PP_shared"), c("PP_shared", "P_unique"),c("BP_shared", "P_unique"),c("BP_shared", "BP_random"), c("PP_shared", "PP_random"),c("P_unique", "P_random"))
        stat_df <- all_data %>% pairwise_wilcox_test(V6 ~ group,alternative = "less",p.adjust.method = "none",comparisons = my_comparisons)
        stat_df$annotations <- ifelse(stat_df$p < 0.05, stat_df$p.adj.signif, paste0("p = ", signif(stat_df$p, 2)))
        signif_data <- data.frame(xmin = c(1, 2, 1, 1, 2, 3),xmax = c(2, 3, 3, 4, 5, 6),annotations = stat_df$annotations, y_position = c(1.1, 1.36, 1.62, 1.88, 2.14, 2.4), textsize = c(10, 9, 10, 10, 9, 10))
        p <- ggplot(all_data, aes(x = group, y = V6, fill = group, alpha = alpha_group)) +
                geom_boxplot(outlier.size = 0.5, width = 0.7) +
                theme_bw() +
                labs(x = NULL, y = "PhyloP30way score") +
                scale_fill_manual(values = mycolor)+
                scale_alpha_identity() +
                ylim(-1,2.5)+
                theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none",panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black"))
        for(i in 1:nrow(signif_data)) {
                p <- p + geom_signif(
                xmin = signif_data$xmin[i],
                xmax = signif_data$xmax[i],
                annotations = signif_data$annotations[i],
                y_position = signif_data$y_position[i],
                textsize = signif_data$textsize[i],
                tip_length = 0.02
                )
                }
}else if (args[1]=="gc5Base"){
        my_comparisons <- list(c("BP_shared", "PP_shared"), c("PP_shared", "P_unique"),c("BP_shared", "P_unique"),c("BP_shared", "BP_random"), c("PP_shared", "PP_random"),c("P_unique", "P_random"))
        stat_df <- all_data %>% pairwise_wilcox_test(V6 ~ group,alternative = "greater",p.adjust.method = "none",comparisons = my_comparisons)
        print(stat_df)
        stat_df$annotations <- ifelse(stat_df$p < 0.05, stat_df$p.adj.signif, paste0("p = ", signif(stat_df$p, 2)))
        signif_data <- data.frame(xmin = c(1, 2, 1, 1, 2, 3),xmax = c(2, 3, 3, 4, 5, 6),annotations = stat_df$annotations, y_position = c(80, 89, 98, 107, 116, 125), textsize = c(10, 9, 10, 9, 9, 9))
        p <- ggplot(all_data, aes(x = group, y = V6, fill = group, alpha = alpha_group)) +
                geom_boxplot(outlier.size = 0.5, width = 0.7) +
                theme_bw() +
                labs(x = NULL, y = "GC content") +
                scale_fill_manual(values = mycolor)+
                scale_alpha_identity() +
                ylim(0,128)+
                theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none",panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black"))
        for(i in 1:nrow(signif_data)) {
                p <- p + geom_signif(
                xmin = signif_data$xmin[i],
                xmax = signif_data$xmax[i],
                annotations = signif_data$annotations[i],
                y_position = signif_data$y_position[i],
                textsize = signif_data$textsize[i],
                tip_length = 0.02
                )
                }
}else{
        my_comparisons <- list(c("BP_shared", "PP_shared"), c("PP_shared", "P_unique"),c("BP_shared", "P_unique"),c("BP_shared", "BP_random"), c("PP_shared", "PP_random"),c("P_unique", "P_random"))
        stat_df <- all_data %>% pairwise_wilcox_test(V6 ~ group,alternative = "greater",p.adjust.method = "none",comparisons = my_comparisons)
        stat_df$annotations <- ifelse(stat_df$p < 0.05, stat_df$p.adj.signif, paste0("p = ", signif(stat_df$p, 2)))
        signif_data <- data.frame(xmin = c(1, 2, 1, 1, 2, 3),xmax = c(2, 3, 3, 4, 5, 6),annotations = stat_df$annotations, y_position = c(1.1, 1.23, 1.36, 1.49, 1.62, 1.75), textsize = c(10, 9, 10, 10, 9, 10))
        p <- ggplot(all_data, aes(x = group, y = V6, fill = group, alpha = alpha_group)) +
                geom_boxplot(outlier.size = 0.5, width = 0.7) +
                theme_bw() +
                labs(x = NULL, y = "Depletion rank") +
                scale_fill_manual(values = mycolor)+
                scale_alpha_identity() +
                ylim(0,1.81)+
                theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none",panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black"))
        for(i in 1:nrow(signif_data)) {
                p <- p + geom_signif(
                xmin = signif_data$xmin[i],
                xmax = signif_data$xmax[i],
                annotations = signif_data$annotations[i],
                y_position = signif_data$y_position[i],
                textsize = signif_data$textsize[i],
                tip_length = 0.02
                )
                }
}
