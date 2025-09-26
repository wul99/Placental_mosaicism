library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(ggbreak)
library(patchwork)
library(ggpubr)


## figure 4a
data <- read.csv("all_mut.gnomad.af.txt", sep='\t', header=F,stringsAsFactors=FALSE)
data$V4[is.na(data$V4)] <- 0
data1 <- na.omit(data)
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))
data1$group <- factor(data1$V6,levels=c("ME","EE","NC"),ordered = TRUE)

groups <- unique(data1$group)
combs <- combn(groups, 2, simplify = FALSE)
wilcox_results <- lapply(combs, function(pair) {
  subset_data <- data1[data1$group %in% pair, ]
  test <- wilcox.test(V4 ~ group, data = subset_data, alternative="greater")
  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = test$p.value
  )
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$p_label <- paste0(wilcox_df$group1, " vs. ", wilcox_df$group2,": p = ", signif(wilcox_df$p_value, 3))
p_text <- paste(wilcox_df$p_label, collapse = "\n")

breaks <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1)
labels <- c("0", "1e-6~1e-5", "1e-5~1e-4", "1e-4~1e-3", "1e-3~1")

data_binned <- data1 %>%
  mutate(
    bin = cut(
      V4,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE,
      right = FALSE  
    )
  ) %>%
  count(group, bin) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))    

mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))

pic <- ggplot(data_binned,aes(x = bin,y = prop, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha=0.6) +  
#   geom_text(
#     aes(label = sprintf("%.1f%%", prop * 100),vjust = ifelse(group == "shared", -0.5, -1.5)), 
#     # position = position_dodge2(width = 0.7, preserve = "single"),
#     position = position_dodge(width = 0.7),
#     # vjust = 1.5, 
#     size = 8
#   ) +
  theme_bw() +
  scale_fill_manual(
    values = c("ME" = brewer.pal(9,"Blues")[5], "EE" = brewer.pal(9,"Greens")[4], "NC" = brewer.pal(9,"YlOrBr")[3]),
    labels = c("Mosaic\nembryos", "Euploid\nembryos", "Natural\nconceptions"),
    name = "Group"
  ) +
  scale_y_break(c(0.2,0.7),space = 0.4)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,0.9),labels = scales::percent) +
  labs(
    x = "Population allele frequency",
    y = "Proportion of mutations"
  ) +
  scale_x_discrete(
  labels = c(
    "0" = "0",
    "1e-6~1e-5" = expression(10^{-6}*"~"*10^{-5}),
    "1e-5~1e-4" = expression(10^{-5}*"~"*10^{-4}),
    "1e-4~1e-3" = expression(10^{-4}*"~"*10^{-3}),
    "1e-3~1"    = expression(10^{-3}*"~"*10^{0})
  ))+
  theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 30,color = "black"),axis.text.y = element_text(size = 30,color = "black"),
    panel.grid.major.x = element_blank(),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black")
  )+guides(fill = guide_legend())+
  theme(legend.position = "top",legend.justification = "left")

group_full <- c(ME = "Mosaic embryos", EE = "Euploid embryos", NC = "Natural conceptions")

table_df <- wilcox_df %>%
  mutate(
    group1 = group_full[group1],
    group2 = group_full[group2],
    p_value = sprintf("%.3f", p_value)
  ) %>%
  select(group1, group2, p_value)

table_df$comparison <- paste0("  ", table_df$group1," vs. ",table_df$group2, "  ")
table_df$p_value <- paste0(" ", table_df$p_value, " ")
table_df <- table_df[ ,c("comparison","p_value")]
colnames(table_df) <- c("Comparison", "P-value")

table_plot <- ggtexttable(table_df, rows = NULL,
                          theme = ttheme(
                            base_style = "classic",
                            base_size = 30,
                            padding = unit(c(4, 6), "mm"),
                            colnames.style = colnames_style(
                              color = "black", face = "bold", size = 30, fill = "gray90",linewidth = 1
                            ),
                            tbody.style = tbody_style(
                              color = "black", size = 30, fill = c("white")
                            )
                          ))

final_plot <- pic / table_plot + plot_layout(heights = c(5,  1))


## figure 4b
data <- read.csv("all_mut.gnomad.af.txt", sep='\t', header=F,stringsAsFactors=FALSE)
data$V4[is.na(data$V4)] <- 0
data1 <- na.omit(data)
data1$group <- factor(data1$V2,levels=c("BP_shared","PP_shared","P_unique"),ordered = TRUE)

groups <- unique(data1$group)
combs <- combn(groups, 2, simplify = FALSE)
wilcox_results <- lapply(combs, function(pair) {
  subset_data <- data1[data1$group %in% pair, ]
  test <- wilcox.test(V4 ~ group, data = subset_data, alternative="greater")
  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = test$p.value
  )
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$p_label <- paste0(wilcox_df$group1, " vs. ", wilcox_df$group2,":\np = ", signif(wilcox_df$p_value, 3))
p_text <- paste(wilcox_df$p_label, collapse = "\n")

breaks <- c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1)
labels <- c("0", "1e-6~1e-5", "1e-5~1e-4", "1e-4~1e-3", "1e-3~1")

data_binned <- data1 %>%
  mutate(
    bin = cut(
      V4,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE,
      right = FALSE  
    )
  ) %>%
  count(group, bin) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))    

pic <- ggplot(data_binned,aes(x = bin,y = prop, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha=0.7) +  
#   geom_text(
#     aes(label = sprintf("%.1f%%", prop * 100),vjust = ifelse(group == "shared", -0.5, -1.5)), 
#     # position = position_dodge2(width = 0.7, preserve = "single"),
#     position = position_dodge(width = 0.7),
#     # vjust = 1.5, 
#     size = 8
#   ) +
  theme_bw() +
  scale_fill_manual(
    values = c("BP_shared" = "#D74B34", "PP_shared" = "#F3AAA4", "P_unique" = "#4BB0C8"),
    labels = c("BP_shared", "PP_shared", "P_unique"),
    name = "Group"
  ) +
  scale_y_break(c(0.2,0.6),space = 0.4)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,0.9),labels = scales::percent) +
  labs(
    x = "Population allele frequency",
    y = "Proportion of mutations"
  ) +
  scale_x_discrete(
  labels = c(
    "0" = "0",
    "1e-6~1e-5" = expression(10^{-6}*"~"*10^{-5}),
    "1e-5~1e-4" = expression(10^{-5}*"~"*10^{-4}),
    "1e-4~1e-3" = expression(10^{-4}*"~"*10^{-3}),
    "1e-3~1"    = expression(10^{-3}*"~"*10^{0})
  ))+
  theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 30,color = "black"),axis.text.y = element_text(size = 30,color = "black"),
    panel.grid.major.x = element_blank(),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black")
  )+guides(fill = guide_legend())+
  theme(legend.position = "top",legend.justification = "left")

format_pval <- function(p) {
  if (p < 0.001) {
    formatted <- formatC(p, format = "e", digits = 1)
    formatted <- gsub("e–", " x 10^–", formatted)
    return(formatted)
  } else {
    return(formatC(p, format = "f", digits = 3))
  }
}
# wilcox_df$p_label <- sapply(wilcox_df$p_value, format_pval)
table_df <- wilcox_df %>%
  mutate(
    p_value = sapply(wilcox_df$p_value, format_pval)
  ) %>%
  select(group1, group2, p_value)

table_df$comparison <- paste0("       ", table_df$group1," vs. ",table_df$group2, "       ")
table_df$p_value <- paste0("     ", table_df$p_value, "     ")
table_df <- table_df[ ,c("comparison","p_value")]
colnames(table_df) <- c("Comparison", "P-value")

table_plot <- ggtexttable(table_df, rows = NULL,
                          theme = ttheme(
                            base_style = "classic",
                            base_size = 30,
                            padding = unit(c(4, 6), "mm"),
                            colnames.style = colnames_style(
                              color = "black", face = "bold", size = 30, fill = "gray90",linewidth = 1
                            ),
                            tbody.style = tbody_style(
                              color = "black", size = 30, fill = c("white")
                            )
                          ))

final_plot <- pic / table_plot + plot_layout(heights = c(5,  1))
