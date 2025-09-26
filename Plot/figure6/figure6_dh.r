library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(optparse)

option_list <- list(
  make_option(c("-t", "--mut_type"), type = "character",
              metavar = "character"),
  make_option(c("-i", "--input"), type = "character", 
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

mut_type <- args$mut_type
input_path <- args$input
output_path <- args$output

df <- read.delim(input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

plot_data <- df %>%
  mutate(annotation = factor(annotation,levels = c('intergenic','intron','exon'))) %>%
  mutate(star = case_when(
    pvalue >= 0.01 & pvalue < 0.05 ~ "*",
    pvalue >= 0.001 & pvalue < 0.01 ~ "**",
    pvalue < 0.001 ~ "***",
    TRUE ~ ""
  )) %>%
  mutate(
    label_pos = ifelse(l2fold > 0, l2fold + 0.1, l2fold - 0.1)
  ) %>%
  mutate(
    CI95_low_ratio = log2(observed / CI95high),
    CI95_high_ratio = log2(observed / CI95low)
  )

point_to_mm <- 30 * 0.3527
p <- ggplot(plot_data, aes(x = annotation, y = l2fold)) +
  geom_col(fill = "#8594B6",width = 0.6) +
  geom_text(aes(y = label_pos, label = star),
            color = "black",size = point_to_mm, vjust = -0.2,hjust = 1) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "black", linewidth = 1) +
  geom_errorbar(
    aes(ymin = CI95_low_ratio, ymax = CI95_high_ratio),  
    width = 0.2,       
    color = "#2c3e50",  
    linewidth = 1,      
    alpha = 1        
  ) +
  labs( 
       x = NULL,
       y = sprintf("Log2 (O/E) of %s", mut_type)) + 
  scale_x_discrete(
    labels = c("intron" = "Intron","exon" = "Exon","intergenic" = "Intergenic")) +
  scale_y_continuous(
    breaks = seq(round(min(plot_data$CI95_low_ratio)),round(max(plot_data$CI95_high_ratio)),by = (round(max(plot_data$CI95_high_ratio))-round(min(plot_data$CI95_low_ratio))/2))
  ) +  
  theme_bw() +  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black"),
        text = element_text(size = 30),
        axis.text = element_text(size = 30),
        legend.position = "none") +
  coord_flip() 

ggsave(output_path, plot = p, 
       width = 10, height = 10, 
       unit = "in",
       dpi = 600)