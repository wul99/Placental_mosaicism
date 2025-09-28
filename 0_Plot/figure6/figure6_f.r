library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(stringr)
library(scales)
library(RColorBrewer)
library(nlme) 
library(multcomp)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "LOH.bed",
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", 
              metavar = "character"),
  make_option(c("-d", "--depth_file"), type = "character", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

input <- args$input
output_path <- args$output
depth_file <- args$depth_file

mut_file <- file.path(input)
raw_lines <- readLines(mut_file)
processed_lines <- gsub("[,:]+","\t", raw_lines)  
mut_data <- read.delim(text = processed_lines, header = FALSE, sep = "\t") %>%
  setNames(c("chr", "start", "end", "family", "type", "B","P1","P2","P3","P4","P5", "origin","sample_group","group")) %>%
  filter(sample_group != "B_unique") %>%
  pivot_longer(
    cols = c("B","P1","P2","P3","P4","P5"),  
    names_to = "sample",           
    values_to = "AF"             
    ) 

mut_data_AF <- mut_data %>%
  group_by(family, sample) %>%  
  filter(!is.na(AF)) %>%
  filter(AF != "0") %>%
  filter(sample != "B") %>%
  summarise(median_AF = median(AF, na.rm = TRUE),group = first(group), .groups = "keep") %>%
  mutate(group = factor(group, levels = c("ME", "EE","NC"), ordered = TRUE)) %>%
  left_join(depth_file, by = c("family", "sample"))

print(mut_data_AF,n=100)

adjust_color <- function(hex, sat_factor, light_factor, alpha) {
  require(colorspace)

  hcl_col <- hex2RGB(hex) %>% 
    as("polarLUV") %>% 
    coords() %>%
    as.data.frame()
  
  adjusted_hcl <- polarLUV(
    H = hcl_col$H,
    C = hcl_col$C * sat_factor,
    L = hcl_col$L * light_factor
  )

  hex_alpha <- hex(adjusted_hcl) %>% 
    paste0(., as.hexmode(round(alpha * 255)) %>% str_pad(2, pad = "0"))
  
  return(hex_alpha)
}
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))

mod <- lme(median_AF~group+depth, random=~1|family, data=mut_data_AF)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
p <- ggplot(mut_data_AF, aes(x = group, y = median_AF)) +
  geom_boxplot(aes(fill = group),alpha = 0.7, width = 0.6) +
  scale_fill_manual(values = mycolor) +
  labs(y = "LOH mosaic rate") +
  geom_signif(
    tip_length = 0.02,
    xmin = c(1, 2, 1),
    xmax = c(2, 3, 3),
    annotations = annotations, 
    y_position = c(max(mut_data_AF$median_AF)+0.05,max(mut_data_AF$median_AF)+0.1,max(mut_data_AF$median_AF)+0.15),
    textsize = 10) +
  scale_x_discrete(
    labels = c("ME" = "Mosaic\nembryos","EE" = "Euploid\nembryos","NC" = "Natural\nconceptions")) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    expand = c(0, 0.05)
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    legend.position = "none",
    axis.text = element_text(size = 30,color = "black"),
    axis.title = element_text(size = 30,color = "black"),
    axis.title.x = element_blank()
  )

ggsave(p,file=output_path, width=8, height=10, dpi = 600)