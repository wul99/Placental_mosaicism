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

mut_data_number <- mut_data %>%
  distinct(family, group) %>%
  bind_rows(tibble(family = "C3", group = "NC")) %>%  
  right_join(
    expand_grid(
      family = c(unique(mut_data$family),"C3"),
      sample = unique(filter(mut_data, sample != "B")$sample)
    ),
    by = "family"
  ) %>%
  left_join(
    mut_data %>%
      filter(sample != "B", !is.na(AF), AF != "0") %>%
      group_by(family, sample) %>%
      summarise(count = n(), .groups = "drop"),
    by = c("family", "sample")
  ) %>%
  mutate(
    count = replace_na(count, 0),
    group = factor(group, levels = c("ME", "EE", "NC"), ordered = TRUE)
  ) %>%
  arrange(family, sample) %>%
  left_join(depth_file, by = c("family", "sample"))

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

print(mut_data_number,n=100)

mod <- lme(count~group+depth, random=~1|family, data=mut_data_number)
anova(mod, type='marginal')
summary(glht(mod, linfct = mcp(group = "Tukey")))
p_values <- as.numeric(summary(glht(mod, linfct = mcp(group = "Tukey")))$test$pvalues)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
p <- ggplot(mut_data_number, aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group),alpha = 0.7,width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = mycolor) +
  labs(x = "", y = "LOH burden") +
  geom_signif(
    tip_length = 0.005,
    xmin = c(1, 2, 1),
    xmax = c(2, 3, 3),
    annotations = annotations, 
    y_position = c(max(mut_data_number$count)-13.4,max(mut_data_number$count)-12.8,max(mut_data_number$count)-12.2),
    textsize = 10) +
  scale_x_discrete(
    labels = c("ME" = "Mosaic\nembryos","EE" = "Euploid\nembryos","NC" = "Natural\nconceptions")) +
  scale_y_continuous(
    breaks = seq(0, 7, by = 1),  
    labels = number_format(accuracy = 1)  
  ) +
  coord_cartesian(ylim = c(0, 5)) +
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