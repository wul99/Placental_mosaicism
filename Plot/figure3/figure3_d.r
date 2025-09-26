library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(stringr)
library(scales)
library(coin)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "snv.bed",
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

input <- args$input
output_path <- args$output

mut_file <- file.path(input)
raw_lines <- readLines(mut_file)
processed_lines <- gsub(":", "\t", raw_lines)  
mut_data <- read.delim(text = processed_lines, header = FALSE, sep = "\t") %>%
    setNames(c("chr", "start", "end", "family", "type", "AF", "origin","sample","group")) %>%
    filter(origin %in% c("F", "M")) %>% 
    filter(sample != "B_unique") 

mut_data_number <- mut_data %>%
  group_by(family, origin) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(
    nesting(family),  
    origin = c("F", "M"),  
    fill = list(count = 0) 
  ) %>%
  mutate(origin = factor(origin, levels = c("F", "M"), ordered = TRUE))

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

print(mut_data_number,n=100)
test_result <- wilcox_test(count ~ origin, data = mut_data_number, paired = TRUE, distribution = "exact", alternative = "greater")
p_values <- pvalue(test_result)
annotations <- ifelse(p_values < 0.001, "***",ifelse(p_values < 0.01, "**",ifelse(p_values < 0.05, "*", "ns")))
p <- ggplot(mut_data_number, aes(x = origin, y = count)) +
    geom_boxplot(
    aes(fill = origin),
    alpha = 0.7,
    width = 0.6,          
    color = "black",       
    outlier.color = "black"
  ) +
  scale_fill_manual(values = c("F" = adjust_color("#317EC2",1.2,1.1,1), "M" = adjust_color("#C35B54",1.1,1.1,1))) +
  labs(y = "Number of mutations") +
  geom_signif(xmin = 1,
    xmax = 2,
    tip_length = 0.02,
    annotations = annotations, 
    y_position = max(mut_data_number$count)+max(mut_data_number$count)/10,
    textsize = 10) +
  scale_y_continuous(
    breaks = pretty_breaks(),  
    labels = number_format(accuracy = 1),  
    expand = c(0.1, 0.1)
  ) +
  scale_x_discrete(
  labels = c("F" = "Paternal","M" = "Maternal")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    legend.position = "none",
    axis.text = element_text(size = 30,color = "black"),
    axis.title = element_text(size = 30,color = "black"),
    axis.title.x = element_blank()
  )

ggsave(p,file=output_path, width=7, height=10, dpi = 600)