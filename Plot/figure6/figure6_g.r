library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "LOH.bed",
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
    filter(sample != "B_unique") %>%
    mutate(length = log10(end - start)) %>%
    mutate(group = factor(group, levels = c("ME", "EE","NC"), ordered = TRUE))
print(mut_data)
mycolor <- c(c(brewer.pal(9,"Blues")[5]),c(brewer.pal(9,"Greens")[4]),c(brewer.pal(9,"YlOrBr")[3]))

pic <- ggplot(mut_data, aes(x = length, fill = group)) +
  geom_histogram(binwidth = 0.5,breaks = seq(0, 10, by = 0.5), color = "white", alpha = 0.7) +
  scale_fill_manual(values = mycolor,
    labels = c("ME" = "Mosaic embryos","EE" = "Euploid embryos","NC" = "Natural conceptions"))+ 
  theme_bw()+
  labs(x = "LOH length (bp)", y = "Number of LOHs") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), labels = parse(text = paste0("10^", 0:10))) +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  coord_cartesian(xlim = c(3, 8)) +
  theme(panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    text = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.position = c(0.75,0.88))

ggsave(pic,file=output_path, width = 10, height = 10, dpi = 600)

