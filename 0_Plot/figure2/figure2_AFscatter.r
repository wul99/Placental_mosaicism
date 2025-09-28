library(optparse)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggrastr)

width_base_px <- 18000
width_base <- width_base_px/1200

option_list <- list(
  make_option(c("-f", "--family"), type = "character", default = "M1",
              metavar = "character"),
  make_option(c("-r", "--region"), type = "character", default = "chr14",
              metavar = "character"),
  make_option(c("-i", "--input"), type = "character",   # ${phased_input} generate by generate_phased_input.sh
              metavar = "character"),            
  make_option(c("-o", "--output"), type = "character", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

family <- args$family
region <- args$region
input_path <- args$input
output_path <- args$output
samples <- list(
      order = c('B','P1','P2','P3','P4','P5','E'),
      number = 7
    )

parse_region <- function(region_str) {
  if (str_detect(region_str, ":")) {
    parts <- str_split(region_str, "[:-]", simplify = TRUE)
    chrom <- parts[1,1]
    start <- as.numeric(parts[1,2])
    end <- as.numeric(parts[1,3])
    length <- end - start
    return(list(
      chrom = chrom,
      start = start - length,
      end = end + length
    ))
  } else {
    return(list(
      chrom = region_str,
      start = 0,
      end = 999999999
    ))
  }
}

read_sample_data <- function(family, region_info, input_path) {
  path <- file.path(input_path)
  
  read.delim(path, header = FALSE, na.strings = "") %>%
    select(where(~ !all(is.na(.)))) %>%
    select(1:!!(samples$number + 2)) %>%
    filter(V2 >= region_info$start & V2 <= region_info$end) %>%
    setNames(c("chr", "pos", samples$order)) %>%
    pivot_longer(
    cols = 3:!!(samples$number + 2),  
    names_to = "sample",           
    values_to = "F_AF"             
    ) %>%
    mutate(M_AF = 1 - F_AF) %>%
    pivot_longer(
    cols = c(F_AF, M_AF),
    names_to = "origin",
    values_to = "AF",
    values_drop_na = TRUE
    ) %>%
    mutate(
    origin = factor(origin, levels = c("F_AF", "M_AF")),
    sample = factor(sample, 
      levels = c('E','B','P1','P2','P3','P4','P5','V'),
      labels = c('TE biopsy','Umbilical cord blood','Placental sample 1','Placental sample 2','Placental sample 3','Placental sample 4','Placental sample 5','Post-abortion chorionic villus')),
    pos_mb = pos / 1e6
    )
}

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

region_info <- parse_region(region)
data <- read_sample_data(family, region_info, input_path)

p <- ggplot(data, aes(x = pos_mb, y = AF, color = origin, group = sample)) +
  geom_point_rast(
      size = 1/sqrt(nrow(data)),  
      raster.dpi = getOption("ggrastr.default.dpi", 600)
  ) +
  scale_size_continuous(
      range = c(0.00001, 0.1),  
      limits = c(0, 1)
  ) +  
  scale_color_manual(
      name = NULL,
      values = c("F_AF" = adjust_color("#317EC2",1,1,0.2), "M_AF" = adjust_color("#C35B54",1,1,0.2)),
      labels = c("Paternal", "Maternal")
  ) +
  guides(
      color = guide_legend(
          override.aes = list(
              title = NULL,
              shape = 15 ,       
              size = width_base/3,         
              alpha = 1         
          )
      )
  ) +
  scale_x_continuous(
      name = sprintf("Position on chromosome %s (Mb)", str_remove(region_info$chrom,"chr")),
      breaks = scales::breaks_pretty(5)
  ) +
  scale_y_continuous(
      name = "AF",
      breaks = seq(0, 1, 0.5)
  ) +
  theme_bw()+
  theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black"),
      text = element_text(size = width_base),
      axis.text = element_text(size = width_base),
      strip.text = element_text(size = width_base),
      legend.text = element_text(size = width_base),
      legend.position = "top",
      legend.box.spacing = margin(t = 0),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "#f5f5f5", linewidth = 0.5)
  ) +
  facet_wrap(~ sample, ncol = round(samples$number/2), nrow = 2)

ggsave(plot = p, 
  file = output_path,
  width = width_base, height = width_base * samples$number / 15,
  dpi = 600,
  bg = "transparent", 
  device = cairo_pdf)
