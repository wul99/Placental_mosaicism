library(tidyverse)
library(ggplot2)

data <- read.table("mutations_from_M3.txt",header=T,stringsAsFactors=FALSE)

data$mut <- factor(data$mut, levels = data$mut)
data$mut_label <- paste0("mut", seq_len(nrow(data)))
data$mut_label <- factor(data$mut_label, levels = data$mut_label)

data_long <- data %>%
  pivot_longer(
    cols = -c(mut,id,mut_label),
    names_to = "Sample",
    values_to = "AF"
  ) %>%
  mutate(mut_label = factor(mut_label, levels = unique(data$mut_label)), AF_label = sprintf("%.2f", AF))


pic <- ggplot(data_long, aes(x = Sample, y = mut_label, fill = AF)) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_text(
    aes(label = AF_label),
    size = 10, 
    show.legend = FALSE,
    color = "gray30"
  )+
  scale_fill_gradientn(
    colors = colorRampPalette(c("white", "#95C0DB"))(100),
    limits = c(0, 0.5), 
    na.value = "blue",     
    guide = guide_colorbar(barheight = unit(10, "cm"), breaks = seq(0, 0.5, by = 0.1), title.vjust = 5)
  ) +
  scale_y_discrete(limits = levels(data$mut_label)) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Allele\nfrequency"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),text = element_text(size = 30),axis.text.y = element_text(size = 30,color = "black"),axis.text.x = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black")
  )
