library(dplyr)
library(ggplot2)
library(ggforce)
library(stringr)
library(ggrepel) 

data <- read.csv("dif-mut-func-num.txt",header=F,sep="\t")
data$group <- factor(data$V5,levels=c("mosaic embryos","euploid embryos", "natural conceptions"),ordered = TRUE)

data <- data %>%
  mutate(
    V5 = str_to_sentence(V5),
    group = factor(V5, levels=c("Mosaic embryos", "Euploid embryos", "Natural conceptions")),
    type = factor(V1)
)

data$type <- factor(data$V1)
df <- data %>% group_by(group, type) %>% summarise(count = sum(V2),.groups = 'drop') %>% group_by(group) %>% mutate(percentage = count/sum(count)*100) %>% ungroup()

df$type <- factor(df$type, levels=c("Intergenic","Gene_flanking","Intron","UTR","Synonymous","Missense","LoF"),ordered = TRUE)

df <- df %>% arrange(group, type)

type_radius <- data.frame(
  type = c("Intergenic", "Gene_flanking", "Intron", "UTR", "Synonymous", "Missense","LoF"),
  radius_mult = c(1,1,1,1,1.5,1.5,1.5),
  dis_x = c(0,-0.08,-0.04,-0.09,0.08,0.08,0.1),
  dis_y = c(0,0,0,0,-0.1,-0.1,-0.5)
)

label_df <- left_join(df,df %>% 
  group_by(group) %>% summarise(Cnt_total = sum(percentage))) %>% left_join(type_radius, by = "type")  %>%  group_by(group) %>%
  mutate(
    end_angle = 2*pi*cumsum(percentage)/Cnt_total,
    start_angle = lag(end_angle, default = 0),
    mid_angle = 0.5*(start_angle + end_angle),
    label = case_when(
      percentage > 0.3 ~ sprintf("%.2f%%", percentage),  
      percentage < 0.3 ~ sprintf("%s\n%.2f%%", type, percentage)),
    # label = sprintf("%s\n%.2f%%", type, percentage),
    xend = radius_mult * sin(mid_angle)+dis_x,
    yend = radius_mult * cos(mid_angle)+dis_y
  ) %>%
  ungroup() 
  # %>% filter(type %in% c("Missense", "LoF"))

color_palette <- colorRampPalette(c("#F5DCD6", "#FF0000"))(nlevels(data$type))
names(color_palette) <- levels(df$type)

pic <- ggplot() +
  geom_arc_bar(
    data = df,
    stat = "pie",
    aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = percentage, fill = type),
    color = "white", linewidth = 0.2
  ) +
  # geom_segment(
  #   data = label_df,
  #   aes(x = 0.5*sin(mid_angle), y = 0.5*cos(mid_angle),
  #       xend = 0.5*xend, yend = 0.5*yend),
  #   color = "gray35", linewidth = 0.4) +
  # geom_label_repel(
  #   data = label_df,
  #   aes(x = 0.5*xend, y = 0.5*yend, label = label),
  #   size = 8, 
  #   min.segment.length = 0,
  #   segment.color = "gray35",
  #   color = "gray35",
  #   box.padding = 0.4,
  #   nudge_y = 0,
  #   direction = "both",
  #   force=0,
  #   fill = NA,
  #   label.size = NA
  #   ) +
  facet_wrap(~ group, strip.position = "top", nrow = 2, ncol = 2) + 
  scale_fill_manual(values = color_palette,breaks = levels(df$type)) +
  theme_void() +
  theme(
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.25),
    text = element_text(size = 30),
    strip.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    strip.text = element_text(size = 30, margin = margin(2,0,3,0))
)

## add labels in AI
