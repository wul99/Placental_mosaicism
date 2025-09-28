library(tidyverse)
library(ggplot2)

data <- read.csv("SBS_signatures_decomposition_profile.csv")

tidy_data <- data %>%
  mutate(Signatures = str_split(Global.NMF.Signatures, " & ")) %>%
  unnest(Signatures) %>%
  mutate(Signature = str_extract(Signatures, "Subs-[0-9]+"), 
         Signature = str_replace(Signature, "Subs-0?", "SBS"),
         Percentage = as.numeric(str_extract(Signatures, "\\d+\\.\\d+")),
         Sample = str_extract(Sample.Names, "(?<=Sample ).*")) %>% 
  select(Sample, Signature, Percentage)

sample_order <- c("mosaic_embryos", "euploid_embryos", "natural_conceptions")
tidy_data <- tidy_data %>%
  mutate(Sample = factor(Sample, levels = sample_order))
tidy_data$Signature <- factor(tidy_data$Signature, levels = c("SBS1", "SBS5", "SBS18"))

pic <- ggplot(tidy_data, aes(x = Sample, y = Percentage, fill = Signature)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.6, width = 0.6) +
    theme_bw()+
    labs(y = "SBS signature proportions",fill = "Signature",x=NULL) +
    scale_x_discrete(labels = c("mosaic_embryos" = "Mosaic\nembryos","euploid_embryos" = "Euploid\nembryos","natural_conceptions" = "Natural\nconceptions"))+
    # theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.position = "top",panel.grid = element_blank(),text = element_text(size = 30),axis.text = element_text(size = 30,color = "black"),legend.text = element_text(size = 30,color = "black"),panel.border = element_rect(color = "black"))
