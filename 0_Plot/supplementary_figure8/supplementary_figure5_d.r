pacman::p_load(tidyverse, readxl, vcfR, VariantAnnotation, rtracklayer, Biostrings, Rsamtools, patchwork, plyranges, genomation, dplyr)


data = read_tsv("all_mut.txt")

make_vcf <- function(data){
    c <- data$chr
    s <- data$pos 
    e <- as.integer(s+sapply(as.character(data$ref),nchar, simplify="array")-1)
    g <- GRanges(c, IRanges(s,e))
    v <- VCF(rowRanges=g)
    ref(v) <- DNAStringSet(data$ref)
    alt(v) <- DNAStringSetList(lapply(data$alt, function(x){x}))
    mcols(v)$subject = data$family
    mcols(v)$VAF = data$AF
    mcols(v)$group = data$group
    return(v)
}

v_all = make_vcf(data)
v_shared = make_vcf(data %>% filter(group=="shared"))
v_unique = make_vcf(data %>% filter(group!="shared"))


## use hs1 gtf
gtf_file <- "$HOME/gene_expression/placenta.gtf"

transcripts.gr <- import(gtf_file)
transcripts.gr <- transcripts.gr[transcripts.gr$type == "transcript"]

seqlevelsStyle(transcripts.gr) <- "UCSC"

ref_genome <- "BSgenome.Hsapiens.UCSC.hs1"
library(ref_genome, character.only = TRUE)
df = as.data.frame(transcripts.gr)
df$transcription_start_site <- df$start
tss.gr = GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$transcription_start_site, end=df$transcription_start_site),strand=df$strand)

## using the generic TSS regions
peak_plot <- function(){
    TSSbs_regions = tss.gr
    TSSbs_regions = TSSbs_regions %>% reduce_ranges()
    TSSbs_regions_center <- resize(TSSbs_regions, width=1, fix="center")
    TSSbs_regions_upstream <- resize(TSSbs_regions_center, width=5000, fix="end")
    TSSbs_regions_downstream <- resize(TSSbs_regions_center, width=5000, fix="start")

    bin_num = 5
    x1 = suppressWarnings(ScoreMatrixList(c(shared=granges(v_shared),unique=granges(v_unique)), bin.num=bin_num, bin.op = "max", windows=TSSbs_regions_upstream,strand.aware=TRUE))
    x2 = suppressWarnings(ScoreMatrixList(c(shared=granges(v_shared),unique=granges(v_unique)), bin.num=bin_num, bin.op = "max", windows=TSSbs_regions_downstream,strand.aware=TRUE))

    total_shared <- sum(rev(colSums(x2$shared))+colSums(x1$shared))
    total_unique <- sum(rev(colSums(x2$unique))+colSums(x1$unique))

    print(rev(colSums(x2$shared))+colSums(x1$shared))
    print(rev(colSums(x2$unique))+colSums(x1$unique))
    
    total_shared1 <- dim(v_shared)[1]
    total_unique1 <- dim(v_unique)[1]

    plot_df = data.frame(shared = ((rev(colSums(x2$shared))+colSums(x1$shared)) / total_shared1), 
                        unique= ((rev(colSums(x2$unique))+colSums(x1$unique)) / total_unique1), 
                        tss_dist = seq(4500,500,length.out=bin_num)) 

    bin_size = 10000/bin_num

    p_values <- sapply(1:bin_num, function(i) {

    observed <- matrix(
        c((rev(colSums(x2$shared))+colSums(x1$shared))[i], total_shared1 - (rev(colSums(x2$shared))+colSums(x1$shared))[i],
        (rev(colSums(x2$unique))+colSums(x1$unique))[i], total_unique1 - (rev(colSums(x2$unique))+colSums(x1$unique))[i]),
        nrow = 2,
        byrow = F
    )
    fisher.test(observed)$p.value
    })

    print(p_values)

    plot_df <- plot_df %>%
    mutate(
      shared_low = pmax(shared - 1.96 * sqrt(shared * (1 - shared) / total_shared1), 0),
      shared_high = pmin(shared + 1.96 * sqrt(shared * (1 - shared) / total_shared1), 1),
      unique_low = pmax(unique - 1.96 * sqrt(unique * (1 - unique) / total_unique1), 0),
      unique_high = pmin(unique + 1.96 * sqrt(unique * (1 - unique) / total_unique1), 1)
    )

    p1 = plot_df %>%
    pivot_longer(
        cols = c(shared, unique),
        names_to = "group",
        values_to = "proportion"
    ) %>%
    mutate(
        low_ci = ifelse(group == "shared", shared_low, unique_low),
        high_ci = ifelse(group == "shared", shared_high, unique_high)
    ) %>%
    ggplot(aes(x = factor(tss_dist), y = proportion, color = group)) + 
    geom_pointrange(
    aes(ymin = low_ci, ymax = high_ci, color=group), size=2, linewidth = 2, position=position_dodge(width=0.6), alpha=1) +
    theme_bw() + labs(color="",x= "Distance from TSS (bp)",y="Proportion of mutations") +
    scale_color_manual(breaks=c("shared", "unique"),labels=c("Shared mutations", "Unique mutations"),values=c("#ee827c", "#89c3eb")) +
    theme(legend.position = "top",
    panel.grid = element_blank(),
    axis.text = element_text(size = 30,color = "black"),
    text = element_text(size = 30),
    legend.text = element_text(size = 30,color = "black"),
    legend.title = element_text(size = 30,color = "black"),
    panel.border = element_rect(color = "black")
    )

    return(p1)
}

p1 = peak_plot()



data = read_tsv("all_mut.txt")

make_vcf <- function(data){
    c <- data$chr
    s <- data$pos 
    e <- as.integer(s+sapply(as.character(data$ref),nchar, simplify="array")-1)
    g <- GRanges(c, IRanges(s,e))
    v <- VCF(rowRanges=g)
    ref(v) <- DNAStringSet(data$ref)
    alt(v) <- DNAStringSetList(lapply(data$alt, function(x){x}))
    mcols(v)$subject = data$family
    mcols(v)$VAF = data$AF
    mcols(v)$group = data$group
    return(v)
}

v_all = make_vcf(data)
v_BP = make_vcf(data %>% filter(dif_stage=="BP_shared"))
v_PP = make_vcf(data %>% filter(dif_stage=="PP_shared"))
v_P = make_vcf(data %>% filter(dif_stage=="P_unique"))

## using the generic TSS regions
peak_plot <- function(){
    TSSbs_regions = tss.gr
    TSSbs_regions = TSSbs_regions %>% reduce_ranges()
    TSSbs_regions_center <- resize(TSSbs_regions, width=1, fix="center")
    TSSbs_regions_upstream <- resize(TSSbs_regions_center, width=5000, fix="end")
    TSSbs_regions_downstream <- resize(TSSbs_regions_center, width=5000, fix="start")

    bin_num = 5
    x1 = suppressWarnings(ScoreMatrixList(c(BP=granges(v_BP),PP=granges(v_PP),P=granges(v_P)), bin.num=bin_num, bin.op = "max", windows=TSSbs_regions_upstream,strand.aware=TRUE))
    x2 = suppressWarnings(ScoreMatrixList(c(BP=granges(v_BP),PP=granges(v_PP),P=granges(v_P)), bin.num=bin_num, bin.op = "max", windows=TSSbs_regions_downstream,strand.aware=TRUE))

    print(rev(colSums(x2$BP))+colSums(x1$BP))
    print(rev(colSums(x2$PP))+colSums(x1$PP))
    print(rev(colSums(x2$P))+colSums(x1$P))

    total_BP <- dim(v_BP)[1]
    total_PP <- dim(v_PP)[1]
    total_P <- dim(v_P)[1]

    plot_df = data.frame(BP = ((rev(colSums(x2$BP))+colSums(x1$BP)) / total_BP), 
                        PP = ((rev(colSums(x2$PP))+colSums(x1$PP)) / total_PP), 
                        P = ((rev(colSums(x2$P))+colSums(x1$P)) / total_P),
                        tss_dist = seq(4500,500,length.out=bin_num)) 


    bin_size = 10000/bin_num

    groups <- c("BP", "PP", "P")
    group_pairs <- combn(groups, 2, simplify = FALSE) 

    results <- lapply(group_pairs, function(pair) {
        group1 <- pair[1]
        group2 <- pair[2]
        
        total1 <- get(paste0("total_", group1))
        total2 <- get(paste0("total_", group2))

        p_values <- sapply(1:bin_num, function(i) {
            observed <- matrix(
            c(
                (rev(colSums(x2[[group1]]))+colSums(x1[[group1]]))[i], total1 - (rev(colSums(x2[[group1]]))+colSums(x1[[group1]]))[i],
                (rev(colSums(x2[[group2]]))+colSums(x1[[group2]]))[i], total2 - (rev(colSums(x2[[group2]]))+colSums(x1[[group2]]))[i]
            ),
            nrow = 2
            )
            fisher.test(observed)$p.value
        })

    list(
        comparison = paste(group1, "vs", group2),
        p_values = p_values,
        p_adjusted = p.adjust(p_values, method = "BH")
    )
    })

    print(results)

    plot_df <- plot_df %>%
    mutate(
      BP_low = pmax(BP - 1.96 * sqrt(BP * (1 - BP) / total_BP), 0),
      BP_high = pmin(BP + 1.96 * sqrt(BP * (1 - BP) / total_BP), 1),

      PP_low = pmax(PP - 1.96 * sqrt(PP * (1 - PP) / total_PP), 0),
      PP_high = pmin(PP + 1.96 * sqrt(PP * (1 - PP) / total_PP), 1),

      P_low = pmax(P - 1.96 * sqrt(P * (1 - P) / total_P), 0),
      P_high = pmin(P + 1.96 * sqrt(P * (1 - P) / total_P), 1)
    )

    p1 = plot_df %>%
    pivot_longer(
        cols = c(BP, PP, P),
        names_to = "group",
        values_to = "proportion"
    ) %>% mutate(group = factor(group, levels = c("BP", "PP", "P"))) %>%
    mutate(
        low_ci = ifelse(group == "BP", BP_low, ifelse(group == "PP", PP_low, P_low)),
        high_ci = ifelse(group == "BP", BP_high, ifelse(group == "PP", PP_high, P_high))
    ) %>%
    ggplot(aes(x = factor(tss_dist), y = proportion, color = group)) + 
    geom_pointrange(
    aes(ymin = low_ci, ymax = high_ci, color=group), size=2, linewidth = 2, position=position_dodge(width=0.6), alpha=1) +
    theme_bw() + labs(color="",x= "Distance from TSS (bp)",y="Proportion of mutations") +
    scale_color_manual(breaks=c("BP", "PP", "P"),labels=c("BP_shared", "PP_shared", "P_unique"),values= c(BP="#D26C5A", PP="#F3AAA4", P="#89c3eb"))+ 
    theme(legend.position = "top",
    panel.grid = element_blank(),
    axis.text = element_text(size = 30,color = "black"),
    text = element_text(size = 30),
    legend.text = element_text(size = 30,color = "black"),
    legend.title = element_text(size = 30,color = "black"),
    panel.border = element_rect(color = "black")
    )

    return(p1)
}

p2 = peak_plot()
