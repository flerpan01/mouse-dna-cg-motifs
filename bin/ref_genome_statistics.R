# Extract statistics from the reference genome, GRCm39
# CG motifs (BSseq)
# CCGG motifs (RRBS)
# Distribution of gene types
# Motifs coverage inside CGIs and genes


library(dplyr)
library(purrr)
library(ggplot2)
library(forcats)
library(scales)

# import mouse ref genome, mm39
library("BSgenome.Mmusculus.UCSC.mm39")
genome <- BSgenome.Mmusculus.UCSC.mm39

# Reorder chromosomes, un-annotated last in the list
unanno_chrom <- grepl("_", names(genome))
#chromosomes <- c(names(genome)[!unanno_chrom], names(genome)[unanno_chrom])
chromosomes <- names(genome)[!unanno_chrom]

foo <- function(chrom, motif){
  #motif <- "CG"
  d <- matchPattern(motif, genome[[chrom]]) %>%
    as.data.frame()
  d$chr <- chrom
  d$pos <- paste0(d$chr, ":", d$start)
  d[, c("chr", "start", "end", "pos")]
}

cg_ref <- map_dfr(chromosomes, foo, motif = "CG")
ccgg_ref <- map_dfr(chromosomes, foo, motif = "CCGG")

tab <- function(dat){
  data <- data.frame(table(dat$chr))
  colnames(data) <- c("chr", "count")
  data <- data[match(chromosomes, data$chr),]
  data$chr <- fct_inorder(data$chr)

  return(data)
}

plotter <- function(data, title = NULL, caption = NULL){
  nudge <- max(data$count)
  data$nudge_y <- ifelse(data$count > nudge / 4, -nudge / 18, nudge / 22)

  ggplot(data, aes(chr, count, label = comma(count))) +
    geom_col(fill = "lightblue", color = "black") +
    geom_text(nudge_y = data$nudge_y) +
    coord_flip() + scale_x_discrete(limits = rev(chromosomes)) +
    scale_y_continuous(label = comma) +
    theme_bw(20) + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Chromosome", y = "", title = title, caption = caption)
}

cg <- tab(cg_ref)
ccgg <- tab(ccgg_ref)

plotter(cg, title = "Number CG-motifs", caption = "Mouse ref. genome, GRCm39")
plotter(ccgg, title = "Number CCGG-motifs", caption = "Mouse ref. genome, GRCm39")

sum(ccgg$count) / sum(cg$count)

