# R-Script to compute different transcriptional output from SLAMSeq data

# Copyright (c) 2020 Tobias Neumann
# Based on nf-core/slamseq pipeline (on April 2023)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(tidyverse)
        library(DESeq2)
        library(gridExtra)
        library(RColorBrewer)
        library(ggrepel)
})


######################
# Snakemake parse
######################
countFiles <- snakemake@input[["countFiles"]]
sample_man <- snakemake@input[["sample_manifest"]]

data_outdir <- snakemake@output[["data_outdir"]] # for contrast tables
fig_outdir <- snakemake@output[["fig_outdir"]] # for plots
ddsRDS <- snakemake@output[["ddsRDS"]]

subgroup <- snakemake@params[["subgroup_filter"]] # group of interest ADD CONDITONAL TO ALL GROUPS
pval_th <- snakemake@params[["pval_threshold"]]

######################
# Functions
######################

# Function to parse count table

parseSample <- function(file) {

  # con <- file(file,"r")
  # name <- readLines(con,n=1)
  # name = sub(".*name:","",name)
  # close(con)
  name = tail(unlist(strsplit(file, "/")), n=1)
  name = unlist(strsplit(name, "[.]"))[1]
  countTab = read_tsv(file)

  countTab = countTab %>%
    mutate(RPM = ReadCount / sum(ReadCount) * 1e6) %>%
    dplyr::select(GeneName, ReadCount, TcReadCount, RPM)


  names(countTab) = c("GeneName", paste0(name,"_total"),paste0(name,"_tc"), paste0(name,"_RPM"))

  return(countTab)
}

# Function to run DESeq2

deAnalysis <- function(counts, sample_info) {

  countData.total <- counts %>%
    dplyr::select(contains("_total")) %>%
    as.matrix

  row.names(countData.total) <- counts$GeneName

  countData <- counts %>%
    dplyr::select(contains("_tc")) %>%
    as.matrix

  row.names(countData) <- counts$GeneName

  sampleOrder = sub("_tc","",colnames(countData))

  sample_info = sample_info[match(sampleOrder, sample_info$Sample_name),]
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = sample_info,
                                design = ~ Treatment)

  row.names(dds) = row.names(countData)

  dds.total <- DESeqDataSetFromMatrix(countData = countData.total,
                                      colData = sample_info,
                                      design = ~ Treatment)

  dds.total <- dds.total[ rowSums(counts(dds)) > 0, ] #  remove uninformative rows

  dds <- dds[ rowSums(counts(dds)) > 0, ] #  remove uninformative rows

  #  Run deseq main command

  dds.total <- DESeq(dds.total)
  sizeFactors(dds) <- sizeFactors(dds.total) # apply size factors to tc counts
  dds <- DESeq(dds)

  return(dds)
}

# Get DESeq2 contrast

getContrast <- function(counts, dds, case, control) {

  results <- results(dds, contrast = c("condition",
                                       case,
                                       control))

  results <- data.frame(GeneName = as.character(rownames(dds)),
                        log2FC_deseq2 = results$log2FoldChange,
                        padj = results$padj)

  results <- results[which(results$padj != "NA"),]

  #  Calculate mean baseline expression (mean RPM of untreated samples).

  ctrl = control
  ctrl = sample_info %>%
    dplyr::filter(condition == ctrl) %>%
    .$Sample_name %>%
    paste0(.,"_RPM")

  avg.RPM <- data.frame(GeneName = as.character(counts$GeneName),
                        avg.RPM.ctrl = counts %>%
                          dplyr::select(all_of(ctrl)) %>%
                          rowMeans
  )

  #  Intersect with baseline expression and gene symbols and export as table.

  export.deseq2 <- left_join(results,avg.RPM) %>%
    dplyr::select(GeneName, everything()) %>%
    arrange(padj,log2FC_deseq2)
}

MAPlot <- function(export.deseq2, case, control, cutoff) {

  ###  Export MA-plots of differential gene expression calling

  #  Extract EntrezIDs for (top 20) deregulated genes for plotting

  dereg <- export.deseq2[which(export.deseq2$padj <= cutoff), ] # all dereg. genes
  down <- nrow(dereg[dereg$log2FC_deseq2 <= 0, ]) # downregulated genes
  up <- nrow(dereg[dereg$log2FC_deseq2 >= 0, ]) # upregulated genes
  top.dereg <- dereg[order(abs(dereg$log2FC_deseq2), decreasing = T)[1:20], 1]
  top.dereg <- top.dereg[!is.na(top.dereg)]

  #  Create generic data frame for plotting by ggplot2.
  #  d - gray-scale density coloring for scatter plots.

  x <- log10(export.deseq2$avg.RPM.ctrl) # x-axis: baseline mRNA expression
  y <- export.deseq2$log2FC_deseq2 # y-axis: log2 fold-change treated/ctrl
  GeneName <- export.deseq2$GeneName
  d <- densCols(x, y, nbin = 100,
                colramp = colorRampPalette((brewer.pal(9,"Greys")[-c(1:5)])))
  df <- data.frame(x, y, d, GeneName)

  ### Generate unlabeled plot with accurate marginal density for use in figures.

  # Generate basic MA-like plot with density coloring.

  p <- ggplot(df, aes(x = x, y = y, label = GeneName)) +
    theme_classic() +
    scale_color_identity() +
    labs(x = "average expression total mRNA (RPM)", y = "log2FC")

  p.basic <- p +
    geom_point(aes(x, y, col = d), size = 1.3, shape =16) +
    geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
    geom_hline(yintercept = 0, size = 0.8) +
    geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)

  #  Generate marginal density plot of fold-changes.

  p.right <- ggplot(df, aes(y)) +
    geom_density(alpha = .5, fill = "gray40", size = 0.8) +
    labs(x = "", y = "") +
    coord_flip() +
    geom_vline(xintercept = 0) +
    theme_classic()

  grid.arrange(p.basic, p.right,
               ncol = 2, widths = c(4, 1)) # assemble plots


  ###  Generate summary of MA-plots with additional information & highlights.
  #  NB: Formatting of axes may shift marginal density plot from scatter-plot.
  #  Only use exported plots for visual inspection.

  #  Generate basic MA-like plot with density coloring.

  p <- ggplot(df, aes(x = x, y = y, label = GeneName)) +
    theme_classic() +
    scale_color_identity() +
    labs(x = "average expression total mRNA (RPM)", y = "log2FC") +
    ggtitle(paste0(case," vs ",control),
            paste("n =", nrow(df),
                  "// n(up) =", up,
                  "// n(down) =", down)) +
    theme(axis.line = element_line(size = 0.5),
          axis.text	 = element_text(size = 12),
          axis.ticks	= element_line(size = 0.5))

  p.basic <- p +
    geom_point(aes(x, y, col = d), size = 1.3, shape =16) +
    geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
    geom_hline(yintercept = 0, size = 0.8) +
    geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)

  #  Generate marginal density plot of fold-changes.

  p.right <- ggplot(df, aes(y)) +
    geom_density(alpha=.5, fill="gray40", size = 0.8) +
    coord_flip() +
    geom_vline(xintercept = 0) +
    theme_classic() +
    ggtitle("","") +
    theme(legend.position = "none",
          text = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.ticks.x = element_line(color = "white"),
          axis.line.x = element_line(color = "white"))

  #  Generate MA-plot with highlights and labeling of significantly dereg. genes.

  if(nrow(dereg) > 0){

    p.highlight <- p +
      geom_point(data = df[!(df$GeneName %in% dereg$GeneName),],
                 aes(x, y, col = "gray60"), size = 1.3, shape =16) +
      geom_point(data = df[df$GeneName %in% dereg$GeneName,],
                 aes(x, y, col = "red1"), size = 1.3, shape =16) +
      geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
      geom_hline(yintercept = 0, size = 0.8) +
      geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)

    p.highlight.2 <- p.highlight +
      geom_label_repel(data = df[df$GeneName %in% top.dereg,])

    #  Export plot.

    grid.arrange(p.basic, p.right,
                 ncol = 2, widths = c(4,1)) # assemble plots
    grid.arrange(p.highlight, p.right,
                 ncol = 2, widths = c(4,1)) # add highlights
    grid.arrange(p.highlight.2,p.right,
                 ncol = 2, widths = c(4,1)) # add labels
  }
}

######################
# Read sample_info + counts
######################
cat("Reading sample counts...",  sep="\n")

sample_info = read.table(sample_man, header=TRUE, sep="\t")

sample_info$Sample_name <- paste(sample_info$Sample_type, "_", sample_info$Treatment, "_Bio-rep_", sample_info$Bio_rep, sep="")

if (nrow(sample_info) == length(unique(sample_info$Treatment))) {
  quit(save = "no", status = 0, runLast = TRUE)
}

counts = parseSample(file.path(countFiles[1]))

if(length(countFiles) > 1) {

  for (i in 2:length(countFiles)) {
    counts = counts %>% inner_join(parseSample(file.path(countFiles[i])), by = "GeneName") # may give error due to ensembl_id
  }
}
cat("\n")
######################
# Run DESeq2
######################

head(counts)
print(sample_info)
cat("Running DESeq2...",  sep="\n")
dds <- deAnalysis(counts, sample_info)
cat("\n")
#  Export PCA plot (default deseq PCA on 500 most variable genes).

cat("Generating PCA plots...",  sep="\n")
if (!dir.exists(fig_outdir)) { dir.create(fig_outdir) }

pdf(file.path(fig_outdir,"PCA.pdf"))

plotPCA(varianceStabilizingTransformation(dds), intgroup = "Treatment")

dev.off()
cat("\n")

cat("Generation contrasts...",  sep="\n")
colData(dds)
colData(dds)$Treatment

ctrl = sample_info %>%
  filter(Control == 1) %>%
  .$Treatment %>%
  unique

cases = sample_info %>%
  filter(Control == 0) %>%
  .$Treatment %>%
  unique

cat("\n")
i <- 1
for (case in cases) {
  
  if (!dir.exists(file.path(data_outdir,case))) { dir.create(file.path(data_outdir,case)) }

  ######################
  # Extract contrasts
  ######################
  cat(paste("Extracting contrast for ", ctrl, "_vs_", case, "...", sep=""), sep="\n")
  export.deseq2 <- getContrast(counts, dds, case, ctrl) # CHANGE NAME SIGNIF

  write_tsv(export.deseq2,
            file.path(data_outdir,case,"DESeq2.txt")
  )
  cat("\n")
  ######################
  # Plot results
  ######################
  cat("\tGenerating MA plot...",  sep="\n")

  pdf(file.path(fig_outdir,case,"MAPlot.pdf"))
  MAPlot(export.deseq2, case, ctrl, pval_th)
  dev.off()
  cat(paste(i, "/", length(cases), " done!", sep=""), sep="\n")
  cat("\n")
  i <- i+1
}

cat("\n")

cat("Saving output data...",  sep="\n")
saveRDS(dds, file = ddsRDS)
cat("\n")

cat("DONE!", sep="\n")