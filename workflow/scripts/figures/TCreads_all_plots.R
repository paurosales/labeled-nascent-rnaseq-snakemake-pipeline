log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(RColorBrewer)
        library(tibble)
        library(scales)
        library(ggrepel)
})

# Parse snakemake arguments
files <- snakemake@input[["geneCounts"]]
tcPercentPDF <- snakemake@output[["tcPercentPDF"]]
sample_manifest <- snakemake@params[["sample_manifest"]]
# indir <- snakemake@params[["indir"]]
genome_build <- snakemake@params[["genome_build"]]

palette <- brewer.pal(n = 12, name = "Paired")

if (genome_build == "m38"){
        gene_symbol<- "mgi_symbol"
} else if (genome_build == "hg38"){
        gene_symbol <- "hgnc_symbol"
}


samples <- read.table(sample_manifest, stringsAsFactors=FALSE, header=T, sep="\t")

# files <- paste(indir, "/", samples$Label, "_tcount_by_gene_extended.tsv", sep="") # particular handle tcount.tsv files, change if needed
order <- c("0", "1", "2", "3-5", "6-10", "11-25", "26-50", "51-100", "101-250", "251-500", "501-1000", ">1000")
cat(paste("Plotting T>C read counts percentages for ", length(files), " samples with...", sep=""), sep="\n")
allLabels <- c()
allRanges <- c()
allFreq <- c()
for (f in 1:length(files)){
        collapsedCountsExtra <- read.table(files[f], stringsAsFactors=FALSE, header=TRUE, sep="\t")
        rangeFreq <- as.data.frame(table(collapsedCountsExtra$TcReadRange))
        rangeFreq <- rangeFreq[match(order, rangeFreq$Var1),]
        # label <- rep(paste(samples$Condition[f], "_", samples$Target[f], sep=""), nrow(rangeFreq))
        label <- rep(paste(samples$Sample_type[f], "_", samples$Treatment[f], "_Bio-rep_", samples$Bio_rep[f], sep=""), nrow(rangeFreq))
        allLabels <- append(allLabels, label)
        allRanges <- append(allRanges, rangeFreq$Var1)
        allFreq <- append(allFreq, rangeFreq$Freq)
}


rangesData <- data.frame(sample=allLabels, ranges=allRanges, freq=allFreq)

p_percent <- ggplot(rangesData, aes(fill= factor(ranges, levels=order), y=freq, x=sample)) + 
                geom_bar(position="fill", stat="identity") + 
                scale_fill_manual(values=palette) +
                scale_y_continuous(labels=percent_format(scale = 100)) +
                labs(x = "Sample", y = "Percentage", fill = "Counts\n", tag = paste("n =", nrow(collapsedCountsExtra))) +
                # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
                # scale_x_discrete(guide = guide_axis(angle = 90)) +
                coord_flip() +
                theme_classic() +
                theme(plot.tag.position = c(0.95, 0.1))

pdf(file=tcPercentPDF, width=8, height=6)
    print(p_percent)     
dev.off() 
cat("\n")


