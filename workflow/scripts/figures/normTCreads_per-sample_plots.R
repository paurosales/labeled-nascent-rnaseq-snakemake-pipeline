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
geneCountsRData <- snakemake@input[["geneCountsRData"]]

tcFreqPDF <- snakemake@output[["tcFreqPDF"]]
tcFreqNoZeroPDF <- snakemake@output[["tcFreqNoZeroPDF"]]
exprCountsPDF <- snakemake@output[["exprCountsPDF"]]
exprCountsNoLimsPDF <- snakemake@output[["exprCountsNoLimsPDF"]]
exprCpmPDF  <- snakemake@output[["exprCpmPDF"]]

genome_build <- snakemake@params[["genome_build"]]


cat("Formating input data...", sep="\n")
load(geneCountsRData) #  collapsedCounts, collapsedCountsExtra

collapsedCountsExtra <- collapsedCountsExtra[order(collapsedCountsExtra$GeneProfile, decreasing=TRUE),]
palette <- brewer.pal(n = 8, name = "Dark2")

if (genome_build == "m39"){
        gene_symbol<- "mgi_symbol"
} else if (genome_build == "hg38"){
        gene_symbol <- "hgnc_symbol"
}

# Asign transparency for scatter plots
collapsedCountsExtra$alpha <- 1
collapsedCountsExtra$alpha[collapsedCountsExtra$GeneProfile == "Other"] <- 1/10 

topUpReads <- collapsedCountsExtra[order(collapsedCountsExtra$normReadCount, decreasing=TRUE),][[gene_symbol]][1:5]
topUpTC <- collapsedCountsExtra[order(collapsedCountsExtra$normTcReadCount, decreasing=TRUE),][[gene_symbol]][1:5]

cat("Plotting T>C reads frequencies...", sep="\n")
p_freq <- ggplot(collapsedCountsExtra, aes(normTcReadRange)) + # p values for genes with mean normalized count larger than 1
        geom_bar(fill = palette[2]) +
        labs(x = "# TC Reads", y = "Frequency",  tag = paste("n =", nrow(collapsedCountsExtra))) +
        scale_y_continuous(limits=c(0,50000), label=comma) +
        theme_classic() +
        theme(plot.tag.position = c(0.5, 1))
         
pdf(file=tcFreqPDF, width=8, height=6)
    print(p_freq)     
dev.off() 
cat("\n")

p_freq <- ggplot(collapsedCountsExtra[collapsedCountsExtra$normTcReadRange!="0",], aes(TcReadRange)) + # p values for genes with mean normalized count larger than 1
        geom_bar(fill = palette[2]) +
        labs(x = "# TC Reads", y = "Frequency",  tag = paste("n =", nrow(collapsedCounts[collapsedCountsExtra$normTcReadRange!="0",]))) +
        scale_y_continuous(limits=c(0,3000),label=comma) +
        theme_classic() +
        theme(plot.tag.position = c(0.5, 1))
         
pdf(file=tcFreqNoZeroPDF, width=8, height=6)
    print(p_freq)     
dev.off() 
cat("\n")

cat("Plotting T>C reads vs. Steady state...", sep="\n")
# collapsedCountsExtra <- filter(collapsedCountsExtra, normReadCount < 100000 ) # remove outlier
p_counts <- ggplot(collapsedCountsExtra, aes(x = normReadCount, y = normTcReadCount)) +
                geom_point(aes(color = GeneProfile, alpha = alpha)) +
                geom_text_repel(aes(label = dplyr::case_when(collapsedCountsExtra[[gene_symbol]] %in% topUpReads ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra[[gene_symbol]] %in% topUpTC ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra$GeneProfile == "Housekeeping" ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra$GeneProfile == "Core OSN" ~ collapsedCountsExtra[[gene_symbol]])), size=2, max.overlaps = 50) +
                guides(alpha = FALSE) +
                labs(x = "# Overall Reads", y = "# T>C Reads", color = "Gene Profile\n") +
                scale_color_manual(values=c("red", "lightblue", "blue", "#999999"))+
                scale_x_continuous(limits=c(0,150000),label=comma) +
                scale_y_continuous(limits=c(0, 5000), label=comma) + # scale limits
                theme_bw()

pdf(file=exprCountsPDF, width=8, height=6)
    print(p_counts)     
dev.off() 
cat("\n")

cat("Plotting T>C reads vs. Steady state...", sep="\n")
# collapsedCountsExtra <- filter(collapsedCountsExtra, normReadCount < 100000 ) # remove outlier
p_counts <- ggplot(collapsedCountsExtra, aes(x = normReadCount, y = normTcReadCount)) +
                geom_point(aes(color = GeneProfile, alpha = alpha)) +
                geom_text_repel(aes(label = dplyr::case_when(collapsedCountsExtra[[gene_symbol]] %in% topUpReads ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra[[gene_symbol]] %in% topUpTC ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra$GeneProfile == "Housekeeping" ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra$GeneProfile == "Core OSN" ~ collapsedCountsExtra[[gene_symbol]])), size=2, max.overlaps = 50) +
                guides(alpha = FALSE) +
                labs(x = "# Overall Reads", y = "# T>C Reads", color = "Gene Profile\n") +
                scale_color_manual(values=c("red", "lightblue", "blue", "#999999"))+
                scale_x_continuous(label=comma) +
                theme_bw()

pdf(file=exprCountsNoLimsPDF, width=8, height=6)
    print(p_counts)     
dev.off() 
cat("\n")


cat("Plotting CPM...", sep="\n")

p_CPM <- ggplot(collapsedCountsExtra, aes(x = normReadCountCPM, y = normTcReadCPM)) +
                geom_point(aes(color = GeneProfile, alpha = alpha)) +
                geom_text_repel(aes(label = dplyr::case_when(collapsedCountsExtra$GeneProfile == "Housekeeping" ~ collapsedCountsExtra[[gene_symbol]],
                                        collapsedCountsExtra$GeneProfile == "Core OSN" ~ collapsedCountsExtra[[gene_symbol]])), size=2, max.overlaps = 50) +
                guides(alpha = FALSE) +
                labs(x = "Steady state (cpm)", y = "T>C Reads (cpm)", color = "Gene Profile\n") +
                scale_color_manual(values=c("red", "lightblue", "blue", "#999999")) +
                scale_x_log10(limits = c(10^-1, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + # collapse log steps to same length
                scale_y_log10(limits = c(10^-2, 10^4), breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + # collapse log steps to same length
                theme_bw()

pdf(file=exprCpmPDF, width=8, height=6)
    print(p_CPM)     
dev.off() 
cat("\n")
   

cat("DONE!")
cat("\n")