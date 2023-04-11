log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(biomaRt)
        library(dplyr)
})


tcounTSV <- snakemake@input[["tcounTSV"]]
genesetTSV <- snakemake@input[["genesetTSV"]]
transcriptsetTSV <- snakemake@input[["transcriptsetTSV"]]

geneCountsRData <- snakemake@output[["geneCountsRData"]]
geneCountsTSV <- snakemake@output[["geneCountsTSV"]]
geneCountsExtraTSV <- snakemake@output[["geneCountsExtraTSV"]]

ensembl_version <- snakemake@params[["ensembl_version"]]
genome_build <- snakemake@wildcards[["genome"]]

if (genome_build == "m38"){
        ensembl_data <- "mmusculus_gene_ensembl"
        gene_symbol<- "mgi_symbol"
} else if (genome_build == "hg38"){
        ensembl_data <- "hsapiens_gene_ensembl"
        gene_symbol <- "hgnc_symbol"
}



cat("Reading Ensembl tables...", sep="\n")
ensembl <- read.table(genesetTSV, header = TRUE, sep="\t")
gene_to_transcript <- read.table(transcriptsetTSV, header = TRUE, sep="\t")

cat("\n")

cat(paste("Writting collapsed count data for ", tcounTSV, "...", sep=""), sep="\n")
countsTable <- read.table(tcounTSV, stringsAsFactors=FALSE, header=TRUE, sep="\t")
collapsedCounts <- data.frame(ensembl_transcript_id = countsTable["Name"])
colnames(collapsedCounts) <- "ensembl_transcript_id"
collapsedCounts <- cbind(collapsedCounts, countsTable[,7:13])
cat("\n")


cat(paste("Original data contains", nrow(collapsedCounts), " transcript entries."), sep="\n")
removed <-  nrow(collapsedCounts)
collapsedCounts <- inner_join(gene_to_transcript, collapsedCounts, by="ensembl_transcript_id")
collapsedCounts$ensembl_transcript_id <- NULL
collapsedCounts <- aggregate(.~ensembl_gene_id, data = collapsedCounts, FUN = sum)
removed <-  removed -  nrow(collapsedCounts)
cat(paste("Collapsed all read counts data contains", nrow(collapsedCounts), " transcript entries. Removed", removed, "duplicated gene entries."), sep="\n")
cat("\n")


cat("Computing new conversion rates and cpm according to gene read counts...", sep="\n")
collapsedCounts <- collapsedCounts %>%
                    mutate(ConversionRateAdj = case_when(
                    CoverageOnTs > 0 ~ ConversionsOnTs/CoverageOnTs,
                    CoverageOnTs <= 0 ~ 0
                    ))

collapsedCounts$TcReadCPM <- (collapsedCounts$TcReadCount * 1000000)/ sum(collapsedCounts$TcReadCount)
collapsedCounts$ReadCountCPM <- (collapsedCounts$ReadCount * 1000000)/ sum(collapsedCounts$ReadCount)
cat("\n")

cat("Computing new conversion rates and cpm according to gene read counts...", sep="\n")
collapsedCounts$normReadCounts <- collapsedCounts$ReadCount/collapsedCounts$gene_length
collapsedCounts$normTCeadCounts <- collapsedCounts$TCReadCount/collapsedCounts$gene_length


cat("Adding extra gene information to counts extenden table...", sep="\n")

# ------- Genes of interest for further analysis -------
# To add your own category:
#       1. Create a list of gene names (Ensembl)
#       2. Use <collapsedCountsExtra[[gene_symbol]]  %in% [YOUR_NEW_LIST] ~ "[YOUR_CATEGORY_NAME]"> line to add your list
#       3. Add <YOUR_NEW_LIST> to the exception category "Other" evaluation

OSN <- c("Pou5f1", "Sox2" , "Nanog")
OSN_target <- c("Max", "Klf4", "Med1", "Rtf1", "Wdr5", "Ash2l", "Eed", "Jarid2", "Lin28")
hkg <- c("Actb", "Atp5f1" , "B2m", "Gapdh", "Hprt1", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc")


collapsedCountsExtra <- inner_join(collapsedCounts, ensembl, by="ensembl_gene_id", multiple = "all")
collapsedCountsExtra$GeneSize <- collapsedCountsExtra$end_position - collapsedCountsExtra$start_position
collapsedCountsExtra$GeneProfile <- case_when(collapsedCountsExtra[[gene_symbol]] %in% hkg ~ "Housekeeping",
                                        collapsedCountsExtra[[gene_symbol]]  %in% OSN_target ~ "OSN Target",
                                        collapsedCountsExtra[[gene_symbol]]  %in% OSN ~ "Core OSN",
                                        !(collapsedCountsExtra[[gene_symbol]]  %in% c(hkg, OSN_target, OSN)) ~ "Other"
                                        )

collapsedCountsExtra$TcReadRange <- cut(collapsedCounts$TcReadCount, c(0,1,2,3,6,11,26,50,101,251,501,1001, 5001), 
                                        labels=c("0", "1", "2", "3-5", "6-10", "11-25", "26-50", "51-100", "101-250", "251-500", "501-1000", ">1000"),
                                        right = FALSE) # FALSE -> [x)   TRUE -> (x]        
cat("\n")

write.table(collapsedCounts, file=geneCountsTSV, sep="\t",quote=FALSE, row.names=FALSE)
rownames(collapsedCountsExtra) <- collapsedCountsExtra$ensembl_gene_id
collapsedCountsExtra$ensembl_gene_id <- NULL

write.table(collapsedCountsExtra, file=geneCountsExtraTSV, sep="\t",quote=FALSE, row.names=FALSE)
rownames(collapsedCounts) <- collapsedCounts$ensembl_gene_id
collapsedCounts$ensembl_gene_id <- NULL

cat("Saving output data...",  sep="\n")
save(collapsedCounts, collapsedCountsExtra, file = geneCountsRData)
cat("\n")

cat("DONE!", sep="\n")