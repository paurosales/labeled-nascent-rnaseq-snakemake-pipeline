log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        # library(biomaRt)
        library(dplyr)
})


tcounTSV <- snakemake@input[["tcounTSV"]]
transcriptsBED <- snakemake@input[[1]]
genesBED <- snakemake@input[[2]]

genecountsRData <- snakemake@output[["genecountsRData"]]
genecountsTSV <- snakemake@output[["genecountsTSV"]]
genecountsExtTSV <- snakemake@output[["genecountsExtTSV"]]

cat("Formating BED files...", sep="\n")
tx_info <- read.table(transcriptsBED, header = FALSE, sep="\t")
names(tx_info) <- c("Chromosome", "Start", "End", "TranscriptID", "Score", "Strand", "GeneID") # names must correspond to BED file colums!
tx_info <- tx_info[,c("TranscriptID", "GeneID")]

gene_info <- read.table(genesBED, header = FALSE, sep="\t")
names(gene_info) <- c("Chromosome", "Start", "End", "GeneID", "Score", "Strand", "GeneSize", "GeneName", "GeneType") # names must correspond to BED file colums!
gene_info <- gene_info[,c("GeneID", "GeneSize", "GeneName", "GeneType")]

norm_factor <- gene_info[,c("GeneID", "GeneSize")]
cat("\n")

cat(paste("Writting collapsed count data for ", tcounTSV, "...", sep=""), sep="\n")
counts_tab <- read.table(tcounTSV, stringsAsFactors=FALSE, header=TRUE, sep="\t")
genecounts <- data.frame(TranscriptID = counts_tab["Name"])
names(genecounts) <- "TranscriptID"
genecounts <- cbind(genecounts, counts_tab[,7:13])

cat("\n")


cat(paste("Original data contains", nrow(genecounts), " transcript entries."), sep="\n")
removed <-  nrow(genecounts)
genecounts <- inner_join(tx_info, genecounts, by="TranscriptID")
genecounts$TranscriptID <- NULL

genecounts <- aggregate(.~GeneID, data = genecounts, FUN = sum)
removed <-  removed -  nrow(genecounts)
cat(paste("Collapsed all read counts data contains", nrow(genecounts), " transcript entries. Removed", removed, "duplicated gene entries."), sep="\n")
cat("\n")




cat("Computing new conversion rates and cpm according to gene read counts...", sep="\n")
genecounts <- genecounts %>%
                    mutate(ConversionRate = case_when(
                    CoverageOnTs > 0 ~ ConversionsOnTs/CoverageOnTs,
                    CoverageOnTs <= 0 ~ 0
                    ))

genecounts$TcReadCPM <- (genecounts$TcReadCount * 1000000)/ sum(genecounts$TcReadCount)
genecounts$ReadCountCPM <- (genecounts$ReadCount * 1000000)/ sum(genecounts$ReadCount)


cat("Normalizing counts by gene length...", sep="\n")
genecounts <- inner_join(norm_factor, genecounts, by="GeneID")

genecounts$NormReadCount <- genecounts$ReadCount/genecounts$GeneSize
genecounts$NormTcReadCount <- genecounts$TcReadCount/genecounts$GeneSize

genecounts$NormTcReadCPM <- (genecounts$TcReadCount * 1000000)/ sum(genecounts$TcReadCount)
genecounts$NormReadCountCPM <- (genecounts$ReadCount * 1000000)/ sum(genecounts$ReadCount)
genecounts$GeneSize <- NULL

cat("\n")



cat("Adding extra gene information to counts extended table...", sep="\n")

# ------- Genes of interest for further analysis -------
# To add your own category:
#       1. Create a list of gene names (Ensembl)
#       2. Use <genecounts_ext[[gene_symbol]]  %in% [YOUR_NEW_LIST] ~ "[YOUR_CATEGORY_NAME]"> line to add your list
#       3. Add <YOUR_NEW_LIST> to the exception category "Other" evaluation

OSN <- c("Pou5f1", "Sox2" , "Nanog")
OSN_target <- c("Max", "Klf4", "Med1", "Rtf1", "Wdr5", "Ash2l", "Eed", "Jarid2", "Lin28")
hkg <- c("Actb", "Atp5f1" , "B2m", "Gapdh", "Hprt1", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc")




genecounts_ext <- inner_join(genecounts, gene_info, by="GeneID", multiple = "all")

genecounts_ext$GeneProfile <- case_when(genecounts_ext$GeneName %in% hkg ~ "Housekeeping",
                                        genecounts_ext$GeneName  %in% OSN_target ~ "OSN Target",
                                        genecounts_ext$GeneName  %in% OSN ~ "Core OSN",
                                        !(genecounts_ext$GeneName  %in% c(hkg, OSN_target, OSN)) ~ "Other"
                                        )


genecounts_ext$TcReadRange <- cut(genecounts$TcReadCount, c(0,1,2,3,6,11,26,50,101,251,501,1001, 5001), 
                                        labels=c("0", "1", "2", "3-5", "6-10", "11-25", "26-50", "51-100", "101-250", "251-500", "501-1000", ">1000"),
                                        right = FALSE) # FALSE -> [x)   TRUE -> (x]  

cat("\n")

write.table(genecounts, file=genecountsTSV, sep="\t",quote=FALSE, row.names=FALSE)

write.table(genecounts_ext, file=genecountsExtTSV, sep="\t",quote=FALSE, row.names=FALSE)

cat("Saving output data...",  sep="\n")
save(genecounts, genecounts_ext, file = genecountsRData)
cat("\n")

cat("DONE!", sep="\n")