suppressPackageStartupMessages({
        library(biomaRt)
        library(dplyr)
})



genesetTSV <- snakemake@output[["genesetTSV"]]
transcriptsetTSV <- snakemake@output[["transcriptsetTSV"]]

ensembl_version <- snakemake@params[["ensembl_version"]]
genome_build <- snakemake@wildcards[["genome"]]

if (genome_build == "m38"){
        ensembl_data <- "mmusculus_gene_ensembl"
        gene_symbol<- "mgi_symbol"
} else if (genome_build == "hg38"){
        ensembl_data <- "hsapiens_gene_ensembl"
        gene_symbol <- "hgnc_symbol"
}

# biomartCacheClear()
cat("Fetching data from Ensembl...", sep="\n")
mart <- useDataset(ensembl_data, useEnsembl(biomart="ensembl", version=ensembl_version))
geneset <- getBM(attributes=c("ensembl_gene_id", gene_symbol, "chromosome_name", 
                                "start_position", "end_position", "gene_biotype"), mart = mart)
geneset <- geneset[!duplicated(geneset$ensembl_gene_id),]
geneset$gene_length <- geneset$end_position - geneset$start_position
write.table(as.data.frame(geneset), file=genesetTSV, sep="\t", quote=FALSE, row.names=FALSE)

gene_to_transcript <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)

#         cat("Filtering non-standard chromosomes and lncRNA genes out...", sep="\n")
#         ensembl <- ensembl[grep("CHR|GL|KI", ensembl$chromosome_name, invert=T),]
#         ensembl <- filter(ensembl, mgi_symbol != "", gene_biotype != "lncRNA") # exclude empty gene symbols and lncRNA genes
        # cat("\n")

write.table(as.data.frame(gene_to_transcript), file=transcriptsetTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")