# ADAPTED FROM SLAMDUNK CODE
#
# This script is an adaptation from : globalRatePlotter.R script
# defined in Slamdunk stats.py (https://github.com/t-neumann/slamdunk/blob/master/slamdunk/dunks/stats.py#L668)
#
# Adaptations: Added more columns regarding gene information

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


raw_files <- snakemake@input[["raw_cts"]]
collapsed_files <- snakemake@input[["collapsed_cts"]]
geneset_TSV <- snakemake@input[["geneset_TSV"]]
evalExpression <- snakemake@params[["eval_expression"]]

outputFile <- snakemake@output[[1]]

columnName <- 2

# collapsed_files = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))



readMeatInfo <- function(fileName) {
  #fileName = collapsed_files[1]
  sampleInfo = read.table(fileName, nrows = 1, comment.char = "")
  version = paste(lapply(sampleInfo[1,1:3], as.character), collapse = '\t')
  sampleID = as.character(sampleInfo[1, ]$V7)
  sampleName = as.character(sampleInfo[1, ]$V6)
  sampleType = as.character(sampleInfo[1, ]$V8)
  sampleTime = as.numeric(sampleInfo[1, ]$V9)
  sampleInfo = read.table(fileName, nrows = 1, skip = 1, comment.char = "")
  annotationMD5 = as.character(sampleInfo[1, ]$V3)
  annotationName = as.character(sampleInfo[1, ]$V2)
  c(sampleID, sampleName, sampleType, sampleTime, annotationName, annotationMD5, version)
}

sampleNumber = length(collapsed_files)
mergedRates = data.frame()

annotationName = ""
annotationMD5 = ""
version = ""
IDs = c()

geneset <- read.table(geneset_TSV, header = FALSE,  sep = '\t',  stringsAsFactors = FALSE)
names(geneset) <- c("Chromosome", "Start", "End", "Name", "Strand", "Length", "Symbol", "Biotype")


# Merge rates from all samples
for(i in 1:length(collapsed_files)) {
  #i = 1
  file = collapsed_files[i]
  meta = readMeatInfo(raw_files[i])
  sampleName = meta[columnName]

  if(i == 1) {
    version = meta[7]
    annotationName = meta[5]
    annotationMD5 = meta[6]
  } else {
    if(annotationMD5 != meta[6]) {

    }
  }

  IDs = c(IDs, as.numeric(meta[1]))
  data = read.table(file, header = T)
  names(data) <- c("Name", "Length", "ReadsCPM", "ConversionRate", "Tcontent", "CoverageOnTs", "ConversionsOnTs", "ReadCount", "TcReadCount", "multimapCount")
    gene_info <- geneset[match(data$Name, geneset$Name),]
    head(gene_info)
    print(nrow(gene_info))
    head(data)
    print(nrow(data))
  if(i == 1) {
    # mergedRates = data[, c(1:6)]
    mergedRates = gene_info
    mergedRates$avgReadsCPM = data$ReadsCPM
    mergedRates$avgMultimapper = data$multimapCount
    mergedRates$avgTcontent = data$Tcontent
    mergedRates$avgCoverageOnTs = data$CoverageOnTs
  } else {
    mergedRates$avgReadsCPM = mergedRates$avgReadsCPM + data$ReadsCPM
    mergedRates$avgMultimapper = mergedRates$avgMultimapper + data$multimapCount
    mergedRates$avgTcontent = mergedRates$avgTcontent + data$Tcontent
    mergedRates$avgCoverageOnTs = mergedRates$avgCoverageOnTs + data$CoverageOnTs
  }
  #if(perRead == T) {
  attach(data)
  #mergedRates[,sampleName] = data$TcReadCount / data$ReadCount
  mergedRates[,sampleName] = eval(parse(text=evalExpression))
  detach(data)
  mergedRates[data$ReadCount == 0,sampleName] = 0
  #} else {
  #  mergedRates[,sampleName] = data$ConversionRate
  #}
}
# compute average CPM and multimapper per UTR
mergedRates$avgReadsCPM = mergedRates$avgReadsCPM / sampleNumber
mergedRates$avgMultimapper = mergedRates$avgMultimapper / sampleNumber
mergedRates$avgTcontent = mergedRates$avgTcontent / sampleNumber
mergedRates$avgCoverageOnTs = mergedRates$avgCoverageOnTs / sampleNumber

#head(mergedRates)
# Sort columns by sample name
colNumber = length(colnames(mergedRates))
firstSampleColumn = (colNumber - sampleNumber + 1)
sampleNames = colnames(mergedRates)[firstSampleColumn:colNumber]
sampleColumnOrder = order(IDs)
mergedRates = mergedRates[, c(1:(firstSampleColumn - 1), (sampleColumnOrder + firstSampleColumn - 1))]

#head(mergedRates)

# Write to output file
con <- file(outputFile, open="wt")
writeLines(version, con)
writeLines(paste0("#Annotation:\t", annotationName, "\t", annotationMD5), con)
writeLines(paste0("#Expression:\t", evalExpression), con)
write.table(mergedRates, con, sep = "\t", quote = F, row.names = F, col.names = T)
close(con)