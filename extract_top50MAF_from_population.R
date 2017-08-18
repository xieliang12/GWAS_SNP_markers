factor_to_numeric <- function(x) {
  as.numeric(levels(x))[as.integer(x)]
}
library(data.table)
allele_freqs <- read.csv("data/hapmap/allele_freqs_chrY_CHB_r28_nr.b36_fwd.txt",header=TRUE, row.names = NULL, sep=" ")
  selected_columns <- c("rs.", "chrom", "pos", "strand", "build", "refallele","refallele_freq",
                      "otherallele", "otherallele_freq")
  allele_freqs <- allele_freqs[, selected_columns]
  allele_freqs$refallele_freq <- factor_to_numeric(allele_freqs$refallele_freq)
  allele_freqs$otherallele_freq <- factor_to_numeric(allele_freqs$otherallele_freq)
  allele_freqs$MAF <- ifelse(allele_freqs$refallele_freq<allele_freqs$otherallele_freq, 
                           allele_freqs$refallele_freq, allele_freqs$otherallele_freq)
  allele_freqs$alleles <- paste0(allele_freqs$refallele, "/", allele_freqs$otherallele, sep="")
  allele_freqs$chrom <- as.character(allele_freqs$chrom)
  allele_freqs <- subset(allele_freqs, allele_freqs$chrom != "chrom")
  allele_freqs_order_by_MAF_unique <- allele_freqs[!duplicated(allele_freqs[c("chrom","MAF")]),]
  allele_freqs_order_by_MAF_unique <- allele_freqs_order_by_MAF_unique[rev(order(allele_freqs_order_by_MAF_unique$MAF)),]
  top2MAF_by_chrom <- data.table(allele_freqs_order_by_MAF_unique, key="chrom")
  top2MAF_all_chrom <- top2MAF_by_chrom[, head(.SD, 200), by=chrom]
  write.table(top2MAF_all_chrom, "data/SNP_first_round/first_panel_control_markers_ChrY_top200.txt", quote=FALSE, row.names=FALSE, sep="\t")


 
control_markers_chrom <- read.table("Downloads/data_for_phenogram_plotting.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
control_markers_chrom <- control_markers_chrom[,c(5,2)]
control_markers_first_panel <- top2MAF_all_chrom[(top2MAF_all_chrom$rs. %in% control_markers_chrom$ANNOTATION),]

write.table(contorl_markers_first_panel, "data/SNP_first_round/first_panel_control_markers.txt", quote=FALSE, row.names=FALSE, sep="\t")

