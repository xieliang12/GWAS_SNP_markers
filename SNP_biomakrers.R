###handling with the output files downloaded from HapMap Biomart tools by limiting the population###
###to Chinese and Caucasian after limiting the marker list to the uploaded file###
chinese <- read.table("data/SNP_first_round/Chinese_HapMart_export.txt", header=TRUE, sep=",")
caucasian <- read.table("data/SNP_first_round/caucasian_HapMart_export.txt", header=TRUE, sep=",")
colnames(chinese)[5:6] <- c("chinese_reference_freq", "chinese_other_freq")
colnames(caucasian)[7:8] <- c("caucasian_other_freq","caucasian_reference_freq")

joined <- merge(chinese, caucasian, by="marker.id", all=TRUE)
joined$chinese_MAF <- ifelse(joined$chinese_reference_freq>joined$chinese_other_freq, joined$chinese_other_freq, joined$chinese_reference_freq)
joined$chinese_MAF_allele <- ifelse(joined$chinese_reference_freq>joined$chinese_other_freq, as.character(joined$other.allele.x), as.character(joined$reference.allele.x))
joined$caucasian_MAF <- ifelse(joined$caucasian_reference_freq>joined$caucasian_other_freq, joined$caucasian_other_freq, joined$caucasian_reference_freq)
joined$caucasian_MAF_allele <- ifelse(joined$caucasian_reference_freq>joined$caucasian_other_freq, as.character(joined$other.allele.y), as.character(joined$reference.allele.y))

###make a annotation file for uploading onto ENSEMBL Variant Effect Predictor tools for getting###
###the annotation for SNP###
select_columns <- c("marker.id","chromosome.x", "position.x", "alleles.x", "chinese_MAF","chinese_MAF_allele", "caucasian_MAF", "caucasian_MAF_allele", "strand.y")
SNP_list <- joined[, select_columns]
SNP_list$chromosome.x <- gsub("^chr", "", SNP_list$chromosome.x)
VEP_annotation_file <- SNP_list[,c(2,3,3,4)]
write.table(VEP_annotation_file, "data/SNP_first_round/VEP_annotation_file.txt", quote=FALSE, row.names = FALSE, sep="\t")

###clean the file content from the output of ENSEMBL Variant Effect Predictor tools###
exon_region <- read.table("data/SNP_first_round/VEP_annotation_for_SNP.txt", fill=TRUE, sep="\t")
SNP_annotation <- exon_region[, c("V2","V4","V5","V6","V7","V13","V18")]
SNP_annotation$V13 <- gsub(",", "", SNP_annotation$V13)
SNP_annotation$V13 <- gsub("CM[0-9]+", "", SNP_annotation$V13)
SNP_annotation$V13 <- gsub("COSM[0-9]+", "", SNP_annotation$V13)
SNP_annotation$V13 <- gsub("C[SR][0-9]+", "", SNP_annotation$V13)

###make a bed file for genome build GRCh38 for extracting flanking sequence### 
colnames(SNP_annotation) <- c("GRCh38_coordinate","Gene", "Feature", "Feature_type", "Consequence","marker.id","HGNC")
SNP_information <- merge(SNP_list, SNP_annotation, by="marker.id", all=TRUE)
SNP_information$chromStart <- SNP_information$position.x - 151
SNP_information$chromEnd <- SNP_information$position.x + 150
SNP_information$chromosome <- paste0("chr", SNP_information$chromosome.x)
SNP_for_build36_bed <- SNP_information[,c("chromosome", "chromStart", "chromEnd")]
write.table(SNP_for_build36_bed, "data/SNP_first_round/SNP_for_build36_bed.txt", quote=FALSE, row.names=FALSE, sep="\t")

###extract exactly 150bp sequence data upstream and downstream of SNP location###
library(stringr)
SNP_list_upto150 <- read.table("data/SNP_first_round/SNP_first_round_150.tab", sep="\t")
SNP_list_upto150$V1 <- as.character(SNP_list_upto150$V1)
SNP_list_upto150$V2 <- as.character(SNP_list_upto150$V2)
for (i in 1:nrow(SNP_list_upto150)) {
  SNP_list_upto150$start[i] <- substr(SNP_list_upto150[i,1], (unlist(str_locate_all(pattern=":", SNP_list_upto150[i,1]))[1]+1), 
                                 (unlist(str_locate_all(pattern="-", SNP_list_upto150[i,1]))[1]-1))
}
SNP_list_upto150$start <- as.numeric(SNP_list_upto150$start)
SNP_list_upto150$first_150 <- toupper(substr(SNP_list_upto150$V2,1,150))
SNP_list_upto150$end_150 <- toupper(substr(SNP_list_upto150$V2,152,301))
SNP_list_upto150 <- SNP_list_upto150[,c("start","first_150","end_150")]
combined$SNP_upto150 <- paste0(combined$first_150, "[", combined$alleles.x, "]", combined$end_150)
final_SNP_information <- combined[,-c(1, 3, 17,19,20)]
colnames(final_SNP_information)[c(1,2,3,8)] <- c("rs_id","Cooridinate","alleles","strand")
final_SNP_information$NCBI_build <- "Build36.1"
write.table(final_SNP_information, "data/SNP_first_round/SNP_first_round.txt", quote=FALSE, row.names=FALSE, sep="\t")

###input a bed file after dealt with Awk for substracting and adding 150bp from the SNP###
###coordinate###
SNP_information <- read.table("data/SNP_first_round/third_round_control_markers_dbSNP_extract_sequence.bed", sep="\t")
colnames(SNP_information) <- c("chr","start","end","rs", "score","strand")
combined <- merge(SNP_information, SNP_list_upto150, by.x="start", by.y="start")
###input GRCh37 bed file get from the ouput file using NCBI assembly converting tools###
SNP_information_37 <- read.table("data/SNP_first_round/remapped_third_round_control_markers_GRCh37.bed", sep="\t")
colnames(SNP_information_37) <- c("chr","coordinate","coordiante_2","rs","score", "strand")
combined <- merge(combined, SNP_information_37, by.x="rs", by.y="rs", all=FALSE)
combined <- combined[,c("rs", "chr.x", "coordiante_2", "flanking_sequence_150bp")]
colnames(combined)[1:3] <- c("rs","chrom","coordinate")
combined$genome_build_version <- "GRCh37(hg19)"

###subset the SNP sequence dataset based on rs found in output file with Allele information###
###and get sequence output like this format ---ATCGTGGGG[G/C]GTCGTGCG---###
top2MAF_all_chrom <- as.data.frame(top2MAF_all_chrom)
top2MAF_all_chrom_clean <- subset(top2MAF_all_chrom, select=c("rs.", "alleles"))
colnames(top2MAF_all_chrom_clean) <- c("rs", "alleles")
chrom_with_alleles <- rbind(top2MAF_all_chrom_clean, obesity)
chrom_with_alleles2 <- merge(combined, chrom_with_alleles, by="rs")
chrom_with_alleles2$upstream150 <- substr(chrom_with_alleles2$flanking_sequence_300bp, 1, 150)
chrom_with_alleles2$downstream150 <- substr(chrom_with_alleles2$flanking_sequence_300bp, 152, 301)
chrom_with_alleles2$flanking_sequence_300bp <- paste0(chrom_with_alleles2$upstream150, 
                                                      "[", chrom_with_alleles2$alleles, "]", chrom_with_alleles2$downstream150)
chrom_with_alleles2 <- chrom_with_alleles2[,c(1:5)]
combined_without_SNPseq <- subset(combined, combined$sh==0)
combined_without_SNPseq <- combined_without_SNPseq[, c(1:5)]
third_round_control_markers <- rbind(chrom_with_alleles2, combined_without_SNPseq)
write.table(third_round_control_markers, "data/SNP_first_round/third_round_control_markers.txt", quote=FALSE, row.names=FALSE, sep="\t")