snps <- read.table("data/SNP_first_round/remapped_SNP_coordinate_hg19.bed", sep="\t", header=F)
colnames(snps) <- c("chr","start","end","id","score","strand")

library(ggplot2)
snpDensity <- ggplot(snps)+geom_histogram(aes(x=start), binwidth=1e7) +
  facet_wrap(~ chr, ncol=2) + ggtitle("Density of SNP327 across hg19") + 
  xlab("Position in the genome") + ylab("SNP density")+theme_bw()

png("snp327_density.png", 500, 750)
print(snpDensity)

snps <- snps[with(snps, order(chr, start)),]
snps$substract_value <- c(NA, snps$end[1:326])
snps$distribution <- NA
for (i in 2:nrow(snps)) {
  if (snps$chr[i-1] == snps$chr[i]) {
    snps$distribution[i] <- snps$start[i]-snps$substract_value[i]
  }
}

snps <- snps[, -c(7)]
w <- which(snps$distribution<150)
selected <-c(w[1]-1, w[1:2], w[3]-1, w[3:7])
snps <- snps[selected,]
snps$distribution[c(1,4)] <- ''
write.table(snps, "data/SNP_first_round/SNP_1st_round_distribution_lessthan150bp", quote=FALSE, row.names=FALSE,sep="\t")
