contributes <- read.csv("ruby/genowise_simplicity/db/contribute_table_postgresql.csv",
                        header=TRUE, sep="\t")
contributes <- contributes[which(contributes$Rs_id!="Pending"),]
contributes$P_value <- factor_to_numeric(contributes$P_value)
contributes$P_value <- ifelse(is.na(contributes$P_value), -1, contributes$P_value)
contributes$OR <- factor_to_numeric(contributes$OR)
contributes$OR <- ifelse(contributes$OR=="NA", -1, contributes$OR)
write.table(contributes, "ruby/genowise_simplicity/db/contribute_table_postgresql.csv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

markers <- read.csv("ruby/genowise_simplicity/db/markers_table.csv", header=FALSE, 
                    stringsAsFactors=FALSE, sep="\t")
colnames(markers) <- markers[1,]
markers <- markers[-1,]
markers <- markers[, c(1:16)]
f2 <- function(x) {x <- as.numeric(x)}
cols <- c("CHBSX","JPT","CEU","EUR","EAS","AMR","SAS")
markers[, cols] <- lapply(markers[,cols], f2)
replacement <- function(x) ifelse(is.na(x), -1, x)
markers[, cols] <- lapply(markers[,cols], replacement)
write.table(markers, "ruby/genowise_simplicity/db/markers_table.csv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

hpo <- read.csv("/Volumes/KINGSTON/201409to201412/database/hpo_merge_postgresql.tsv",
                           sep="\t", stringsAsFactors=FALSE, header=FALSE)
colnames(hpo) <- hpo[1,]
hpo <- hpo[-1, ]
Rs_Hpo <- hpo[, c(2,6)]
Rs_Hpo <- Rs_Hpo[!duplicated(Rs_Hpo),]
colnames(Rs_Hpo) <- c("Rs_id","HPO_id")
hpo <- hpo[, -2]
hpo <- hpo[!duplicated(hpo),]
colnames(hpo) <- c("Hpo_Term","DiseaseID", "Gene_symbol","Gene_id","Hpo_id")
write.table(Rs_Hpo, "/Volumes/KINGSTON/201409to201412/database/Rs_Hpo_postgresql.csv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(hpo, "/Volumes/KINGSTON/201409to201412/database/hpo_postgresql.csv",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")