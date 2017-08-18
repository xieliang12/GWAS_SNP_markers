gwas_snp_trait <- read.csv("/Volumes/KINGSTON/201409to201412/database/gwasdb_20140908_snp_trait_filtered.tsv", 
                    header=TRUE, sep="\t", fill=TRUE)
hpo_term <- gwas_snp_trait[,c("SNPID","HPO_TERM")]

##use splitstackshape package to split HPO_TERM column into multiple HPO terms separated by
"|"##
library(splitstackshape)
hpo_term <- concat.split(hpo_term, "HPO_TERM", sep ="|", drop=TRUE)
library(reshape2)
hpo_term_long <- melt(hpo_term, id.vars="SNPID", value.name="HPO_TERM")
hpo_term_long <- hpo_term_long[,c(1,3)]
hpo_term_long <- hpo_term_long[!duplicated(hpo_term_long),]

##load humna phenotypes data to match with dataset hpo_term_long by HPO_TERM##
hpo_all_source <- read.csv("/Volumes//KINGSTON/201409to201412/database/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt",
                           sep="\t", header=FALSE)
colnames(hpo_all_source) <- c("diseaseID","gene_symbol","gene_id", "HPO_ID","HPO_TERM")
hpo_term_long$HPO_TERM <- tolower(hpo_term_long$HPO_TERM)
hpo_all_source$HPO_TERM <- tolower(hpo_all_source$HPO_TERM)
hpo_merge <- merge(hpo_term_long, hpo_all_source, by="HPO_TERM", all.x=TRUE)
hpo_merge <-  hpo_merge[!is.na(hpo_merge$HPO_ID),]
hpo_merge <- hpo_merge[!duplicated(hpo_merge),]
write.table(hpo_merge, "/Volumes//KINGSTON/201409to201412/database/hpo_merge_postgresql.tsv",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

##load the disease ontology data and split the DO_TERM for each column 
do_term <- gwas_snp_trait[,c("SNPID","DO_TERM")]
do_term <- concat.split(do_term, "DO_TERM", sep = "|", drop=TRUE)
do_term_long <- melt(do_term, id.vars="SNPID", value.name="DO_TERM")
do_term_long <- do_term_long[,c(1,3)]
do_term_long <- do_term_long[!duplicated(do_term_long),]
do_term_long <- do_term_long[!is.na(do_term_long$DO_TERM),]
do_all_source <- read.csv("/Volumes//KINGSTON/201409to201412/database/Human_parsed_output.tsv",
                                            sep="\t", header=FALSE)
do_all_source <- do_all_source[which(do_all_source$V1 != "nil"),]
colnames(do_all_source) <- c("DO_ID","DO_TERM","Description")
do_term_long$DO_TERM <- tolower(do_term_long$DO_TERM)
do_all_source$DO_TERM <- tolower(do_all_source$DO_TERM)
do_merge <- merge(do_term_long, do_all_source, by="DO_TERM", all.x=TRUE)
do_merge <-  do_merge[!is.na(do_merge$DO_ID),]
do_merge <- do_merge[!duplicated(do_merge),]
write.table(do_merge, "/Volumes//KINGSTON/201409to201412/database/do_merge_postgresql.tsv",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

##population information about publication
population <- read.csv("/Volumes/KINGSTON/201409to201412/database/population_table_postgresql.csv",
                       header=TRUE, sep="\t")
pop_sub <- population[,c("PUBMEDID", "SUB_POPULATION")]
pop_super <- population[,c("PUBMEDID","SUPER_POPULATION")]
pop_sub$SUB_POPULATION <- gsub("\\|", ",", pop_sub$SUB_POPULATION)
pop_super$SUPER_POPULATION <- gsub("\\|", ",", pop_super$SUPER_POPULATION)
pop_sub <- concat.split(pop_sub, "SUB_POPULATION", sep =",", drop=TRUE)
pop_super <- concat.split(pop_super, "SUPER_POPULATION", sep =",", drop=TRUE)

pop_sub <- melt(pop_sub, id.vars="PUBMEDID", value.name="SUB_POPULATION")
pop_sub <- pop_sub[which(pop_sub$SUB_POPULATION != ""),]
pop_sub$SUB_POPULATION <- gsub("\\)", "", pop_sub$SUB_POPULATION)
pop_sub <- concat.split(pop_sub, "SUB_POPULATION", sep ="(", drop=TRUE)
pop_sub <- pop_sub[, c(1,3,4)]

pop_super <- melt(pop_super, id.vars="PUBMEDID", value.name="SUPER_POPULATION")
pop_super <- pop_super[!is.na(pop_super$SUPER_POPULATION), ]
pop_super$SUPER_POPULATION <- gsub("\\)", "", pop_super$SUPER_POPULATION)
pop_super <- concat.split(pop_super, "SUPER_POPULATION", sep ="(", drop=TRUE)
pop_super <- pop_super[, c(1,3,4)]
pop_super <- pop_super[!duplicated(pop_super),]
population$PUBMEDID <- as.character(population$PUBMEDID)
pop_sub$PUBMEDID <- as.character(pop_sub$PUBMEDID)
pop_pubmed <- data.frame(population[, c("PUBMEDID")])
colnames(pop_pubmed) <- "PUBMEDID"
population_merge <- merge(pop_pubmed, pop_sub, by="PUBMEDID", all.x=TRUE)
population_merge <- population_merge[which(population_merge$SUB_POPULATION_1 != "ALL"),]
population_merge <- merge(population_merge, pop_super, by="PUBMEDID", all=TRUE)
population_merge <- population_merge[!is.na(population_merge$SUPER_POPULATION_2),]
population_merge <- population_merge[which(population_merge$SUPER_POPULATION_1!="ALL"),]
population_merge <- population_merge[!is.na(population_merge$SUB_POPULATION_2), ]
population_merge <- population_merge[!duplicated(population_merge),]
colnames(population_merge) <- c("PUBMEDID","Sub_Pop","Sub_Size","Super_Pop","Super_Size")
write.table(population_merge, "/Volumes/KINGSTON/201409to201412/database/population_merge_postgresql.tsv",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")