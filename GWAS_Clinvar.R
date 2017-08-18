###read Clinical information table filedata=20141105 dbSNP_BUILD_ID=142 reference=GRCh37.p13###
wd <- read.csv("GWAS_disease_information//clinvar_20141105.vcf",stringsAsFactors=FALSE, sep="\t", skip=70)

###extract information between CLNDBN and CLNREVSTAT from table wd and column INFO###
Clinvar <- ""
for (i in 1:nrow(wd)) {
   Clinvar[i] = substring(wd[,8][i], regexpr('CLNDBN',wd[,8][i])[[1]]+7, regexpr('CLNREVSTAT',wd[,8][i])[[1]]-2)
}
Clinvar.data.frame <- as.data.frame(Clinvar)
wd <- cbind(wd[,1:7],Clinvar.data.frame)


###read GWAS data downloaded from www.genome.gov/gwastudies/ filedata=20141113###
gwas <- read.csv("Downloads/gwascatalog.txt", sep="\t", header=TRUE)
gwas <- gwas[,c("Disease.Trait","Mapped_gene", "SNPs","Snp_id_current", "Initial.Sample.Description", "Chr_id",
                "Context", "Strongest.SNP.Risk.Allele", "p.Value", "OR.or.beta",
                "X95..CI..text.", "PUBMEDID", "Date")]
gwas$Initial.Sample.Description <- as.character(gwas$Initial.Sample.Description)

population <- as.data.frame(matrix(NA, ncol=12, nrow=18725))
colnames(population) <- c("population", "european", "east_asian", "chinese", "hispanic", "africa",
                          "south_asian", "latin", "american_indian", "jewish", "middle_eastern", "pacific")

for (i in 1:nrow(gwas)) {
  ifelse((grep('Korculan|Romanian|Finland|Finnish|European|Friuli|Hutterite|Val Borbera|Turkish|Amish|Sardinian|French|Bulgarian', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1),
  population$european[i] <- "European", population$european[i])  
  ifelse((grepl('Japanese|Korean|East Asian', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1), 
    population$east_asian[i] <- "East Asian", population$east_asian[i])
  ifelse((grepl('Chinese', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1), 
          population$chinese[i] <- "Chinese", population$east_asian[i]) 
  ifelse((grepl('Hispanic', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1), 
         population$hispanic[i] <- "Hispanic", population$hispanic[i])
  ifelse((grepl('Indian Asian|Indonesian|Filipino|Punjabi Sikh|Bangladeshi|South Asian|Sri Lankan|Thai|Asian Indian|Malay|Singaporean', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1),
    population$south_asian[i] <- "South Asian", population$south_asian[i])
  ifelse((grepl('Black|African', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1),
    population$africa[i] <- "African", population$africa[i])
  ifelse((grepl('Mexican|Latin|Brazilian|Afro-Caribbean', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1), 
         population$latin[i] <- "Latin", population$latin[i]) 
  ifelse((grepl('American Indian', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1),
    population$american_indian[i] <- "American Indian", population$american_indian[i]) 
  ifelse((grepl('Jewish', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1), 
          population$jewish[i] <- "Jewish", population$jewish[i])  
  ifelse((grepl('Mongolian|Middle Eastern', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1), 
         population$middle_eastern[i] <- "Middle Eastern", population$middle_eastern[i])
  ifelse((grepl('Kosraen|Micronesian', gwas$Initial.Sample.Description[i], ignore.case=TRUE)==1),
         population$pacific[i] <- "Pacific", population$pacific[i])
  population$population[i] <- gsub(",NA|NA,|NA", "", paste(population$european[i], population$east_asian[i], population$chinese[i],
                                    population$hispanic[i], population$south_asian[i], population$africa[i],
                                    population$latin[i], population$american_indian[i], population$jewish[i],
                                    population$middle_eastern[i], sep=","))
} 

gwas <- cbind(gwas[,c(1:4,6:13)], population$population)
colnames(gwas)[13] <- "Ethnicity"
gwas$Mapped_gene <- with(gwas, gsub(" - ", ",", gwas$Mapped_gene))
gwas$Snp_id_current <- with(gwas, paste("rs",gwas$Snp_id_current, sep=""))
gwas$SNPs <- as.character(gwas$SNPs)
gwas$match <- ""
for (i in 1:nrow(gwas)) {ifelse((gwas$SNPs[i]==gwas$Snp_id_current[i]), gwas$match[i] <- 1, gwas$match[i] <- 0)}
gwas_SNP_id_same <- gwas[which(gwas$match==1),]
gwas_SNP_id_different <- gwas[which(gwas$match==0),]

write.table(gwas_SNP_id_same[,3], "GWAS_disease_information/gwas_SNP_id_same_table.tsv", quote=FALSE, row.names=FALSE,
            col.names=TRUE, sep="\t")

##library(biomaRt)
##rs_number <- as.data.frame(gwas_SNP_id_same$SNPs)
##colnames(rs_number) <- "rs_number"
##mart.snp <- useMart("snp","hsapiens_snp")
##get_snp_information <- function(x) {getBM(attributes=c("refsnp_id", "refsnp_source", "chr_name","chrom_start", 
##                                  "allele", "minor_allele_freq"), filters="snp_filter", values=x, mart=mart.snp)}

##read data file downloaded from http://jjwanglab.org/gwasdb and handled by awk##
gwasdb_201409_snp_trait <- read.csv("GWAS_disease_information/gwasdb_20140908_snp_trait_filtered.tsv", stringsAsFactors=FALSE,
                          sep="\t")
gwasdb_201409_snp_annotation <- read.csv("GWAS_disease_information/gwasdb_20140908_annotation_filtered.tsv", stringsAsFactors=FALSE,
                                    sep="\t")
gwasdb_201409_snp <- merge(gwasdb_201409_snp_trait, gwasdb_201409_snp_annotation, by.x="SNPID", by.y="SNP_Id", all.x=TRUE)
gwasdb_201409_snp <- gwasdb_201409_snp[!duplicated(gwasdb_201409_snp[c(1,7,8)]),]

##delete (number) and ALL(number) from subpopulation or superpopulation and replace "|" with ","##
clean_population_info <- function(x) {
  x <- gsub("\\|", ",", x, perl=TRUE)
  x <- gsub("(,)*ALL\\([0-9]*\\)", "", x)
  return(x)
}
gwasdb_201409_snp_pop <- as.data.frame(lapply(gwasdb_201409_snp[,c("SUB_POPULATION","SUPER_POPULATION")], clean_population_info))
gwasdb_201409_snp <- cbind(gwasdb_201409_snp[,c(-11,-12)], gwasdb_201409_snp_pop)
gwasdb_201409_snp$DO_TERM <- gsub("`s","", gwasdb_201409_snp$DO_TERM)

##HPO_TERM and DO_TERM columns both are not NA and if the first word of DO_TERM matches with any## 
##word in HPO_TERM then assign NA to column HPO_TERM##
for (i in 1:nrow(gwasdb_201409_snp)) {
  if (!is.na(gwasdb_201409_snp$DO_TERM[i]) && !is.na(gwasdb_201409_snp$HPO_TERM[i])) {
    if (length(strsplit(gwasdb_201409_snp$DO_TERM[i], " ")[[1]])==1 && grepl(gwasdb_201409_snp$DO_TERM[i], gwasdb_201409_snp$HPO_TERM[i],ignore.case=TRUE)) {
      gwasdb_201409_snp$HPO_TERM[i] <- NA
    } else if (length(strsplit(gwasdb_201409_snp$DO_TERM[i], " ")[[1]])>1 && grepl(strsplit(gwasdb_201409_snp$DO_TERM[i]," ")[[1]][1], gwasdb_201409_snp$HPO_TERM[i],ignore.case=TRUE)) {
      gwasdb_201409_snp$HPO_TERM[i] <- NA
    }
  }
}

##paste2 function combined two values from two columns and ignore NA##
paste2 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  gsub(paste0("(^",sep,"|",sep,"$)"),"",
       gsub(paste0(sep,sep),sep,
            do.call(paste,c(L,list(sep=sep)))))
}

gwasdb_201409_snp$HPO_DO_TERM <- paste2(gwasdb_201409_snp$HPO_TERM, gwasdb_201409_snp$DO_TERM)
gwasdb_201409_snp$ORI_SNPID <- ifelse(gwasdb_201409_snp$SNPID==gwasdb_201409_snp$ORI_SNPID, NA, gwasdb_201409_snp$ORI_SNPID)
gwasdb_201409_snp$GWAS_TRAIT <- gsub("`s","", gwasdb_201409_snp$GWAS_TRAIT)
gwasdb_201409_snp$HPO_TERM <- gsub("`s","", gwasdb_201409_snp$HPO_TERM)

gwas_SNP_id_same_gwasdb201409 <- merge(gwas_SNP_id_same, gwasdb_201409_snp, by.x="SNPs", by.y="SNPID", all.x=TRUE)
gwas_SNP_id_same_gwasdb201409 <- gwas_SNP_id_same_gwasdb201409[!duplicated(gwas_SNP_id_same_gwasdb201409[c(1,2)]),]
gwas_SNP_id_same_gwasdb201409 <- gwas_SNP_id_same_gwasdb201409[,c("Disease.Trait","Mapped_gene","Snp_id_current","CHR","POS","REF","ALT","Context",
                                                                  "Ethnicity", "Strongest.SNP.Risk.Allele","P_VALUE","OR.BETA", "CI95_TEXT",
                                                                  "PMID", "ORI_SNPID")]
gwas_SNP_id_same_gwasdb201409 <- gwas_SNP_id_same_gwasdb201409[!duplicated(gwas_SNP_id_same_gwasdb201409[c(1,3)]),]

library(splitstackshape)
gwas_SNP_id_different_gwasdb201409 <- as.data.frame(cSplit(as.data.table(gwas_SNP_id_different), "SNPs", ","))
gwas_SNP_id_different_gwasdb201409_trait_SNP <- gwas_SNP_id_different_gwasdb201409[,c(1,14:33)]
gwas_SNP_id_different_gwasdb201409_no_SNP <- gwas_SNP_id_different_gwasdb201409[,c(1:12)]
library(reshape2)
gwas_SNP_id_different_gwasdb201409_trait_SNP <- melt(gwas_SNP_id_different_gwasdb201409_trait_SNP,
                                                     id.vars="Disease.Trait")
gwas_SNP_id_different_gwasdb201409_trait_SNP <- gwas_SNP_id_different_gwasdb201409_trait_SNP[grepl('rs', gwas_SNP_id_different_gwasdb201409_trait_SNP$value, ignore.case=TRUE)==1,]
gwas_SNP_id_different_gwasdb201409_SNP <- merge(gwas_SNP_id_different_gwasdb201409_no_SNP, gwas_SNP_id_different_gwasdb201409_trait_SNP,
                                            by="Disease.Trait", all.y=TRUE)
gwas_SNP_id_different_gwasdb201409_SNP <- gwas_SNP_id_different_gwasdb201409_SNP[!duplicated(gwas_SNP_id_different_gwasdb201409_SNP[c(1,14)]),]

gwas_SNP_id_different_gwasdb201409 <- merge(gwas_SNP_id_different_gwasdb201409_SNP, gwasdb_201409_snp, by.x="value", by.y="ORI_SNPID", all.x=TRUE)
gwas_SNP_id_different_gwasdb201409 <- gwas_SNP_id_different_gwasdb201409[!duplicated(gwas_SNP_id_different_gwasdb201409[c(1,2)]),]


colnames(gwas_SNP_id_different_gwasdb201409)[1] <- "ORI_SNPID"
gwas_SNP_id_different_gwasdb201409 <- gwas_SNP_id_different_gwasdb201409[,c("Disease.Trait","GENE_SYMBOL", "SNPID","CHR", "POS", "REF","ALT","Context",
                                                                            "Ethnicity", "Strongest.SNP.Risk.Allele","P_VALUE","OR.BETA", "CI95_TEXT",
                                                                            "PMID", "ORI_SNPID")]
gwas_SNP_id_different_gwasdb201409 <- gwas_SNP_id_different_gwasdb201409[!is.na(gwas_SNP_id_different_gwasdb201409$CHR),]
colnames(gwas_SNP_id_different_gwasdb201409)[1:15] <- colnames(gwas_SNP_id_same_gwasdb201409)[1:15]
##combined the SNP id is the same for Snp_id_current and Ori_SNPID and SNP is different
##for Snp_id_current and Ori_SNPID.
gwas_SNP_id_different_gwasdb201409 <- gwas_SNP_id_different_gwasdb201409[which(!gwas_SNP_id_different_gwasdb201409$Snp_id_current %in% gwas_SNP_id_same_gwasdb201409$Snp_id_current),]
gwas_SNP_id_gwasdb201409_all <- rbind(gwas_SNP_id_same_gwasdb201409,gwas_SNP_id_different_gwasdb201409)
gwas_SNP_id_gwasdb201409_all <- merge(gwas_SNP_id_gwasdb201409_all, wd[, c("ID","Clinvar")], by.x="Snp_id_current", by.y="ID", all.x=TRUE)
gwas_SNP_id_gwasdb201409_all$ORI_SNPID <- ifelse(gwas_SNP_id_gwasdb201409_all$Snp_id_current!=gwas_SNP_id_gwasdb201409_all$ORI_SNPID, 
                                                 gwas_SNP_id_gwasdb201409_all$ORI_SNPID, NA)
colnames(gwas_SNP_id_gwasdb201409_all) <- c("rs","disease_trait","gene","chrom","pos","reference_allele",
                                            "alternative_allele","region_to_gene","ethnicity","risk_allele","p_value","odds_ratio",
                                            "odds_ratio_95CI","pubmed_id","ori_rs","clinvar")
##combine gwas_SNP_id_gwasdb201409_all and gwasdb_201409_snp##
gwas_SNP_all <- merge(gwasdb_201409_snp, gwas_SNP_id_gwasdb201409_all[,c("rs","gene","ethnicity","risk_allele","clinvar")],
                      by.x="SNPID", by.y="rs", all.x=TRUE)
for (i in 1:nrow(gwas_SNP_all)) {
  if (is.na(gwas_SNP_all$GENE_SYMBOL[i])) {
    gwas_SNP_all$GENE_SYMBOL[i]<-gwas_SNP_all$gene[i]
  }
  if (gwas_SNP_all$SUPER_POPULATION[i]=="") {
    gwas_SNP_all$SUB_POPULATION[i]<-gwas_SNP_all$ethnicity[i]
  }
}

gwas_SNP_all <- gwas_SNP_all[,c(-18,-19)]
gwas_SNP_all_chinese <- gwas_SNP_all[grepl('Chinese', gwas_SNP_all$SUB_POPULATION,ignore.case=TRUE), ]
gwas_SNP_all_east_asian <- gwas_SNP_all[grepl('Chinese|Japanese|Korean|East Asian', gwas_SNP_all$SUB_POPULATION,ignore.case=TRUE), ]
##add subtype of disease_trait##
gwas_SNP_all$disease_type <- ""
for (i in 1:nrow(gwas_SNP_all)) {
  if (grepl('cancer|carcinoma|metalloproteinase|Osteosarcoma|sarcoma|cell-free DNA|Glioblastoma|Essential tremor|Melanoma|microglubulin|Uterine fibroids|carcinoids|Glioma|Multiple myeloma|Meningioma|myeloma|endothelial|Mammographic|lymphoma|Telomere length|tumor|adenocarcinoma|PCA3|leukemia|phramacokinetics', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) { 
    gwas_SNP_all$disease_type[i] <- "cancer"
  } else if (grepl('color|Freckles|eyes|freckling|Acne|Facial|Psoriatic|pigmentation|Breast size|Keloid|baldness|Vitiligo|Hair morphology|Tanning', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "appearance"
  } else if (grepl('Heart|Cardiac|Coronary|stiffness|Anticoagulant|Atrioventricular|Major CVD|Arterial|Urinary albumin|hypotension|vasoactive peptide|Homoarginine|arrhythmia|Tetralogy|Echocardiographic|QT interval|Retinal vascular|QRS|Myocardial infarction|Total ventricular|Atrial fibrillation|Electrocardiographic|pressure|Morbidity-free|Stroke|Plasma homocysteine|Pulmonary|Aortic|cardiomyopathy|Cardiovascular|Brugada syndrome|PR interval', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "cardiovascular"
  } else if (grepl('Alopecia areata|Alcoholism|Nicotine|Nicotine|Self-rated|Gambling|Methamphetamine|Behcet|Coffee|Job-related|Pain|Addiction|Nicotine|Substance dependence|Caffeine|Hypersomnia|Sunburns|sleep|Cocaine|Exercise|Temperament|Insomnia|alcohol|Drinking|Molar-incisor|Smoking|Dental|Pit-and-Fissure|Longevity|Smooth-surface|Skin sensitivity|Longevity|Aging', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "health life"
  } else if (grepl('Ovarian reserve|Sex hormone-binding|Estradiol|Sexual|Premature ovarian|Testosterone|Erectile|Recombination|Soluble leptin|Polycystic ovary|Androgen|Male infertility|azoospermia', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "human fertilization"
  } else if (grepl('Cholesterol|Lipid|HDL|LDL|Metabolite|Lp (a) levels|Paraoxonase|phosphatase|plasminogen|Apolipoprotein|Angiotensin|arteriolar caliber|Hepcidin|myeloperoxidase|Serum albumin|aminotransferase|N-glycan|Glycemic|Hypertriglyceridemia|ceruloplasmin|Urinary uromodulin|myeloperoxidase|gluatamyl transferase|Cystatin C|Butyrylcholinesterase|uric acid|Palmitoleic|Oleic acid|transpeptidase|Glycemic|Biochemical|prostate-specific antigen|C3 and C4|Proinsulin|insulin|enzyme|Thyroid|Metabolic|C-reactive|diabetes|Triglycerides|Digestive|glucose|Serum dimethylarginine|Natriuretic peptide|Adiponectin|Barrett|Stearic acid|Phospholipid|Creatinine|Protein|Bilirubin|beta-cell|Diabetic|Renal|Urate', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "metabolism and diet"
  } else if (grepl('Inflammatory|Immune|antibodies|antibody|Chemerin|erythematosus|Tuberculosis|microbiota|Helicobacter|seropositivity|Leishmaniasis|Dermatomyositis|hypersecretion|Stevens-Johnson|E-selectin|YKL-40|Atopy|Meningococcal|lymphocyte|adhesion|Monocyte|Sarcoidosis|Interleukin-18|ICAM-1|Sclerosing cholangitis|Inflammatory|Gout|Neutrophil|Hepatitis B|Periodontal disease|Immunoglobulin|Psoriasis|Dialysis|Graves|White blood|Endometriosis|Kawasaki|IgA|IgG|IgE|AIDS|HIV|IgM|f6gren|Ulcerative|allergy|Allergic|Tripanosoma|Celiac disease', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "autoimmune"
  } else if (grepl('chemotherapy|induced|hydrochloride|Radiation|therapy|Opioid|Bleomycin|Treatment|Lumiracoxib-related|Capecitabine|Schizophrenia|radiation|Functional MRI|Vaccine|fenofibrate|Drug|antidepressant|Bronchodilator response|Aspirin|platinum-based agents|Response to', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "drug or therapy response"
  } else if (grepl('lung function|lung disease|kidney|Prion|complex|hormone|lumbar|methylation|Glomerular|Gallstones|Nephrolithiasis|Duodenal ulcer|Nephropathy|Dehydroepiandrosterone sulphate|Emphysema-related|Preeclampsia|Nephrotic|Malaria|Crohn|Migraine with aura|forced vital capacity|Airflow|biliary cirrhosis|Rheumatoid|spondylitis|Asthma', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "common disease"
  } else if (grepl('Red blood|Platelet|Fibrinogen|Hematocrit|Hypertension|F-cell|cytopenia|HbA2|Moyamoya|hemorrhage|Prothrombin|anemia|Glomerulosclerosis|erythrocyte|Haptoglobin|granulomatosis|Peripheral artery|Thrombin|thromboplastin|Hematology|atherosclerosis|Homocysteine|Plasminogen|Carotid intima|D-dimer|Myeloproliferative|corpuscular|Coagulation|sedimentation|aneurysm|Hematological|Factor VII|FVIII|Venous thromboembolism|Hemoglobin|aneurysms', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "blood system"
  } else if (grepl('Parkinson|Amyotrophic|Gray matter|dyskinesia|Hirschsprung|Autism|Sphingolipid|Myasthenia|Leprosy|Neuroticism|Rhegmatogenous retinal|Restless legs|Creutzfeldt-Jakob|Narcolepsy|Tourette|supranuclear palsy|Epilepsy|Brain lesion|dystonia|White matter hyperintensity|Neuroblastoma|Narcolepsy|Conduct disorder|cortical|Migraine|Sensory disturbances|Lentiform nucleus|Migraine|sclerosis|Alzheimer|Vascular dementia|Cerebrospinal|MRI atrophy|amyloid', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "neurological disorder"
  } else if (grepl('Waist-to-hip|ratio|Body mass index|Waist|mass|BMI|Height|circumference|Anthropometric|Quantitative|Ankle-brachial|Hip|Axial length|Waist-hip|Fasting glucose-related', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "body mass index"
  } else if (grepl('fat|fatty|Obesity|Weight|Adiposity|Gaucher|Palmitic|Resistin', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "obesity"
  } else if (grepl('Dupuytren|Adolescent idiopathic|Bone properties|Scoliosis|Temporomandibular|Clubfoot|Osteonecrosis|Intracranial|Otosclerosis|Osteoporosis|brain volume|cytoarchitecture|geometry|Paget|structure|Periodontitis|Bone mineral|Pancreatitis|Left ventricular|Osteoarthritis|Cleft|Subcutaneous adipose', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "body development"
  } else if (grepl('Corneal|hand skill|Hippocampal|Vertical cup-disc|Behavioural|Myopia|Retinopathy|Cutaneous|hypokalemic|Handedness|repetition|Caudate nucleus|Partial epilepsies|Anterior chamber|oscillations|Odorant|reading|Contrast|Intraocular|Refractive|Hearing|macular degeneration|Tonometry|Optic|Glaucoma|Bitter taste|Cognitive', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "cognitive"
  } else if (grepl('craniosynostosis|Dengue shock|tooth|Bronchopulmonary|Hypospadias|Biliary|Neonatal|newborns|language|Tooth agenesis|juvenile idiopathic|Cystic fibrosis|Puberty|Infantile|childhood|young onset|Primary tooth|dermatitis|Asperger|Otitis media|Orofacial|Menopause|Menarche|infant|Birth weight|Common traits', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "child grwoth and development"
  } else if (grepl('Intelligence|Cannabis|Educational|Panic|Personality|ideology|internalizing|Memory|hyperactivity|brain serotonin|hyperresponsiveness|mood|behavioral|Obsessive-compulsive|Post-traumatic|Hoarding|Anger|Psychosis|Information processing|dyslexia|Self-employment|communication|Callous-unemotional|Suicide|depressive|Social autistic-like|Bipolar|preferences|Inattentive|Iris characteristics|Attention deficit|thought|Hyperactive-impulsive', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "mental"
  } else if (grepl('Eating disorders|Anorexia nervosa|Eosinophil|Iron|Phytosterol|Cu levels|Carotenoid|Phosphorus|Dietary|selenium|Bulimia nervosa|Acenocoumarol|Vitamin|Retinol|Magnesium|Vaspin|Zn levels|Calcium levels|Se levels', gwas_SNP_all$GWAS_TRAIT[i], ignore.case=TRUE)==1) {
    gwas_SNP_all$disease_type[i] <- "nutrition"
  } else gwas_SNP_all$disease_type[i] <- "unknown"
}

write.table(gwas_SNP_all, "GWAS_disease_information/gwas_snp_all.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

gwas_snp_maf <- read.table("GWAS_disease_information/gwas_snp_maf.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
gwas_snp_maf <- gwas_snp_maf[,c(-2,-3,-4,-5)]
gwas_SNP_all_maf <- merge(gwas_SNP_all, gwas_snp_maf, by.x="SNPID", by.y="Rs", all.x=TRUE)
gwas_SNP_all_maf <- gwas_SNP_all_maf[!duplicated(gwas_SNP_all_maf[c(1,2,3)]),]
gwas_SNP_all_maf <- gwas_SNP_all_maf[,c(-6)]
colnames(gwas_SNP_all_maf)[27] <- "ORI_SNPID"
gwas_SNP_all_ori_snp <- gwas_SNP_all[!is.na(gwas_SNP_all$ORI_SNPID),]
gwas_SNP_all_maf2 <- merge(gwas_SNP_all_ori_snp, gwas_snp_maf, by.x="ORI_SNPID", by.y="Ori_rs", all.x=TRUE)
gwas_SNP_all_maf2 <- gwas_SNP_all_maf2[!duplicated(gwas_SNP_all_maf2[c(1,2,3,4)]),]
gwas_SNP_all_maf2 <- gwas_SNP_all_maf2[,c(2:19, 21:28, 1)]
gwas_SNP_all_maf <- rbind(gwas_SNP_all_maf, gwas_SNP_all_maf2)

gwas_SNP_all_chinese <- gwas_SNP_all_maf[grepl('Chinese', gwas_SNP_all_maf$SUB_POPULATION,ignore.case=TRUE), ]
gwas_SNP_all_east_asian <- gwas_SNP_all_maf[grepl('Chinese|Korean|Japanese|East Asian', gwas_SNP_all_maf$SUB_POPULATION,ignore.case=TRUE), ]
write.table(gwas_SNP_all_maf, "GWAS_disease_information/gwas_snp_maf_all.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
write.table(gwas_SNP_all_chinese, "GWAS_disease_information/gwas_snp_maf_chinese.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
write.table(gwas_SNP_all_east_asian, "GWAS_disease_information/gwas_snp_maf_east_asian.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

wellness_markers <- read.table("GWAS_disease_information/replacement.txt")
replacement <- gwas_SNP_all_maf[which(gwas_SNP_all_maf$SNPID %in% wellness_markers$V1),]
replacement <- replacement[,c(12,13,1,4,5,16,6,2,3,19,20,21,22,23,24,25)]
write.table(replacement, "GWAS_disease_information/replacement.txt", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

####make data for postgresql relational database###
gwas_snp_full <- merge(gwas, gwasdb_201409_snp, by.x="SNPs", by.y="SNPID", all=TRUE)
gwas_snp_full$p.Value <- factor_to_character(gwas_snp_full$p.Value)
gwas_snp_full$X95..CI..text.<- factor_to_character(gwas_snp_full$X95..CI..text.)
gwas_snp_full$Ethnicity <- factor_to_character(gwas_snp_full$Ethnicity)
gwas_snp_full$OR.or.beta <- factor_to_character(gwas_snp_full$OR.or.beta)
gwas_snp_full$Disease.Trait <- factor_to_character(gwas_snp_full$Disease.Trait)

gwas_snp_full$P_VALUE <- ifelse(is.na(gwas_snp_full$P_VALUE), gwas_snp_full$p.Value, gwas_snp_full$P_VALUE)
gwas_snp_full$OR.BETA <- ifelse(is.na(gwas_snp_full$OR.BETA), gwas_snp_full$OR.or.beta, gwas_snp_full$OR.BETA)
gwas_snp_full$CI95_TEXT <- ifelse(is.na(gwas_snp_full$CI95_TEXT), gwas_snp_full$X95..CI..text., gwas_snp_full$CI95_TEXT)
gwas_snp_full$SUB_POPULATION <- ifelse(is.na(gwas_snp_full$SUB_POPULATION), gwas_snp_full$Ethnicity, gwas_snp_full$SUB_POPULATION)
gwas_snp_full$GWAS_TRAIT <- ifelse(is.na(gwas_snp_full$GWAS_TRAIT), gwas_snp_full$Disease.Trait, gwas_snp_full$GWAS_TRAIT)
gwas_snp_full$PUBMEDID <- ifelse(is.na(gwas_snp_full$PMID), gwas_snp_full$PUBMEDID, gwas_snp_full$PMID)
gwas_snp_full$GENE_SYMBOL <- ifelse(is.na(gwas_snp_full$GENE_SYMBOL), gwas_snp_full$Mapped_gene, gwas_snp_full$GENE_SYMBOL)
gwas_snp_maf_full$Context <- factor_to_character(gwas_snp_maf_full$Context)
gwas_snp_maf_full$TG_EUR <- ifelse(is.na(gwas_snp_maf_full$TG_EUR), gwas_snp_maf_full$Context, gwas_snp_maf_full$TG_EUR)

contribute_by <- gwas_snp_full[,c(1,11,21,22,23)]
contribute_by <- contribute_by[!duplicated(contribute_by),]
write.table(contribute_by,"Desktop/contribute_table_postgresql.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
publication <- gwas_snp_full[, c("PUBMEDID","GWAS_TRAIT","Date")]
publication <- publication[!duplicated(publication),]
write.table(publication, "Desktop/publication_table_postgresql.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
population <- gwas_snp_full[, c(11,24,25)]
population <- population[!duplicated(population),]
write.table(population, "Desktop/population_table_postgresql.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
gwas_snp_maf_full <- merge(gwas_snp_full, gwas_snp_maf, by.x="SNPs", by.y="Rs", all=TRUE)
gwas_snp_maf_full <- gwas_snp_maf_full[!duplicated(gwas_snp_maf_full),]
markers <- gwas_snp_maf_full[,c("SNPs","CHR","POS","REF","ALT","GENE_SYMBOL","TG_EUR","Strongest.SNP.Risk.Allele","CHBSX",
                               "JPT","CEU","EUR","EAS","AMR","SAS","Ori_rs","HPO_TERM","DO_TERM")]
markers <- markers[!duplicated(markers),]
write.table(markers, "Desktop/markers_table_postgresql.csv", quote=FALSE,  
            row.names=FALSE, col.names=TRUE, sep="\t")