require(tidyverse)

dataframe = read.table("test-data/sampleExport/TCRG_MRD_Day29_Case1.tsv",sep="\t",header=T)
x <- read.table("test-data/QC/ImmunoSEQ/ES8_TSNLQEQIGW_3.tsv",sep="\t",header=T)
unique(x$frame_type)

x2 <- x %>%
  select_if(~ !any(is.na(.)))

x2 <- x2 %>% mutate_all(na_if,"")

x2 <- subset(x2, x2$frame_type=="In")
x2 <- data.frame(cloneCount = x2[names(x2) %in% "templates"], x2)
x3 <- x2


x3$TRJ <- x3[names(x3) %in% "j_family"]
x3$TRJ

gsub("^TCR","",x3$j_gene)

x3$TRV <- gsub("^TCR","",x3$v_gene)
x3$TRJ <- gsub("^TCR","",x3$j_gene)
x3$TRD <- gsub("^TCR","",x3$d_gene)

x3 <- x3[!is.na(x3[names(x3) %in% c("TRV")]),]

x3$TRVJ <- paste(x3$TRV,x3$TRJ,sep=".")
x3$TRVDJ <- paste(x3$TRV,x3$TRD,x3$TRJ,sep=".")
x3$TRVDJ <- gsub(".NA.",".",x3$TRVDJ)
x3$TRVJ_CDR3 <- paste(x3$TRVJ,x3$amino_acid,sep="_")

x3 <- x3[!names(x3) %in% c("product_subtype","frame_type","total_dj_reads",
                           "productive_entropy","rearrangement_type",
                           "order_name","release_date",
                           "upload_date","primer_set",
                           "total_outofframe_reads","sample_catalog_tags","sample_rich_tags_json",
                           "fraction_productive","sample_tags","sku","total_templates",
                           "sequence_result_status","productive_clonality","stop_rearrangements",
                           "outofframe_rearrangements","total_rearrangements","total_reads","sample_cell",
                           "productive_rearrangements","counting_method","v_allele_ties","v_gene_ties","antibody",
                           "sample_clonality", "max_productive_frequency","sample_entropy","sample_simpson_clonality",
                           "max_frequency","productive_simpson_clonality","total_stop_reads","total_productive_reads",
                           "v_deletions",	"d5_deletions",	"d3_deletions",	"j_deletions",	"n2_insertions",
                           "n1_insertions",	"v_index",	"n1_index",	"n2_index",	"d_index",	"j_index",	"v_family_ties",	"d_family_ties",	"d_gene_ties",	"d_allele_ties",	"j_gene_ties"
)]

x3[is.na(x3$TRD),]
x3[is.na(x3)] <- "-"
unique(x3$j_family)
x3$v_gene <- gsub(" ","",x3$v_gene)
unique(x3$v_gene)
x3$TRJ <- gsub(" ","",x3$TRJ)

unique(x3$TRJ)
x3
