require(tidyverse)

dataframe = read.table("test-data/sampleExport/TCRG_MRD_Day29_Case1.tsv",sep="\t",header=T)
x <- read.table("test-data/sampleExport/TCRG_MRD_Day29_Case1.tsv",sep="\t",header=T)
head(x)

x2 <- x %>%
  select_if(~ !any(is.na(.)))

unique(x2$frame_type)

x2 <- subset(x2, x2$frame_type=="In" & x2$j_gene!="unresolved")

x2 <- data.frame(cloneCount = x2$seq_reads, x2)

x3 <- x2[!names(x2) %in% c("product_subtype","frame_type","total_dj_reads",
                           "productive_entropy","rearrangement_type",
                           "order_name","release_date",
                           "upload_date","primer_set",
                           "total_outofframe_reads",
                           "fraction_productive","seq_reads",
                           "sequence_result_status","productive_clonality","stop_rearrangements",
                           "outofframe_rearrangements","total_rearrangements","total_reads",
                           "productive_rearrangements","counting_method","v_allele_ties","v_gene_ties",
                           "sample_clonality", "max_productive_frequency","sample_entropy","sample_simpson_clonality",
                           "max_frequency","productive_simpson_clonality","total_stop_reads","total_productive_reads"
                           )]

x3$TRJ <- gsub("^TCR","",x3$j_family)
x3$TRV <- gsub("^TCR","",x3$v_family)
x3$TRD <- gsub("^TCR","",x3$d_gene)
x3$TRVJ <- paste(x3$TRV,x3$TRJ,sep=".")
x3$TRVDJ <- paste(x3$TRV,x3$TRD,x3$TRJ,sep=".")
x3$TRVDJ <- gsub(".unresolved.",".",x3$TRVDJ)
x3$TRVJ_CDR3 <- paste(x3$TRVJ,x3$amino_acid,sep="_")

x3
