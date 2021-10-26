df1 <- read_excel("~/Downloads/vquest.xls")
df2 <- read_excel("~/Downloads/vquest.xls", sheet = 2)

head(df1)

df3 <- df1[names(df1) %in% c("Sequence number","Sequence ID","V-DOMAIN Functionality", "V-GENE and allele","V-REGION identity %","J-GENE and allele","J-REGION identity %","D-GENE and allele","JUNCTION frame","JUNCTION (with frameshift)","CDR3-IMGT (with frameshift)")]

df4 <- df2[names(df2) %in% c("Sequence number","Sequence ID","JUNCTION","JUNCTION (AA)","JUNCTION (with frameshift)","JUNCTION (AA) (with frameshift)","CDR3-IMGT","CDR3-IMGT (AA)")]

dim((merge(df3,df4,by=c("Sequence number","Sequence ID"))))


df_chain1 <- merge(df3,df4,by=c("Sequence number","Sequence ID"))


write.csv(df_chain1,"df_chain1.csv")
head(df_chain1)
df_chain1 <- as.data.frame(df_chain1)
df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
df_chain1$`J-GENE and allele` <- gsub('[(]','',df_chain1$`J-GENE and allele`)
df_chain1$`J-GENE and allele` <- gsub('[)]','',df_chain1$`J-GENE and allele`)
df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
df_chain1$`D-GENE and allele` <- gsub(' F','',df_chain1$`D-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('[(]','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('[)]','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub(' F,','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub(' F','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('[[]','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('[]]','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('F','',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub(' or ',', ',df_chain1$`V-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub(' ','',df_chain1$`V-GENE and allele`)
df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
df_chain1$`J-GENE and allele` <- gsub('TR','',df_chain1$`J-GENE and allele`)
df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)

df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION) 
df_chain1$V.sequence.quality.check <- ifelse(df_chain1$`V-REGION identity %`>=90,"quality sequence alignment","check chromatogram")
df_chain1$clone_quality <- NA 
df_chain1$comments <- NA
df_chain1
