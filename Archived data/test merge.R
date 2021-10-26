df1 <- read.csv("~/Desktop/4905_Vd1_CBZ.csv")
df1 <- as.data.frame(df1)
df <- subset(df1,df1$clone_quality=="pass")
df <- as.data.frame(df)
df2 <- df[!names(df) %in% c("V.sequence.quality.check","clone_quality","comments","JUNCTION..with.frameshift.","CDR3.IMGT..with.frameshift.","JUNCTION..AA...with.frameshift.")]

df.Vgene <- as.data.frame(do.call(rbind, strsplit(as.character(df2$V.GENE.and.allele), ",")))
df2$V.GENE <- df.Vgene$V1
y = dim(df2)[2]
y
df2$V.GENE <- gsub(" ","",df2$V.GENE)
df2$cloneCount <- 1


df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
head(df_name2)
df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "-")))
df_name4 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), "-")))
df_name3$V1 <- gsub("G","",df_name3$V1)
df_name3$V1 <- gsub("D","",df_name3$V1)

df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
head(df_name5)
df2$ID <- df_name2$V1

df2$Indiv.group <- df_name3$V1
df2$Indiv <-df_name5$V1
df2$group <- df_name5$V2
df2$clone <- df_name4$V2

chain1 <- df2[grep("G",df2$V.GENE.and.allele),]
chain2 <- df2[grep("D",df2$V.GENE.and.allele),]
chain1$ID <- gsub("G-","-",chain1$ID)
names(chain1)[1:y] <- paste(names(chain1)[1:y],"G",sep="_")
head(chain1)
chain2$ID <- gsub("D-","-",chain2$ID)
names(chain2)[1:y] <- paste(names(chain2)[1:y],"D",sep="_")
head(chain2)

z = y+1
x <- names(chain2)[z:dim(chain2)[2]]
merged_chain <- merge(chain1,chain2,by =x)
head(merged_chain)
merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_G","Sequence.ID_D","V.DOMAIN.Functionality_G","V.DOMAIN.Functionality_D","D.GENE.and.allele_G","JUNCTION.frame_G","JUNCTION.frame_D"))]
dat <- merged_chain2
dat$GJ <- paste(dat$V.GENE_G,".",dat$J.GENE.and.allele_G,sep="")
dat$DJ <- paste(dat$V.GENE_D,".",dat$J.GENE.and.allele_D,sep="")
dat$GJ <- gsub("[*]0.","",dat$GJ)
dat$DJ <- gsub("[*]0.","",dat$DJ)

dat$GJ <- gsub(", or GJ","/",dat$GJ)
head(dat)
dat$GJDJ <- paste(dat$GJ,".",dat$DJ,sep="")
dat$GJ_gCDR3 <- paste(dat$GJ,dat$JUNCTION..AA._G,sep="_")
dat$DJ_dCDR3 <- paste(dat$DJ,dat$JUNCTION..AA._D,sep="_")

dat$GJ_gCDR3_DJ_dCDR3 <- paste(dat$DJ_dCDR3,dat$GJ_gCDR3,sep=" & ")
dat
