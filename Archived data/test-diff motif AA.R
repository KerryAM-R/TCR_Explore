df <- read.csv("test-data/Group/paired_unsummarised2021.09.22.csv",header=T) 


colorset(alphabet="AA",
         colorScheme="chemistry")

ASN$cols <- colorset(alphabet="AA",
                     colorScheme="chemistry")
ASN
names(df3) <- paste(names(df3),"length", sep="_")
df4 <- cbind(df,df3)
df4
df <- as.data.frame(df)
names(df)
df_unique <- as.data.frame(ddply(df,(c("group","JUNCTION_A")),numcolwise(sum)))
names(df_unique) <- c("group","chain","cloneCount")

df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% "chain"])

df_unique

x <- AAStringSet(df_unique$chain)

aln <- muscle(x)

df1 <- as.data.frame(aln@unmasked)

df_unique$chain1 <- df1$x

motif <- as.data.frame(t(as.data.frame(strsplit(df_unique[,grep("chain1",names(df_unique))], ""))))
z=dim(motif)[2]
z
df_unique1 <- subset(df_unique,df_unique$group=="IFN")
df_unique
motif1 <- as.data.frame(t(as.data.frame(strsplit(df_unique1[,grep("chain1",names(df_unique1))], ""))))
motif_count1 <- Nucleotide(cbind(x=1,y=2,motif1), z)
motif_count1<-pcm2pfm(motif_count1)
motif_count1

df_unique2 <- subset(df_unique,df_unique$group=="CD8")
motif2 <- as.data.frame(t(as.data.frame(strsplit(df_unique2[,grep("chain1",names(df_unique2))], ""))))
motif_count2 <- Nucleotide(cbind(x=1,y=2,motif2), z)

motif_count2<-pcm2pfm(motif_count2)
?diffLogoTableConfiguration

diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1), pwm2 = as.data.frame(motif_count2), alphabet = DNA)
diffLogoObj
diffLogo(diffLogoObj)

?diffLogo()

df <- read.csv("~/Desktop/4905 Vd1.csv")
df <- as.data.frame(df)
names(df)
df_unique <- as.data.frame(ddply(df,(c("group","JUNCTION_G")),numcolwise(sum)))
names(df_unique) <- c("group","chain","cloneCount")

df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% "chain"])
x <- AAStringSet(df_unique$chain)
aln <- muscle(x,quiet = T)
df1 <- as.data.frame(aln@unmasked)
as.data.frame(aln@unmasked)

df_unique$chain1 <- df1$x
df_unique
