dataframe = read.csv("test-data/Group/paired_unsummarised2021.09.22.csv",header=T) 
head(dataframe)

df.names <-  dataframe[ , -which(names(dataframe) %in% c("cloneCount","clone"))]
df1 <- ddply(dataframe,names(df.names) ,numcolwise(sum))
df1 <- df1[order(df1$cloneCount, decreasing = T),]

names(df1)

df.group <- unique(df1[names(df1) %in% "Indiv.group"])
names(df.group) <- "V1"

column.length <- length(df.group$V1)
column.length

row.length <- length(unique(df1$AJBJ))
row.length

m = matrix(NA,ncol=column.length, nrow=row.length)
samps <- df.group$V1

for (j in 1:column.length){
  df2 <- subset(df1,df1$Indiv.group==samps[j])
  m[,j] <- c(df2$cloneCount, rep(NA, row.length - length(df2$cloneCount)))
}
m <- as.data.frame(m)
names(m) <- samps
head(m)
m

a <- matrix(nrow=1,ncol=dim(m)[2])
b <- matrix(nrow=1,ncol=dim(m)[2])
d <- matrix(nrow=1,ncol=dim(m)[2])

for( i in 1:dim(m)[2]) {
  
  samp <- m[,i]
  samp <- na.omit(samp)
  a[,i] <- diversity(samp,"invsimpson")
  b[,i] <- sum(samp)
  d[,i] <- nrow(as.data.frame(samp))
}

a1 <- rbind(a,b,d)  
a1 <- as.data.frame(a1)
names(a1) <- names(m)
a1 <- rbind(a,b,d)  
a1 <- as.data.frame(a1)
names(a1) <- names(m)
a1
df_name <- as.data.frame(do.call(rbind, strsplit(as.character(names(m)), "\\.")))
head(df_name) 

a2 <- as.data.frame(t(a1))
names(a2) <- c("inv.simpson.index","total # clones","unique # clones")
a2

both <- cbind(a2,df_name)
both$Indiv_group <- paste(both$V1,both$V2,sep = "_")
both
