df <- read.csv("~/Desktop/Vd1 updated paired_TCR_file2022.05.12.csv")
df <- as.data.frame(df)
head(df)


df1 <-  df[names(df) %in% c("cloneCount",
                            "Indiv",
                            "group",
                            "GVJ",
                            "DVDJ",
                            "AVJ",
                            "BVDJ",
                            names(df[grep("JUNCTION",names(df))]),
                            names(df[grep("IMGT",names(df))])
)]  

df2 <- ddply(df1,names(df1[-c(1)]),numcolwise(sum))
df2

df3 <-  df2[names(df2) %in% c(names(df2[grep("JUNCTION",names(df2))]),
                            names(df2[grep("IMGT",names(df2))])
)]  

df3
for (i in 1:dim(df3)[1]) {
  df3[i,] <- nchar(df3[i,])
  df3
}
df3
names(df3) <- paste(names(df3),"length", sep="_")
df4 <- cbind(df2,df3)
df4
