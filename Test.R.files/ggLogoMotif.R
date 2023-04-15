dataframe = read.csv("~/Desktop/TCR_Explore data_KM220413/All data/All paired TCR_210422_KM.csv")
dataframe <- read.csv("test-data/Group/paired_unsummarised2021.09.22.csv")
head(dataframe)
# 
# dataframe = read.csv("test-data/Group/paired_unsummarised2021.09.22.csv")

df.names <-   dataframe[ , -which(names(dataframe) %in% c("cloneCount","clone","Sequence_A","Sequence_B"))]

names(df.names) %in% c("Indiv","group","AVJ_aCDR3_BVDJ_bCDR3")

df1 <- ddply(dataframe,c("Indiv","group","AVJ_aCDR3_BVJ_bCDR3"),numcolwise(sum))

df1 <- df1[order(df1$cloneCount, decreasing = T),]

names(df1)

V1 <- df1[names(df1) %in% "Indiv"]
V1
V2 <- df1[names(df1) %in% "group"]
V1V2 <- cbind(V1,V2)
names(V1V2) <- c("V1","V2")

df1$selected.groups <- paste(V1V2$V1,V1V2$V2,sep="_")

df.group <- unique(df1[names(df1) %in% "selected.groups"])
names(df.group) <- "V1"

column.length <- length(df.group$V1)
column.length

df.group2 <- unique(df1[names(df1) %in% "AVJ_aCDR3_BVJ_bCDR3"])
names(df.group2) <- "V1"

row.length <- length(df.group2$V1)
row.length

m = matrix(NA,ncol=column.length, nrow=row.length)
samps <- df.group$V1

for (j in 1:column.length){
  df2 <- subset(df1,df1$selected.groups==samps[j])
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
a
f <- a/d
a1 <- rbind(a,b,d,f)  
a1 <- as.data.frame(a1)
names(a1) <- names(m)
a1
a1 <- as.data.frame(a1)
names(a1) <- names(m)
a1
df_name <- as.data.frame(do.call(rbind, strsplit(as.character(names(m)), "_")))
df_name$both <- paste(df_name$V1,df_name$V2,sep = "_")
head(df_name) 
names(df_name) <- c("Patient_ID","Drug",,"both")
a2 <- as.data.frame(t(a1))
names(a2) <- c("inv.simpson.div.index","total # clones","unique # clones","norm.inv.simpson.div.index")
a2

both <- cbind(a2,df_name)

both
dat <- both
conf <- 0.95
dat

dat <- dat[order(dat$both),]
dat
ve <- ifelse(input$varequal == 'y', TRUE, FALSE)
pair_samp <- ifelse(input$paired == 'y', TRUE, FALSE)
group1 <- subset(dat, dat$Patient_ID=="CD8") # group 1
group2 <- subset(dat, dat$Patient_ID=="IFN") # group 2

group1
group2
t.test(group1$inv.simpson.div.index, group2$inv.simpson.div.index, paired = T, alternative = "two.sided",conf.level = 0.95)
t.test(group1$norm.inv.simpson.div.index, group2$norm.inv.simpson.div.index, paired = T, alternative = "two.sided",conf.level = 0.95)

devtools::install_github("omarwagih/ggseqlogo")
require("ggseqlogo")
?ggseqlogo
ggseqlogo(mat, method='custom', seq_type='aa') + 
  ylab('JS divergence') + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) 


