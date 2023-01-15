?rpois
a <- matrix(rpois(100,100),nrow=10) 

b <- matrix(rpois(100,10000),nrow=10) *-1
a2 <- rbind(a,b)
a3 <- a2

log(1/a2*-1,10)


df <- as.data.frame(
  replicate(10, sample(1:100,100))
)

a2 <- as.matrix(a2)
as.data.frame(ifelse(a2 < 0, 1/a2*-1, a2))


a2 <- as.data.frame(a2)
a2[a2< -10000] <- 0.0001
a2[a2< -9000] <- 0.0002
a2[a2< -8000] <- 0.0003
a2[a2< -7000] <- 0.0004
a2[a2< -6000] <- 0.0005
a2[a2< -5000] <- 0.0006
a2[a2< -4000] <- 0.0007
a2[a2< -3000] <- 0.0008
a2[a2< -2000] <- 0.0009
a2[a2< -1000] <- 0.0010
a2[a2< -900] <- 0.0020
a2[a2< -800] <- 0.0030
a2[a2< -700] <- 0.0040
a2[a2< -600] <- 0.0050
a2[a2< -500] <- 0.0060
a2[a2< -400] <- 0.0070
a2[a2< -300] <- 0.0080
a2[a2< -200] <- 0.0090
a2[a2< -190] <- 0.0091
a2[a2< -180] <- 0.0092
a2[a2< -170] <- 0.0093
a2[a2< -160] <- 0.0094
a2[a2< -150] <- 0.0095
a2[a2< -140] <- 0.0096
a2[a2< -130] <- 0.0097
a2[a2< -120] <- 0.0098
a2[a2< -110] <- 0.0099
a2[a2< -100] <- 0.0100
a2[a2< -90] <- 0.020
a2[a2< -80] <- 0.030
a2[a2< -70] <- 0.040
a2[a2< -60] <- 0.050
a2[a2< -50] <- 0.060
a2[a2< -40] <- 0.070
a2[a2< -30] <- 0.080
a2[a2< -20] <- 0.090
a2[a2< -10] <- 0.1
a2[a2< -9] <- 0.2
a2[a2< -8] <- 0.3
a2[a2< -7] <- 0.4
a2[a2< -6] <- 0.5
a2[a2< -5] <- 0.6
a2[a2< -4] <- 0.7
a2[a2< -3] <- 0.8
a2[a2< -2] <- 0.9
a2[a2< 0] <- 1
a2

df <- read.csv("test-data/Index/TCR_Explore_index.clonal.2021.11.19.csv",header = T, fileEncoding = "UTF-8")
validate(
  need(nrow(df)>0,
       error_message_val2)
)
df <- as.data.frame(df)
head(df)


df2 <- df[,c("cloneCount",c("Indiv","group","TRBV","CDR3b.Sequence","TRBJ","TRAV","CDR3a.Sequence", "TRAJ","AJ", "BJ","AJBJ"))] 
df2
df3 <- as.data.frame(ddply(df2,input$string.data,numcolwise(sum)))
df3
df1 <- subset(df3,df3$cloneCount>input$numeric.cloneCount)
df1$clonal <- "yes"
colnames(df1)[which(names(df1) == "cloneCount")] <- "# of clones"
a <- df1[names(df1) %in% input$V.gene.1]
names(a) <- "V1"
b <- df1[names(df1) %in% input$CDR3.1]
names(b) <- "V1"
group.CDR <- df1[names(df1) %in% input$group.col.dot]
names(group.CDR) <- "V1"
d <- df1[names(df1) %in% input$V.gene.2]
names(d) <- "V1"
e  <- df1[names(df1) %in% input$CDR3.2]
names(e) <- "V1"
df1$gene.CDR3.1 <- paste(a$V1,b$V1,group.CDR$V1,sep="_")
df1$gene.CDR3.2 <- paste(d$V1,e$V1,group.CDR$V1,sep = "_")
df1$gene.CDR3.both <- paste(a$V1,b$V1,d$V1,e$V1,group.CDR$V1,sep = "_")
df1
a2 <- merge(df,df1,by=input$string.data,all=T)

df.mat1 <- a2[ , which(names(a2) %in% c("cloneCount","Indiv","group","TRBV","CDR3b.Sequence","TRBJ","TRAV","CDR3a.Sequence", "TRAJ","AJ", "BJ","AJBJ","well",    "Vb",  "Jb","TRBD",  "Va","Ja"))]
df.mat1
df.mat <- a2[ , -which(names(a2) %in% c("cloneCount","Indiv","group","TRBV","CDR3b.Sequence","TRBJ","TRAV","CDR3a.Sequence", "TRAJ","AJ", "BJ","AJBJ","well",    "Vb",  "Jb","TRBD",  "Va","Ja"))]


a2
df.mat <- as.matrix(df.mat)
df.mat <- as.data.frame(ifelse(df.mat < 0, 1/df.mat*-1, df.mat))
cbind(df.mat1,df.mat)
