# chao1, hill, div, gini.simp, inv.simp, gini, raref, d50, dxx.

require(vegan)
require("iNEXT")
require("fossil")
require (hillR)

require(rasterdiv)


df_chain1$clone_quality <- ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="Very high",'pass',
                                  ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="High",'pass',
                                         ifelse(df_chain1$V.sequence.quality.check=="No alignment" |
                                                  df_chain1$chromatogram_check=="Poor" |
                                                  df_chain1$V.sequence.quality.check=="No arrangement" |
                                                  df_chain1$chromatogram_check=="Low" ,'fail','fail')))

                                                

                                                      
                                                               
                                                                             ifelse(df_chain1$chromatogram_check=="High" & df_chain1$V.sequence.quality.check=="J Identity issue","pass",
                                                                                    ifelse(df_chain1$chromatogram_check=="High" & df_chain1$V.sequence.quality.check=="V Identity issue","pass",


                                                                                           "fail"
                                                #                                     )
                                                #                                     
                                                #                                     
                                                #                              )
                                                #                       )
                                                #                )
                                                #        )))))



dataframe = read.csv("test-data/Group/paired_TCR_file2022.05.24.csv",header=T) 

as.data.frame(do.call(rbind, strsplit(as.character(dataframe$`Sequence ID`), "_")))[1]

df1 <- ddply(dataframe,c("group","Indiv","AVJ_aCDR3_BVJ_bCDR3"),numcolwise(sum))
head(df1)
# df1 <- ddply(dataframe,c(input$group_column_simp,input$group_column_simp2,input$group_column_simp3),numcolwise(sum))
df1 <- df1[order(df1$cloneCount, decreasing = T),]

V1 <- df1[names(df1) %in% "group"]
V1
V2 <- df1[names(df1) %in% "Indiv"]
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
m
names(m) <- samps
head(m)
m

richness_samp <- matrix(nrow=1,ncol=dim(m)[2])
inverse_simp <- matrix(nrow=1,ncol=dim(m)[2])
shannon_div <- matrix(nrow=1,ncol=dim(m)[2])
chao1_div <- matrix(nrow=1,ncol=dim(m)[2])
Pielou_even <- matrix(nrow=1,ncol=dim(m)[2])
b <- matrix(nrow=1,ncol=dim(m)[2])
d <- matrix(nrow=1,ncol=dim(m)[2])

# evenness <- H/log(richness)


# samp <- m[,1]
# samp <- na.omit(samp)
# specnumber(samp)


# ?diversity
for( i in 1:dim(m)[2]) {
  
  samp <- m[,i]
  samp <- na.omit(samp)
  richness_samp[,i] <-specnumber(samp)
  inverse_simp[,i] <- diversity(samp,"invsimpson")
  shannon_div[,i] <- diversity(samp,"shannon")
  Pielou_even[,i] <- diversity(samp,"shannon")/log(specnumber(samp))
  chao1_div[,i] <-chao1(samp, taxa.row = TRUE)
  b[,i] <- sum(samp)
  d[,i] <- nrow(as.data.frame(samp))
}




a1 <- rbind(b, richness_samp, shannon_div, Pielou_even, inverse_simp,chao1_div)  
a1 <- as.data.frame(a1)
a1
names(a1) <- names(m)

df_name <- as.data.frame(do.call(rbind, strsplit(as.character(names(m)), "_")))
df_name$both <- paste(df_name$V1,df_name$V2,sep = "_")
head(df_name) 
names(df_name) <- c("group","Indiv","both")

a2 <- as.data.frame(t(a1))
names(a2) <- c("total # clones","richness_samp", "shannon_div", "Pielou_even", "inverse_simp","chao1_div")
both <- cbind(a2,df_name)
both$inv.simpson.index_div_unique.samp <- both$inv.simpson.index/both$`total # clones`
as.data.frame(both)

aov(Indiv ~ chao1_div, both)

oneway.test(group ~ chao1_div,
            data = both,
            var.equal = TRUE # assuming equal variances
)




both <- as.data.frame(both)

cols <- unlist(colors_inv.simp())

selected.col <- both[names(both) %in% "group"]
names(selected.col) <- "V1"
both[names(both) %in% "group"] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
head(both)

unique.col <- as.data.frame(unique(both[names(both) %in% input$group2.index]))
names(unique.col) <- "V1"
unique.col$simp.inv_palette <- cols
head (both)

ggplot(both,aes(x=group,y=Pielou_even))+
  stat_summary(fun.data = quantiles_95, geom = "errorbar", position = position_dodge(1)) +       
  stat_summary(fun.data = middle, geom = "boxplot", position = position_dodge(1))  


dat <- both
conf <- 0.95
dat <- dat[order(dat$both),]
ve <- T
pair_samp <- T
group1 <- subset(dat, dat$group=="CD8") # group 1
group2 <- subset(dat, dat$group=="IFN") # group 1
head(group2)
group1[,names(group1) %in% "shannon_div"]

t.test(group1[,names(group1) %in% "shannon_div"],group2[,names(group2) %in% "shannon_div"], paired = pair_samp, var.equal = ve,conf.level = conf)
