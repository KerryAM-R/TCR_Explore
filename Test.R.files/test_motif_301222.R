devtools::install_github('wleepang/shiny-directory-input')

# motif -----
df <- read.csv("test-data/Group/paired_unsummarised2021.09.22.csv")
validate(
  need(nrow(df)>0,
       error_message_val1)
)
df <- as.data.frame(df)
names(df)
df_unique2 <- as.data.frame(unique (nchar(df[,names(df) %in% "JUNCTION..AA._B"])))
df_unique2
names(df_unique2) <- "len"

df_unique2$len<- df_unique2[order(df_unique2$len),]
df_unique2


df_unique <- nchar(df_unique[,names(df_unique) %in%" JUNCTION..AA._B"])
df_unique$Unique_clones <- 1
df_unique2 <- as.data.frame(unique(df_unique$len1))
names(df_unique2) <- "len"



subset(df_unique,df_unique$len1==c(15,14))

df_subset <- subset(df_unique,df_unique$len1==15)
df_subset <- subset(df_subset,df_subset$group=="CD8")

motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep("JUNCTION..AA._B",names(df_subset))], ""))))
motif
cbind(x=1,y=2,motif)


motif_count <- aa.count.function(cbind(x=1,y=2,motif), 15)
motif_count<-pcm2pfm(motif_count)
motif_count
motif<-new("pfm", mat=motif_count, name="",
           color=colorset(alphabet="AA",
                          colorScheme="chemistry"))
motif



motif_count1_aa <- motif

df_subset <- subset(df_unique,df_unique$len1==15)
df_subset <- subset(df_subset,df_subset$group=="IFN")

motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep("JUNCTION..AA._B",names(df_subset))], ""))))
motif
cbind(x=1,y=2,motif)


motif_count <- aa.count.function(cbind(x=1,y=2,motif), 15)
motif_count<-pcm2pfm(motif_count)
motif_count
motif<-new("pfm", mat=motif_count, name="",
           color=colorset(alphabet="AA",
                          colorScheme="chemistry"))



motif_count2_aa <- motif
motif_count2_aa

diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1_aa@mat), 
                                   pwm2 = as.data.frame(motif_count2_aa@mat), 
                                   alphabet = ASN
                                   
)

tail(diffLogoObj$letters$x)

diffLogo(diffLogoObj, diffLogoConfiguration = list(showSequenceLogosTop=T))
df_subset <- subset(df_unique,df_unique$len1==15)
                            
x <- createDiffLogoObject(as.data.frame(motif_count1_aa@mat), as.data.frame(motif_count2_aa@mat), 
                     alphabet = ASN, baseDistribution = differenceOfICs,
                     )
                     


mat <- (x$pwm1 - x$pwm2)

require(ggseqlogo)
names(mat) <- 1:dim(mat)[2]
mat
ggseqlogo(mat, seq_type='aa')
ggseqlogo(mat, seq_type='aa') + 
  ylab('JS diversity')+
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  annotate(geom="text",x=1,y=Inf,vjust=2,label="IFN",size=10,face="plain",family="serif")+
  annotate(geom="text",x=1,y=-Inf,vjust=-2,label="CD8",size=10,face="plain",family="serif")+
theme(
      axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
      axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
      legend.title  =element_blank(),
      legend.position = "bottom",
      legend.text = element_text(colour="black", size=8,family="serif"),
      ) 

ggseqlogo(mat, method='custom', seq_type='aa') + 
  ylab('JS diversity')+ 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme(
    axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
    axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
    axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
    axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
    legend.title  =element_blank(),
    legend.position = "bottom",
    legend.text = element_text(colour="black", size=8,family="serif")) +
  

  
  motif_count<-pcm2pfm(motif_count)
motif_count
motif<-new("pfm", mat=motif_count, name="",
           color=colorset(alphabet="AA",
                          colorScheme="chemistry"))
motif  
  


longData<-longData[longData$value!=0,]


motif<-new("pfm", mat=mat, name="",
           color=colorset(alphabet="AA",
                          colorScheme="chemistry"))

motif@mat

#levels must be names
head(longData)
ggplot(longData, ggfortify(x = Var1, y = value)) + 
  geom_logo() 


data(sequences)

library(ggplot2)
library(gglogo)

df_subset <- subset(df_unique,df_unique$len1==15)
x <- ggfortify(df_subset, JUNCTION..AA._B,treatment = group)

df.tmp <- merge(df.g1,df.g2,by=c("element","position","Polarity","Water"),all=T)
df.tmp

df.tmp[is.na(df.tmp)] <- 0

df.tmp$result <- abs(df.tmp$bits.x - df.tmp$bits.y)

head(df.tmp,10)


type <- df.tmp[names(df.tmp) %in% c("element","Polarity","Water")]
type

unique(df.tmp)
type <- unique(type)
type
head(longData)
head(type)
df.1 <- merge(as.data.frame(longData),as.data.frame(type),by.x="Var1",by.y="element")

names(df.1)[1:2] <- c("element","position")


subset(x,x$position==5)

?geom_logo

ggplot(data = subset(x,x$group=="CD8")) +  
  geom_logo(aes(x=position, y=bits,  
                label=element, fill=interaction(Polarity, Water),
                ),
            inherit.aes = F,
            alpha = 1) 

scale_fill_brewer("Amino-acids properties",palette="Paired") +
  geom_hline(yintercept=0)+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid  = element_blank(),
        )+
  guides(color = "none", size = "none")+
  ylab("Shannon information in bits")


longData

ggplot(data = longData) +     
  geom_logo(aes(x=Var1, y=value, group=element, 
                label=element, fill=interaction(Polarity, Water)),
            alpha = 1)+
  
  scale_fill_brewer("Amino-acids properties",palette="Paired") +
  geom_hline(yintercept=0)+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid  = element_blank(),
  )+
  guides(color = "none", size = "none")+
  ylab("Shannon information in bits")

diffLogo(diffLogoObj, diffLogoConfiguration = list(showSequenceLogosTop=T))
