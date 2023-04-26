require(scRepertoire)
require(ggalluvial)
# 
# 
compareClonotypes
dat = read.csv("test-data/Group/paired_TCR_file2022.05.24.csv",header=T) 

head(dat)

df1 <- dat[names(dat) %in% c("cloneCount","Indiv.group","AV", "AJ","BV","BJ","BD","CDR3.IMGT..AA._B","CDR3.IMGT..AA._A")]
df2 <- as.data.frame(ddply(dat,names(df1)[-c(1)],numcolwise(sum)))
df2_sub <- subset(df2,df2$Indiv.group=="E10630.IFN")
ggplot(data = df2,
       aes(axis1 = BV, axis2 = BD, axis3 = BJ,
           y = cloneCount)) +
  scale_x_discrete(limits = c("BV", "BD", "BJ"), expand = c(.2, .05)) +
  xlab("") +
  geom_alluvium(aes(fill = CDR3.IMGT..AA._B)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() 


ggplot(data = df2,
       aes(axis1 = BV, axis2 = BD, axis3 = BJ, axis4 = AV, axis5 = AJ,
           y = cloneCount)) +
  scale_x_discrete(limits = c("BV", "BD", "BJ","AV","AJ"), expand = c(.2, .05)) +
  xlab("") +
  geom_alluvium(aes(fill = CDR3.IMGT..AA._B)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() 

titanic_long <- to_lodes_form(data.frame(Titanic),
                              key = "Demographic",
                              axes = 1:3)
head(titanic_long)
ggplot(data = titanic_long,
       aes(x = Demographic, stratum = Demographic, alluvium = alluvium,
           y = Freq, label = Demographic)) +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")
