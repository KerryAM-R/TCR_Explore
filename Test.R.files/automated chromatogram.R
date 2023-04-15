dataframe_ab1 = read.csv("test-data/QC/SJS.TEN/AP026/AP026 IFN .ab1 QC 2022-12-31.csv") 
dataframe_ab1$pa_score_base

df_chain1 = read_excel("test-data/QC/SJS.TEN/AP026/vquest.xls") 

name_temp2  <- as.data.frame(do.call(rbind, strsplit(as.character(df_chain1$`Sequence ID`), ".seq")))
df_chain1$name_temp <- name_temp2$V1

df_chain1$name_temp
dataframe_ab1$name_temp

merge(df_chain1,dataframe_ab1, by = "name_temp")

head(dataframe_ab1)

dataframe_ab1$pa_score_base

ifelse(as.numeric(dataframe_ab1$pa_score_base) < 0,"Poor","other")

dataframe_ab1$chromatogram_check <- ifelse(dataframe_ab1$pa_score_base < 0.2,"Poor")
                                       ifelse(dataframe_ab1$pa_score_base > 1,"Very high",
                                              ifelse(dataframe_ab1$pa_score_base > 0.9, "High",
                                                     ifelse(dataframe_ab1$pa_score_base >0.7,"Moderate","Low")))
