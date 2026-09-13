md <- read.csv("~/Desktop/Meta_data_2026-02-03_noBAL_noLTR34.4.csv")

library(Matrix)
library(data.table)

dt <- as.data.table(md)
head(dt)

dt$clone_case <- paste0(dt$vdj_gene_cdr3_AG_BD," ",dt$Case)

mat <- dcast(
  dt,
  vdj_gene_cdr3_AG_BD ~ Case_Time_Status,
  value.var = "cloneCount",
  fun.aggregate = sum,
  fill = 0
)

mat_df <- as.data.frame(mat)
rownames(mat_df) <- mat_df$vdj_gene_cdr3_AG_BD
head(mat_df)
mat_df$vdj_gene_cdr3_AG_BD <- NULL

mat <- as.matrix(mat_df)
mat
library(data.table)
meta <- data.frame(sample_id = colnames(mat))
setDT(meta)
meta[, event_code := tstrsplit(sample_id, "-", fixed=TRUE)[[2]]]
meta[, patient_time := tstrsplit(sample_id, "-", fixed=TRUE)[[1]]]
meta[, patient := tstrsplit(patient_time, "\\.")[[1]]]
meta[, time_raw := tstrsplit(patient_time, "\\.")[[2]]]
meta[, timepoint := as.numeric(gsub("_", ".", time_raw))]
meta[, event := as.numeric(gsub("A", "", event_code))]

meta[, has_event := any(event != 0), by = patient]
meta[, active_event := any(event != 0),by = event_code]
meta


#####


meta_time <- meta[has_event == TRUE]
mat_time <- mat[, meta_time$sample_id]
design <- model.matrix(~active_event+patient, data = meta_time)
design
dge <- DGEList(counts = mat_time)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)



res_time <- glmLRT(fit,coef = "active_eventTRUE")
DEG <- res_time$table
DEG$ID <- paste0("clone_",1:dim(DEG)[1])
DEG <- DEG[order(DEG$PValue),]

DEG

DEG <- DEG[c(1,4,5)]
head(DEG)
DEG$Case <- sapply(rownames(DEG), function(x) {
  paste(unique(dt[dt$vdj_gene_cdr3_AG_BD == x, "Case"]), collapse = ";")
})

DEG$ID <- paste0(DEG$ID,"-",DEG$Case)
DEG <- DEG[c(1,2,3)]
names(DEG) <- c("logFC","Pvalue","ID")
DEG$TCR <- 
write.csv(DEG,"~/Downloads/DEG_active_event.csv",row.names = T)


## alternative stratergy
clone_presence <- rowSums(mat > 0)
clone_presence
keep <- clone_presence >= 3   # or 5 as you said
keep
event_samples <- meta$event != 0

event_freq <- rowSums(mat[, event_samples] > 0)
non_event_freq <- rowSums(mat[, !event_samples] > 0)

score <- event_freq - non_event_freq

score <- as.data.frame(score)
subset(score,score>1)
summary(score)
library(stats)

mat



