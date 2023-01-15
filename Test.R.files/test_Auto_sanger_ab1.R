

hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-A10_C07.ab1") 
hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-A2_A07.ab1") 
hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-A5_B07 low.ab1") 
hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-B4_D07.ab1") 
hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-B7_E07.ab1") 
hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-B9_F07.ab1") 
hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-C3_H07.ab1") 
hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
Primaryseq <- primarySeq(hetcalls, string = TRUE)
secondary_seq <- secondarySeq(hetcalls, string = TRUE)

pa <- pairwiseAlignment(primarySeq(hetcalls), secondarySeq(hetcalls))
pa@score/length(hetcalls@secondarySeq)[1]


chromatogram(hetcalls, width = 100, height = 2, showcalls = "both", trim5 =20, trim3 = 20)

inFile.seq <- input$file1_seq.file
  
  numfiles = nrow(inFile.seq) 
  
  
  df_total = data.frame()
  
  
  for (i in input$lower.seq:input$upper.seq) { 
    tryCatch({
      temp<- readsangerseq(input$file1_seq.file[[i, 'datapath']],header = F) 
      temp
      df_total = data.frame()
      
      hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
      
      print(hetcalls)
      Primaryseq <- primarySeq(hetcalls, string = TRUE)
      secondary_seq <- secondarySeq(hetcalls, string = TRUE)
      pa <- pairwiseAlignment(primarySeq(hetcalls), secondarySeq(hetcalls))
      sc <- pa@score
      name_temp <- paste(inFile.seq[i,1])
        
      
      
      
      

      
      
      df <- cbind(name_temp,sc)
      df_total <- rbind(df_total,df)
      
    }, error=function(e){}) 
  }
  df_total 
