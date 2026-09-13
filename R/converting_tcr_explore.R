library(shiny)
library(dplyr)
library(DT)

options(shiny.maxRequestSize = 500 * 1024^2)

# ── UI helpers ────────────────────────────────────────────────────────────────

section <- function(title, ...) {
  tags$details(
    style = "margin-bottom:6px;",
    tags$summary(
      style = paste(
        "cursor:pointer; font-size:12px; font-weight:600;",
        "color:#555; padding:6px 0; list-style:none;",
        "border-top:1px solid #e0e0e0; user-select:none;"
      ),
      title
    ),
    tags$div(style = "padding:8px 0 4px 0;", ...)
  )
}

sb_label <- function(txt) {
  tags$label(
    style = "font-size:11px; color:#777; display:block; margin-bottom:2px;",
    txt
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  
  tags$head(tags$style(HTML("

    /* Page */
    body { background:#f7f8fa; font-family:'Helvetica Neue',Arial,sans-serif; }

    /* Sidebar */
    .well {
      background:#ffffff !important;
      border:1px solid #e4e4e4 !important;
      border-radius:8px !important;
      box-shadow:none !important;
      padding:16px !important;
    }

    .sidebar-title {
      font-size:13px; font-weight:700; color:#222;
      margin-bottom:14px; padding-bottom:8px;
      border-bottom:2px solid #4a90d9;
    }

    /* Tighten sidebar inputs */
    .well .form-group { margin-bottom:8px; }
    .well label { font-size:12px; color:#555; margin-bottom:2px; }
    .well .form-control { font-size:12px; height:30px; padding:4px 8px; }
    .well select[multiple] { height:auto; }
    .well .radio label,
    .well .checkbox label { font-size:12px; color:#555; }

    /* details/summary collapse */
    details > summary { outline:none; }
    details > summary::-webkit-details-marker { display:none; }
    details > summary::before {
      content:'▸ '; font-size:10px; color:#4a90d9;
    }
    details[open] > summary::before { content:'▾ '; }

    /* Download button */
    #downloadTABLE\\.Immunoseq {
      width:100%; margin-top:10px;
      background:#4a90d9; color:#fff;
      border:none; border-radius:5px;
      font-size:12px; padding:7px 0;
      font-weight:600;
    }
    #downloadTABLE\\.Immunoseq:hover { background:#357abd; }

    /* Main panel card */
    .main-card {
      background:#fff;
      border:1px solid #e4e4e4;
      border-radius:8px;
      padding:18px 20px;
    }

    /* Status bar */
    .status-bar {
      display:flex; align-items:center; gap:14px;
      padding:8px 12px; border-radius:6px;
      background:#f0f5ff; border:1px solid #c9d9f5;
      margin-bottom:14px; flex-wrap:wrap;
    }
    .status-chip {
      font-size:11px; font-weight:600; padding:2px 9px;
      border-radius:10px; background:#4a90d9; color:#fff;
    }
    .status-text { font-size:12px; color:#555; }

    /* Section heading in main */
    .main-section-head {
      font-size:11px; font-weight:700; color:#888;
      text-transform:uppercase; letter-spacing:.05em;
      margin:0 0 6px 0;
    }

    /* Page title */
    h2.shiny-title {
      font-size:18px; font-weight:700; color:#222;
      margin-bottom:18px;
    }

    /* DT table */
    .dataTables_wrapper { font-size:12px; }
    table.dataTable thead th {
      background:#f5f6f8; font-size:11px;
      font-weight:600; color:#444;
    }
  "))),
  
  titlePanel(
    tags$span(style="font-size:18px;font-weight:700;color:#222;",
              "TCR_Explore · File format converter")
  ),
  
  sidebarLayout(
    
    # ── Sidebar ───────────────────────────────────────────────────────────────
    sidebarPanel(
      width = 3,
      style = "overflow-y:auto; max-height:92vh; position:sticky; top:10px;",
      
      tags$div(class = "sidebar-title", "Import settings"),
      
      # 1. Input type
      sb_label("Input type"),
      selectInput("datasource", NULL,
                  choices = c("ImmunoSEQ", "MiXCR", "10x_scSeq",
                              "TIRDLE-seq", "Other"),
                  width = "100%"),
      
      # 2. File upload (always shown — no demo default)
      sb_label("Upload file (.tsv · .csv · .txt)"),
      fileInput("file_TSV.Immunoseq", NULL,
                accept = c(".tsv", ".csv", ".txt"),
                width  = "100%"),
      
      fluidRow(
        column(6,
               sb_label("Separator"),
               selectInput("sep.imm", NULL,
                           choices  = c(Tab = "\t", Comma = ",", Semicolon = ";"),
                           selected = "\t", width = "100%")
        ),
        column(6,
               sb_label("Quote"),
               selectInput("quote.imm", NULL,
                           choices  = c(`"` = '"', `'` = "'", None = ""),
                           selected = '"', width = "100%")
        )
      ),
      
      tags$hr(style="margin:10px 0;"),
      
      # 3. Sample metadata
      fluidRow(
        column(6,
               sb_label("Group ID"),
               textInput("group.imm", NULL, "Group", width = "100%")
        ),
        column(6,
               sb_label("Individual ID"),
               textInput("indiv.imm", NULL, "ID", width = "100%")
        )
      ),
      
      tags$hr(style="margin:10px 0;"),
      
      # 4. Column mapping — collapsible
      section("Column mapping",
              
              sb_label("Count column"),
              selectInput("countcolumn", NULL, choices = "", width = "100%"),
              
              conditionalPanel(
                condition = "input.datasource == 'ImmunoSEQ' || input.datasource == 'Other'",
                sb_label("In-frame column"),
                selectInput("inframe_immseq", NULL, choices = "", width = "100%")
              ),
              
              sb_label("CDR3 amino acid"),
              selectInput("CDR3.gene.clean", NULL, choices = "", width = "100%"),
              
              sb_label("V gene"),
              selectInput("V.GENE.clean", NULL, choices = "", width = "100%"),
              
              conditionalPanel(
                condition = "input.D_chain_present == 'yes'",
                sb_label("D gene"),
                selectInput("D.GENE.clean", NULL, choices = "", width = "100%")
              ),
              
              sb_label("J gene"),
              selectInput("J.GENE.clean", NULL, choices = "", width = "100%"),
              
              sb_label("D chain present?"),
              selectInput("D_chain_present", NULL,
                          choices = c("no", "yes"), selected = "no", width = "100%")
      ),
      
      # 5. Columns to remove — collapsible
      section("Columns to remove",
              selectInput("col.to.remove", NULL, choices = "",
                          multiple = TRUE, width = "100%",
                          selectize = FALSE, size = 6)
      ),
      
      downloadButton("downloadTABLE.Immunoseq", "Download CSV")
    ),
    
    # ── Main panel ────────────────────────────────────────────────────────────
    mainPanel(
      width = 9,
      
      tabsetPanel(
        
        tabPanel("Converted output",
                 tags$br(),
                 tags$div(class = "main-card",
                          
                          # Status bar
                          uiOutput("status_bar"),
                          
                          # Table
                          tags$p(class = "main-section-head", "Preview — converted file"),
                          DT::dataTableOutput("ImmunoSeq.table")
                 )
        ),
        
        tabPanel("Video walkthrough",
                 tags$br(),
                 uiOutput("video7")
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  # ── Data loading ─────────────────────────────────────────────────────────────
  
  input.data.Immunoseq <- reactive({
    inFile <- input$file_TSV.Immunoseq
    if (is.null(inFile)) return(NULL)
    read.table(inFile$datapath,
               sep    = input$sep.imm,
               quote  = input$quote.imm,
               header = TRUE)
  })
  
  TSV.col.names <- reactive({
    x <- as.data.frame(input.data.Immunoseq())
    x %>% select_if(~ !any(is.na(.)))
  })
  
  # ── Status bar ───────────────────────────────────────────────────────────────
  
  output$status_bar <- renderUI({
    src   <- input$datasource
    inFile <- input$file_TSV.Immunoseq
    
    chip_col <- switch(src,
                       "ImmunoSEQ"  = "#2e7d32",
                       "MiXCR"      = "#6a1b9a",
                       "10x_scSeq"  = "#bf360c",
                       "TIRDLE-seq" = "#e65100",
                       "#1565c0"
    )
    
    file_txt <- if (is.null(inFile)) {
      tags$span(class = "status-text",
                style = "color:#c0392b;",
                "\u26a0 No file uploaded")
    } else {
      nrow_txt <- tryCatch({
        df <- input.data.Immunoseq()
        paste0(format(nrow(df), big.mark = ","), " rows · ",
               ncol(df), " columns")
      }, error = function(e) "file loaded")
      tags$span(class = "status-text",
                paste0("\u2713 ", basename(inFile$name), " — ", nrow_txt))
    }
    
    tags$div(class = "status-bar",
             tags$span(class = "status-chip",
                       style = paste0("background:", chip_col, ";"),
                       src),
             file_txt,
             tags$span(class = "status-text",
                       style = "margin-left:auto; color:#888; font-size:11px;",
                       paste("Group:", input$group.imm,
                             "\u00b7 ID:", input$indiv.imm))
    )
  })
  
  # ── Auto column-mapping observers ────────────────────────────────────────────
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ"  = "count..templates.reads.",
                  "MiXCR"      = "cloneCount",
                  "10x_scSeq"  = "number_clonotypes",
                  "TIRDLE-seq" = "wij",
                  "duplicate_count"
    )
    updateSelectInput(session, "countcolumn", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ"  = "vFamilyName",
                  "MiXCR"      = "allVHitsWithScore",
                  "10x_scSeq"  = "v_gene_B",
                  "TIRDLE-seq" = "va",
                  "v_call"
    )
    updateSelectInput(session, "V.GENE.clean", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ"  = "dFamilyName",
                  "MiXCR"      = "allDHitsWithScore",
                  "10x_scSeq"  = "d_gene_B",
                  "TIRDLE-seq" = "",
                  "d_call"
    )
    updateSelectInput(session, "D.GENE.clean", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ"  = "jFamilyName",
                  "MiXCR"      = "allJHitsWithScore",
                  "10x_scSeq"  = "j_gene_B",
                  "TIRDLE-seq" = "ja",
                  "j_call"
    )
    updateSelectInput(session, "J.GENE.clean", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ"  = "aminoAcid",
                  "MiXCR"      = "aaSeqCDR3",
                  "10x_scSeq"  = "cdr3_amino_acid_B",
                  "TIRDLE-seq" = "cdr3a",
                  "junction_aa"
    )
    updateSelectInput(session, "CDR3.gene.clean", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ" = "sequenceStatus",
                  "Other"     = "vj_in_frame",
                  ""
    )
    updateSelectInput(session, "inframe_immseq", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  # ── Columns-to-remove observer ───────────────────────────────────────────────
  
  observe({
    cols <- names(TSV.col.names())
    sel <- switch(input$datasource,
                  "ImmunoSEQ"  = c("count..templates.reads.", "frequencyCount...."),
                  "TIRDLE-seq" = c("wi", "wj", "pval", "pval_adj", "r", "ts",
                                   "loss_a_frac", "loss_b_frac", "score", "wa", "wb"),
                  character(0)
    )
    updateSelectInput(session, "col.to.remove", choices = cols,
                      selected = intersect(sel, cols))
  })
  
  # ── Core processing reactive ──────────────────────────────────────────────────
  
  TSV.file.Immunoseq <- reactive({
    
    x <- as.data.frame(input.data.Immunoseq())
    validate(need(nrow(x) > 0, "Upload a file to begin"))
    
    x2 <- x %>% select_if(~ !any(is.na(.)))
    
    # ── ImmunoSEQ ──────────────────────────────────────────────────────────────
    if (input$datasource == "ImmunoSEQ") {
      
      x2 <- subset(x2, x2[[input$inframe_immseq]] == "In")
      x2 <- data.frame(cloneCount = x2[[input$countcolumn]], x2)
      names(x2)[1] <- "cloneCount"
      x2$group       <- input$group.imm
      x2$Indiv       <- input$indiv.imm
      x2$Indiv.group <- paste(x2$group, x2$Indiv, sep = ".")
      x3 <- x2
      
      x3$TRAV <- x3[[input$V.GENE.clean]]
      x3$TRAV[x3$TRAV == ""] <- "TRAV"
      x3$TRAJ <- x3[[input$J.GENE.clean]]
      x3$TRAJ[x3$TRAJ == ""] <- "TRAJ"
      
      if (input$D_chain_present == "yes") {
        x3$TRAD <- x3[[input$D.GENE.clean]]
        x3$TRAD[x3$TRAD == ""] <- "-"
      }
      
      x3 <- x3[!is.na(x3$TRAV), ]
      x3 <- x3[!is.na(x3$TRAJ), ]
      
      x3$TRAVJ      <- paste(x3$TRAV, x3$TRAJ, sep = ".")
      x3$TRAVJ_CDR3 <- paste(x3$TRAVJ, x3[[input$CDR3.gene.clean]], sep = "_")
      
      if (input$D_chain_present == "yes") {
        x3$TRAVDJ      <- paste(x3$TRAV, x3$TRAD, x3$TRAJ, sep = ".")
        x3$TRAVDJ_CDR3 <- paste(x3$TRAVDJ, x3[[input$CDR3.gene.clean]], sep = "_")
      }
      
      x3 <- x3[!names(x3) %in% input$col.to.remove]
      x3[x3 == ""] <- "Missing"
      x3
      
      # ── MiXCR ──────────────────────────────────────────────────────────────────
    } else if (input$datasource == "MiXCR") {
      
      x2$group       <- input$group.imm
      x2$Indiv       <- input$indiv.imm
      x2$group.indiv <- paste(x2$group, x2$Indiv, sep = ".")
      x2 <- x2 %>%
        select(all_of(c(input$countcolumn, "group", "Indiv", "group.indiv")),
               everything())
      names(x2)[1] <- "cloneCount"
      x3 <- x2
      
      df_v <- as.data.frame(t(as.data.frame(
        strsplit(x3[[input$V.GENE.clean]], "[*]"))))
      x3$TRAV <- df_v$V1
      
      df_j <- as.data.frame(t(as.data.frame(
        strsplit(x3[[input$J.GENE.clean]], "[*]"))))
      x3$TRAJ <- df_j$V1
      
      if (input$D_chain_present == "yes") {
        df_d <- as.data.frame(t(as.data.frame(
          strsplit(x3[[input$D.GENE.clean]], "[*]"))))
        x3$TRAD <- df_d$V1
      }
      
      cdr3_col <- input$CDR3.gene.clean
      if (any(grepl("[_]", x3[[cdr3_col]]))) x3 <- x3[!grepl("[_]", x3[[cdr3_col]]), ]
      if (any(grepl("[*]", x3[[cdr3_col]]))) x3 <- x3[!grepl("[*]", x3[[cdr3_col]]), ]
      
      x3$TRAVJ      <- paste(x3$TRAV, x3$TRAJ, sep = ".")
      x3$TRAVJ_CDR3 <- paste(x3$TRAVJ, x3[[cdr3_col]], sep = "_")
      
      if (input$D_chain_present == "yes") {
        x3$TRAVDJ      <- paste(x3$TRAV, x3$TRAD, x3$TRAJ, sep = ".")
        x3$TRAVDJ      <- gsub("[.]NA[.]", ".", x3$TRAVDJ)
        x3$TRAD        <- gsub("NA", "-", x3$TRAD)
        x3$TRAVDJ_CDR3 <- paste(x3$TRAVDJ, x3[[cdr3_col]], sep = "_")
      }
      
      if (any(is.na(x3$TRAV))) x3 <- x3[!is.na(x3$TRAV), ]
      if (any(is.na(x3$TRAJ))) x3 <- x3[!is.na(x3$TRAJ), ]
      
      x3 <- x3[!names(x3) %in% input$col.to.remove]
      x3
      
      # ── 10x scSeq ──────────────────────────────────────────────────────────────
    } else if (input$datasource == "10x_scSeq") {
      
      x2$cloneCount <- 1
      contigs <- subset(x2, productive == TRUE)
      
      drop_cols <- c("is_cell", "contig_id", "high_confidence",
                     "raw_consensus_id", "exact_subclonotype_id",
                     "umis", "reads", "length", "cdr3_nt",
                     grep("fwr|cdr1|cdr2", names(contigs), value = TRUE))
      contigs_lim <- contigs[!names(contigs) %in% drop_cols]
      
      contig_AG <- subset(contigs_lim, chain %in% c("TRA", "TRG"))
      ag_key    <- c(grep("gene", names(contig_AG), value = TRUE),
                     grep("cdr3", names(contig_AG), value = TRUE), "chain")
      contig_AG <- contig_AG %>% select(all_of(ag_key), everything())
      names(contig_AG)[seq_along(ag_key)] <-
        paste0(names(contig_AG)[seq_along(ag_key)], "_AG")
      
      contig_BD <- subset(contigs_lim, chain %in% c("TRB", "TRD"))
      bd_key    <- c(grep("gene", names(contig_BD), value = TRUE),
                     grep("cdr3", names(contig_BD), value = TRUE), "chain")
      contig_BD <- contig_BD %>% select(all_of(bd_key), everything())
      names(contig_BD)[seq_along(bd_key)] <-
        paste0(names(contig_BD)[seq_along(bd_key)], "_BD")
      
      contig_paired <- merge(contig_AG, contig_BD,
                             by  = c("barcode", "full_length", "productive",
                                     "raw_clonotype_id"),
                             all = TRUE)
      
      contig_paired$pairing <- ifelse(
        contig_paired$chain_BD == "TRB" & contig_paired$chain_AG == "TRA",
        "abTCR Paired",
        ifelse(contig_paired$chain_BD == "TRD" & contig_paired$chain_AG == "TRG",
               "gdTCR Paired", NA))
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% "d_gene_AG"]
      
      cp <- contig_paired
      cp$d_gene_BD <- sub("^$", "NA", cp$d_gene_BD)
      
      cp$vj_gene_AG  <- gsub("NA.NA", "",
                             paste(cp$v_gene_AG, cp$j_gene_AG, sep = "."))
      cp$vj_gene_BD  <- gsub("NA.NA", "",
                             gsub(".NA.", ".",
                                  paste(cp$v_gene_BD, cp$j_gene_BD, sep = ".")))
      cp$vdj_gene_BD <- gsub("NA.NA", "",
                             gsub(".NA.", ".",
                                  paste(cp$v_gene_BD, cp$d_gene_BD,
                                        cp$j_gene_BD, sep = ".")))
      
      cp$vj_gene_cdr3_AG  <- gsub("_NA", "",
                                  paste(cp$vj_gene_AG, cp$cdr3_AG, sep = "_"))
      cp$vj_gene_cdr3_BD  <- gsub("_NA", "",
                                  paste(cp$vj_gene_BD, cp$cdr3_BD, sep = "_"))
      cp$vdj_gene_cdr3_BD <- gsub("_NA", "",
                                  paste(cp$vdj_gene_BD, cp$cdr3_BD, sep = "_"))
      
      cp$vj_gene_AG_BD  <- paste(cp$vj_gene_AG, cp$vj_gene_BD, sep = " & ")
      
      cp$vdj_gene_AG_BD <- paste(cp$vj_gene_AG, cp$vdj_gene_BD, sep = " & ")
      cp$vdj_gene_AG_BD <- gsub("^ & |& $", "", cp$vdj_gene_AG_BD)
      
      cp$vj_gene_cdr3_AG_BD  <- paste(cp$vj_gene_cdr3_AG,
                                      cp$vj_gene_cdr3_BD, sep = " & ")
      cp$vj_gene_cdr3_AG_BD  <- gsub("^ & | & $", "", cp$vj_gene_cdr3_AG_BD)
      
      cp$vdj_gene_cdr3_AG_BD <- paste(cp$vj_gene_cdr3_AG,
                                      cp$vdj_gene_cdr3_BD, sep = " & ")
      cp$vdj_gene_cdr3_AG_BD <- gsub("^ & | & $", "", cp$vdj_gene_cdr3_AG_BD)
      
      names(cp)[names(cp) == "barcode"] <- "Cell_Index"
      cp <- cp[!duplicated(cp$Cell_Index), ]
      
      cp$group       <- input$group.imm
      cp$indiv       <- input$indiv.imm
      cp$Indiv.group <- paste(cp$group, cp$indiv, sep = ".")
      
      cp <- cp[!names(cp) %in% input$col.to.remove]
      cp
      
      # ── TIRDLE-seq ─────────────────────────────────────────────────────────────
    } else if (input$datasource == "TIRDLE-seq") {
      
      paired <- x2
      paired$ID <- input$indiv.imm
      
      paired <- paired[!duplicated(paired$alpha_beta), ]
      
      paired <- paired %>%
        dplyr::rename(cloneCount = wij,
                      v_gene_AG  = va,
                      j_gene_AG  = ja,
                      v_gene_BD  = vb,
                      j_gene_BD  = jb,
                      cdr3_AG    = cdr3a,
                      cdr3_BD    = cdr3b)
      
      stat_cols <- c("wi", "wj", "pval", "pval_adj", "r", "ts",
                     "loss_a_frac", "loss_b_frac", "score", "wa", "wb")
      paired <- paired[!names(paired) %in%
                         union(stat_cols, input$col.to.remove)]
      
      # Remove frameshifts (_) and stop codons (*) from both CDR3 sequences
      paired <- paired[!grepl("[_]", paired$cdr3_AG), ]
      paired <- paired[!grepl("[_]", paired$cdr3_BD), ]
      paired <- paired[!grepl("[*]", paired$cdr3_AG), ]
      paired <- paired[!grepl("[*]", paired$cdr3_BD), ]
      
      paired$vj_gene_AG <- gsub("NA.NA", "",
                                paste(paired$v_gene_AG, paired$j_gene_AG, sep = "."))
      paired$vj_gene_BD <- gsub("NA.NA", "",
                                paste(paired$v_gene_BD, paired$j_gene_BD, sep = "."))
      
      paired$vj_gene_cdr3_AG <- gsub("_NA", "",
                                     paste(paired$vj_gene_AG, paired$cdr3_AG, sep = "_"))
      paired$vj_gene_cdr3_BD <- gsub("_NA", "",
                                     paste(paired$vj_gene_BD, paired$cdr3_BD, sep = "_"))
      
      paired$vj_gene_AG_BD <- paste(paired$vj_gene_AG, paired$vj_gene_BD, sep = " & ")
      paired$vj_gene_AG_BD <- gsub("^ & | & $", "", paired$vj_gene_AG_BD)
      
      paired$vj_gene_cdr3_AG_BD <- paste(paired$vj_gene_cdr3_AG,
                                         paired$vj_gene_cdr3_BD, sep = " & ")
      paired$vj_gene_cdr3_AG_BD <- gsub("^ & | & $", "", paired$vj_gene_cdr3_AG_BD)
      
      paired$group       <- input$group.imm
      paired$Indiv       <- input$indiv.imm
      paired$Indiv.group <- paste(paired$group, paired$Indiv, sep = ".")
      
      paired
      
      # ── AIRR / Other ───────────────────────────────────────────────────────────
    } else {
      
      x2$InFrame <- tolower(as.character(x2[[input$inframe_immseq]]))
      if (nrow(x2[!grepl("false", x2$InFrame), ]) > 0)
        x2 <- x2[!grepl("false", x2$InFrame), ]
      
      if ("junction_aa" %in% names(x2)) {
        if (any(grepl("[*]", x2$junction_aa)))
          x2 <- x2[!grepl("[*]", x2$junction_aa), ]
        x2 <- x2[grepl(".", x2$junction_aa, fixed = TRUE), ]
        if (any(is.na(x2$junction_aa)))
          x2 <- x2[!is.na(x2$junction_aa), ]
      }
      
      x2 <- data.frame(cloneCount = x2[[input$countcolumn]], x2)
      names(x2)[1] <- "cloneCount"
      x3 <- x2
      x3 <- x3[!names(x3) %in% input$col.to.remove]
      
      x3$TRAV <- x3[[input$V.GENE.clean]]
      x3$TRAJ <- x3[[input$J.GENE.clean]]
      
      if (input$D_chain_present == "yes") {
        x3$TRAD <- x3[[input$D.GENE.clean]]
      }
      
      x3$TRAVJ      <- paste(x3$TRAV, x3$TRAJ, sep = ".")
      x3$TRAVJ_CDR3 <- paste(x3$TRAVJ, x3[[input$CDR3.gene.clean]], sep = "_")
      
      if (input$D_chain_present == "yes") {
        x3$TRAVDJ      <- paste(x3$TRAV, x3$TRAD, x3$TRAJ, sep = ".")
        x3$TRAVDJ      <- gsub("[.]NA[.]", ".", x3$TRAVDJ)
        x3$TRAD        <- gsub("NA", "-", x3$TRAD)
        x3$TRAVDJ_CDR3 <- paste(x3$TRAVDJ, x3[[input$CDR3.gene.clean]], sep = "_")
      }
      
      x3$group       <- input$group.imm
      x3$Indiv       <- input$indiv.imm
      x3$Indiv.group <- paste(x3$group, x3$Indiv, sep = ".")
      
      x3 <- x3[!names(x3) %in% c(input$col.to.remove, "cloneCount.1")]
      x3 <- subset(x3, x3$TRAV != "None")
      x3
    }
  })
  
  # ── Table output ──────────────────────────────────────────────────────────────
  
  output$ImmunoSeq.table <- DT::renderDataTable(
    escape  = FALSE,
    filter  = "top",
    options = list(
      lengthMenu = c(5, 10, 20, 50, 100),
      pageLength  = 10,
      scrollX     = TRUE,
      dom         = "lftip"
    ),
    { TSV.file.Immunoseq() }
  )
  
  # ── Download handler ──────────────────────────────────────────────────────────
  
  output$downloadTABLE.Immunoseq <- downloadHandler(
    filename = function() {
      paste0(input$group.imm, ".", input$indiv.imm,
             "_TCR_Explore.analysis.file.csv")
    },
    content = function(file) {
      write.csv(as.data.frame(TSV.file.Immunoseq()), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
