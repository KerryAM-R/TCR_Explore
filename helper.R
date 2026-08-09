options(
  shiny.maxRequestSize = 200 * 1024^2
)

font <- fontHelper::register_fonts(which = "common")
font <- as.data.frame(font)
names(font) <- "Fonts"


#####
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


###### 

DNA <- list(
  chars = c(
    "A", "C", "G", "T"
  ),
  cols = c(
    "green4",
    "blue",
    "orange",
    "red"
  ),
  size = 4,
  supportReverseComplement = TRUE
)

class(DNA) <- "Alphabet"

ASN <- list(
  chars = c(
    "A", "C", "D", "E", "F",
    "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R",
    "S", "T", "V", "W", "Y"
  ),
  
  cols = unname(
    motifStack::colorset(
      alphabet = "AA",
      colorScheme = "chemistry"
    )
  ),
  
  size = 20,
  
  supportReverseComplement = FALSE
)

class(ASN) <- "Alphabet"
ASN

test_fun <- function() {
  for (i in 1:15) {
    shiny::incProgress(1 / 15)
    sum(stats::runif(1000000, 0, 1))
  }
}

draw_colnames_rotate <- function(coln, gaps, ...) {
  coord <- pheatmap::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  
  res <- grid::textGrob(
    coln,
    x = x,
    y = grid::unit(1, "npc") - grid::unit(3, "bigpts"),
    vjust = 1,
    hjust = 0.5,
    rot = 0,
    gp = grid::gpar(...)
  )
  
  return(res)
}

gg_fill_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

Nucleotide <- function(Nucleotide, seqlength) {
  nt <- c("A", "C", "G", "T")
  spec.no <- nrow(Nucleotide)
  
  count <- base::mat.or.vec(
    nr = 4,
    nc = seqlength
  )
  
  for (i in 1:seqlength) {
    count[1, i] <- length(which(Nucleotide[, i + 2] == nt[1]))
    count[2, i] <- length(which(Nucleotide[, i + 2] == nt[2]))
    count[3, i] <- length(which(Nucleotide[, i + 2] == nt[3]))
    count[4, i] <- length(which(Nucleotide[, i + 2] == nt[4]))
  }
  
  rownames(count) <- nt
  
  return(count)
}



# 95% confidence interval

quantiles_95 <- function(x) {
  r <- stats::quantile(
    x,
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
  )
  
  names(r) <- c(
    "ymin",
    "lower",
    "middle",
    "upper",
    "ymax"
  )
  
  r
}

middle <- function(x) {
  r <- stats::quantile(
    x,
    probs = c(0.25, 0.25, 0.5, 0.75, 0.75)
  )
  
  names(r) <- c(
    "ymin",
    "lower",
    "middle",
    "upper",
    "ymax"
  )
  
  r
}

o <- function(x) {
  base::subset(
    x,
    x < quantiles_95(x)[1] |
      quantiles_95(x)[5] < x
  )
}

graph_type <- c("histogram", "density")
axis_density_group <- c("x-axis", "y-axis")

angle <- c(0, 90, 180, 270)

error_message_val1 <- "No data found"
error_message_val2 <- "Uploading file"
error_message_val3 <- "Upload clone file"
error_message_val4 <- "no own list found\n\nSuggest uploading file\nheaders=ID"

simp.index.names <- c(
  "total # clones",
  "unique # clones"
)


# ---------------------------------------------------------
# 1) PLAIN HELPER FUNCTIONS
# ---------------------------------------------------------

#' Build the four per-link style matrices
build_link_style <- function(
    mat,
    selected,
    thickness,
    colour,
    lty,
    alpha_selected,
    alpha_unselected
) {
  
  lwd_mat <- mat
  lwd_mat[lwd_mat > 0] <- "x"
  lwd_mat[
    rownames(lwd_mat) %in% selected &
      lwd_mat == "x"
  ] <- thickness
  
  lwd_mat[
    !rownames(lwd_mat) %in% selected &
      lwd_mat == "x"
  ] <- 0
  
  lwd_mat[lwd_mat == 0] <- 1
  
  
  border_mat <- mat
  border_mat[border_mat > 0] <- 1
  
  border_mat[
    rownames(border_mat) %in% selected &
      border_mat == 1
  ] <- colour
  
  border_mat[
    !rownames(border_mat) %in% selected &
      border_mat == 1
  ] <- 0
  
  border_mat[border_mat == 0] <- NA
  
  
  lty_mat <- mat
  lty_mat[lty_mat > 0] <- lty
  
  
  alpha_mat <- mat
  alpha_mat[alpha_mat > 0] <- 1
  
  alpha_mat[
    rownames(alpha_mat) %in% selected &
      alpha_mat == 1
  ] <- alpha_selected
  
  alpha_mat[
    !rownames(alpha_mat) %in% selected &
      alpha_mat == 1
  ] <- alpha_unselected
  
  
  list(
    lwd = lwd_mat,
    border = border_mat,
    lty = lty_mat,
    alpha = alpha_mat
  )
}


#' Draw sector text labels on track 1
draw_sector_labels <- function() {
  
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      
      circlize::circos.par(
        track.margin = c(0, 0)
      )
      
      xlim <- circlize::get.cell.meta.data("xlim")
      
      sector.index <- circlize::get.cell.meta.data(
        "sector.index"
      )
      
      theta <- circlize::circlize(
        mean(xlim),
        1.3
      )[1, 1] %% 360
      
      dd <- ifelse(
        theta < 90 || theta > 270,
        "clockwise",
        "reverse.clockwise"
      )
      
      aa <- if (
        theta < 90 || theta > 270
      ) {
        c(0, 0.5)
      } else {
        c(1, 0.5)
      }
      
      circlize::circos.text(
        x = mean(xlim),
        y = 0.1,
        labels = sector.index,
        facing = dd,
        adj = aa
      )
    },
    bg.border = NA
  )
}


#' Create track 1 but draw no text
draw_empty_track <- function() {
  
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      
      circlize::circos.par(
        track.margin = c(0, 0)
      )
    },
    bg.border = NA
  )
}


#' Plain chord diagram with sector labels
plot_chord_labelled <- function(
    mat,
    grid.col,
    order,
    transparency
) {
  
  circlize::circos.clear()
  
  on.exit(
    circlize::circos.clear()
  )
  
  circlize::circos.par(
    "canvas.xlim" = c(-1, 1),
    "canvas.ylim" = c(-1, 1)
  )
  
  circlize::chordDiagram(
    mat,
    annotationTrack = "grid",
    grid.col = grid.col,
    order = order,
    transparency = transparency,
    preAllocateTracks = list(
      track.height = max(
        graphics::strwidth(
          unlist(dimnames(mat))
        )
      )
    )
  )
  
  draw_sector_labels()
}


#' Chord diagram with selected clones highlighted
plot_chord_highlighted <- function(
    mat,
    grid.col,
    order,
    style,
    show_labels
) {
  
  circlize::circos.clear()
  
  circlize::circos.par(
    "canvas.xlim" = c(-1, 1),
    "canvas.ylim" = c(-1, 1)
  )
  
  circlize::chordDiagram(
    mat,
    annotationTrack = "grid",
    grid.col = grid.col,
    order = order,
    link.lty = style$lty,
    link.lwd = style$lwd,
    link.border = style$border,
    
    # transparency = style$alpha,
    
    preAllocateTracks = list(
      track.height = max(
        graphics::strwidth(
          unlist(dimnames(mat))
        )
      )
    )
  )
  
  if (show_labels) {
    draw_sector_labels()
  } else {
    draw_empty_track()
  }
}


#' Bare chord diagram -----
plot_chord_plain <- function(
    mat,
    grid.col,
    order
) {
  
  circlize::circos.clear()
  
  circlize::chordDiagram(
    mat,
    annotationTrack = "grid",
    grid.col = grid.col,
    order = order,
    preAllocateTracks = list(
      track.height = max(
        graphics::strwidth(
          unlist(dimnames(mat))
        )
      )
    )
  )
}


# ---------------------------------------------------------
# MOTIF PLOTTING HELPERS
# ---------------------------------------------------------

motif_theme <- function(input) {
  theme(
    axis.text.x = element_text(
      colour = "black",
      size = 20,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      face = "plain",
      family = input$font_type
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 20,
      angle = 0,
      hjust = 1,
      vjust = 0,
      face = "plain",
      family = input$font_type
    ),
    axis.title.x = element_text(
      colour = "black",
      size = 20,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      face = "plain",
      family = input$font_type
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 20,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      face = "plain",
      family = input$font_type
    ),
    legend.title = element_blank(),
    legend.position = input$legend_position,
    legend.text = element_text(
      colour = "black",
      size = input$legend_text_size,
      family = input$font_type
    )
  )
}


# ---------------------------------------------------------
# AMINO ACID MOTIF PLOT
# ---------------------------------------------------------

plot_aa_motif <- function(motif, input) {
  
  ggseqlogo::ggseqlogo(
    motif,
    seq_type = "aa",
    method = "p"
  ) +
    ylab("bits") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    motif_theme(input)
}


# ---------------------------------------------------------
# DNA MOTIF PLOT
# ---------------------------------------------------------

plot_dna_motif <- function(motif, input) {
  
  ggseqlogo::ggseqlogo(
    motif,
    seq_type = "dna",
    method = "p"
  ) +
    ylab("bits") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    motif_theme(input)
}


# ---------------------------------------------------------
# DIFFERENCE MOTIF PLOT
# ---------------------------------------------------------

plot_diff_motif <- function(motif1, motif2, alphabet, seq_type, input) {
  
  diffLogoObj <- DiffLogo::createDiffLogoObject(
    pwm1 = as.data.frame(motif1),
    pwm2 = as.data.frame(motif2),
    alphabet = alphabet
  )
  
  mat <- diffLogoObj$pwm1 - diffLogoObj$pwm2
  colnames(mat) <- seq_len(ncol(mat))
  
  ggseqlogo::ggseqlogo(
    mat,
    method = "custom",
    seq_type = seq_type
  ) +
    ylab("JS divergence") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    motif_theme(input)
}


# ---------------------------------------------------------
# GET GROUP-SPECIFIC AMINO ACID MOTIF
# ---------------------------------------------------------

get_group_aa_motif <- function(df, group, chain_col) {
  
  validate(
    need(
      nrow(df) > 0,
      error_message_val1
    )
  )
  
  df_group <- subset(
    df,
    df$group == group
  )
  
  validate(
    need(
      nrow(df_group) > 0,
      error_message_val1
    )
  )
  
  # Extract aligned amino-acid sequences
  sequences <- df_group[[chain_col]]
  
  motif <- as.data.frame(
    t(
      as.data.frame(
        strsplit(sequences, "")
      )
    )
  )
  
  # Alignment length
  seqlength <- ncol(motif)
  
  # Count amino acids
  motif_count <- VLF::aa.count.function(
    cbind(
      x = 1,
      y = 2,
      motif
    ),
    seqlength
  )
  
  # Convert counts to probabilities
  motifStack::pcm2pfm(motif_count)
}


