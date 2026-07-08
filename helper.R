font <- fontHelper::register_fonts(which = "common")
font <- as.data.frame(font_families())
names(font) <- "Fonts"

test_fun <- function() {
  for (i in 1:15) {
    incProgress(1/15)
    sum(runif(1000000,0,1))
  }
}



ASN$cols <- colorset(alphabet="AA",
                     colorScheme="chemistry")

draw_colnames_rotate <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = .5, rot = 0, gp = gpar(...)) # rot = rotation for # degrees
  return(res)}
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Nucleotide <- function (Nucleotide, seqlength) {
  nt <- c("A", "C", "G", "T")
  spec.no <- nrow(Nucleotide)
  count <- mat.or.vec(nr = 4, nc = seqlength)
  for (i in 1:seqlength) {
    count[1, i] <- length(which(Nucleotide[, i + 2] == nt[1]))
    count[2, i] <- length(which(Nucleotide[, i + 2] == nt[2]))
    count[3, i] <- length(which(Nucleotide[, i + 2] == nt[3]))
    count[4, i] <- length(which(Nucleotide[, i + 2] == nt[4]))
    
  }
  rownames(count) <- nt
  return(count)
}

options(shiny.maxRequestSize=200*1024^2)

# 95% confidence interval
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
middle <- function(x) { r <- quantile(x, probs=c(0.25, 0.25, 0.5, 0.75, 0.75))
names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
r
}
o <- function(x) {
  subset(x, x < quantiles_95(x)[1] | quantiles_95(x)[5] < x)
}

graph_type <- c("histogram","density")
axis_density_group <- c("x-axis","y-axis")

angle <- c(0,90,180,270)

error_message_val1 <- "No data found"
error_message_val2 <- "Uploading file"
error_message_val3 <- "Upload clone file"
error_message_val4 <- "no own list found\n \nSuggest uploading file\nheaders=ID"

simp.index.names <- c("total # clones","unique # clones")


# ---------------------------------------------------------
# 1) PLAIN HELPER FUNCTIONS (not reactive — define outside server())
# ---------------------------------------------------------

#' Build the four per-link style matrices (thickness/border/type/alpha)
#' used to highlight a selected set of clones on the chord diagram.
build_link_style <- function(mat, selected, thickness, colour, lty,
                             alpha_selected, alpha_unselected) {
  
  lwd_mat <- mat
  lwd_mat[lwd_mat > 0] <- "x"
  lwd_mat[rownames(lwd_mat) %in% selected & lwd_mat == "x"] <- thickness
  lwd_mat[!rownames(lwd_mat) %in% selected & lwd_mat == "x"] <- 0
  lwd_mat[lwd_mat == 0] <- 1
  
  border_mat <- mat
  border_mat[border_mat > 0] <- 1
  border_mat[rownames(border_mat) %in% selected & border_mat == 1] <- colour
  border_mat[!rownames(border_mat) %in% selected & border_mat == 1] <- 0
  border_mat[border_mat == 0] <- NA
  
  lty_mat <- mat
  lty_mat[lty_mat > 0] <- lty
  
  alpha_mat <- mat
  alpha_mat[alpha_mat > 0] <- 1
  alpha_mat[rownames(alpha_mat) %in% selected & alpha_mat == 1] <- alpha_selected
  alpha_mat[!rownames(alpha_mat) %in% selected & alpha_mat == 1] <- alpha_unselected
  
  list(lwd = lwd_mat, border = border_mat, lty = lty_mat, alpha = alpha_mat)
}

#' Draw sector text labels on track 1 (the "with labels" variant).
draw_sector_labels <- function() {
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.par(track.margin = c(0, 0))
    xlim <- get.cell.meta.data("xlim")
    sector.index <- get.cell.meta.data("sector.index")
    theta <- circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    aa <- if (theta < 90 || theta > 270) c(0, 0.5) else c(1, 0.5)
    circos.text(x = mean(xlim), y = 0.1, labels = sector.index, facing = dd, adj = aa)
  }, bg.border = NA)
}

#' Create track 1 but draw no text (the "no labels" variant still
#' needs the track allocated for consistent spacing).
draw_empty_track <- function() {
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.par(track.margin = c(0, 0))
  }, bg.border = NA)
}

#' Plain chord diagram with sector labels, uniform transparency.
plot_chord_labelled <- function(mat, grid.col, order, transparency) {
  circos.clear()                 # reset any leftover state before starting
  on.exit(circos.clear())        # guaranteed cleanup after this call, even on error
  circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
  chordDiagram(
    mat, annotationTrack = "grid", grid.col = grid.col,
    order = order, transparency = transparency,
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat)))))
  )
  draw_sector_labels()
}

#' Chord diagram with selected clones highlighted (thicker/coloured
#' links), with or without sector labels.
plot_chord_highlighted <- function(mat, grid.col, order, style, show_labels) {
  circos.clear()
  circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
  chordDiagram(
    mat, annotationTrack = "grid", grid.col = grid.col, order = order,
    link.lty = style$lty, link.lwd = style$lwd, link.border = style$border,
    # NOTE: original code computed `alpha_mat` but never passed it in
    # (the `transparency = alpha_mat` line was commented out), so the
    # per-selected/unselected transparency had no visible effect.
    # Uncomment the line below to actually apply it:
    # transparency = style$alpha,
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat)))))
  )
  if (show_labels) draw_sector_labels() else draw_empty_track()
}

#' Bare chord diagram — no labels, no highlighting.
plot_chord_plain <- function(mat, grid.col, order) {
  circos.clear()
  chordDiagram(
    mat, annotationTrack = "grid", grid.col = grid.col, order = order,
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat)))))
  )
}
