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
