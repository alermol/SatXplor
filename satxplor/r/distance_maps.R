suppressPackageStartupMessages({
  require(pheatmap)
  require("ComplexHeatmap")
  require(circlize)
  require(ape)
  require(data.table)
  require(scales)
  require(optparse)
  require(dplyr)
  require(stringr)
})

ht_opt$message <- FALSE

option_list <- list(
  make_option(c("-a", "--distance_matrix"), type = "character",
              help = "Binary float64 Mash distance matrix (column-major)"),
  make_option(c("-o", "--output_path"), type = "character", help = "Output files path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# The matrix is written by satxplor/kmer_distances.py as raw little-endian
# float64 values in R's column-major order.  Its size determines the order.
size <- file.info(opt$distance_matrix)$size
if (is.na(size) || size <= 0 || size %% 8 != 0) {
  stop("Distance matrix is missing or not a raw float64 matrix")
}
n <- as.integer(round(sqrt(size / 8)))
if (as.numeric(n) * as.numeric(n) * 8 != size) {
  stop("Distance matrix is not square")
}
if (n < 2) {
  warning("Fewer than two sequences in distance matrix; skipping heatmap")
  quit(save = "no", status = 0)
}

values <- readBin(opt$distance_matrix, what = "double", n = n * n,
                  size = 8, endian = "little")
dim(values) <- c(n, n)
x <- values

# Mash can report non-finite values for sketches with no shared hashes.
x[!is.finite(x)] <- NA
if (all(is.na(x))) {
  stop("Distance matrix contains no finite values")
}

# Same transformation as the original script: scale by the maximum distance,
# then rescale to [0, 1].  Mash distances are already roughly in [0, 1].
max_x <- max(x, na.rm = TRUE)
if (!is.finite(max_x) || max_x <= 0) max_x <- 1
x <- x / max_x

min_x <- min(x, na.rm = TRUE)
max_x <- max(x, na.rm = TRUE)
if (is.finite(min_x) && is.finite(max_x) && max_x > min_x) {
  x <- (x - min_x) / (max_x - min_x)
}
x[is.na(x)] <- min(x, na.rm = TRUE)

col_fun <- colorRamp2(c(0, 0.5, 1), c("#AD3212", "#F4ED7E", "#1A5276"))

# Keep clustering, but reuse the distance matrix already in memory instead of
# letting ComplexHeatmap call dist() on the N x N matrix (O(N^3)).
hc <- hclust(as.dist(x), method = "complete")
rm(values)

filename <- basename(opt$distance_matrix) %>%
  str_remove("\\.kmer_dist\\.bin$") %>%
  str_remove("_flanks$")

# Rendering the full N x N matrix as a raster is O(N^2) in memory and time:
# a 19,648 x 19,648 matrix used 66 GB and did not finish.  The final PNG is
# only 480 px, so draw at most DISPLAY_MAX cells per side in the full hclust
# order.  Distances and clustering remain full-resolution.
DISPLAY_MAX <- 2000L
display_idx <- seq_len(n)
display_note <- NULL
if (n > DISPLAY_MAX) {
  display_idx <- hc$order[seq(1, n, length.out = DISPLAY_MAX)]
  x_display <- x[display_idx, display_idx, drop = FALSE]
  display_note <- paste0("display downsampled to ", DISPLAY_MAX,
                         " of ", n, " sequences (full hclust order)")
} else {
  x_display <- x
}

png(paste0(opt$output_path, filename, "_flank_distances.png"), width = 480, height = 480)
if (n > DISPLAY_MAX) {
  h1 <- Heatmap(1 - x_display,
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_dend = FALSE, show_column_dend = FALSE,
                show_row_names = FALSE, show_column_names = FALSE,
                use_raster = TRUE,
                col = col_fun, name = " ",
                column_title = display_note,
                column_title_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(
                  title = "similarity", at = c(0, 0.5, 1)
                ))
} else {
  h1 <- Heatmap(1 - x_display, show_column_names = FALSE, col = col_fun, name = " ",
                cluster_rows = hc, cluster_columns = hc,
                use_raster = TRUE,
                heatmap_legend_param = list(
                  title = "similarity", at = c(0, 0.5, 1)
                ))
}
draw(h1)
dev.off()
