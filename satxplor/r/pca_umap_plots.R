options(conflicts.policy = list(warn = FALSE))
library(umap, quietly = TRUE)
library(optparse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(RSpectra, quietly = TRUE)
library(ape, quietly = TRUE)
library(stringr, quietly = TRUE)
library(Biostrings, quietly = TRUE)
library(ggplot2, quietly = TRUE)

# Define command-line options
option_list <- list(
  make_option(c("-a", "--sequence_alignment"),
  type = "character", help = "Sequence alignment file"),
  make_option(c("-d", "--dim_red"),
  type = "character", help = "Dimensionality reduction option ('umap' or 'pca' or 'both')"),
  make_option(c("-o", "--output_path"),
  type = "character", help = "Output files path"),
  make_option(c("-m", "--save_matrix"), default = TRUE,
  type = "logical", help = "Path to save matrices, if empty matrices wont be saved.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(opt)

if (is.null(opt$sequence_alignment) || is.null(opt$dim_red) || is.null(opt$output_path)) {
  cat("Usage: Rscript dimensionality_reduction_optparse.R --sequence_alignment <file> --dimensionality_reduction <umap_or_pca> --output_file <file>\n")
  quit(save = "no", status = 1)
}

msa <- readDNAMultipleAlignment(opt$sequence_alignment, format = "fasta")

first_id <- str_remove(rownames(msa)[1], "^_R_")
filename <- str_match(first_id, "^.*_([^_]+)_\\d+$")[, 2]
if (is.na(filename)) {
  filename <- str_split(first_id, "_", simplify = TRUE)[, 2]
}

x <- dist.dna(as.DNAbin(msa), model = "F81", as.matrix = TRUE, pairwise.deletion = TRUE)

cat("Finished distance calculating: \n", opt$sequence_alignment)

matrix <- as.data.table(x)

if (opt$save_matrix == TRUE) {
  fwrite(matrix, paste0(opt$output_path, filename, "_aligned_matrix.csv.gz"))
}

if (nrow(x) < 3) {
  cat("Fewer than three sequences for", filename, "- skipping PCA/UMAP\n")
  quit(save = "no", status = 0)
}

# FactoMineR::PCA computes a full N x N SVD, which is O(N^3) and was the main
# cost at 10k monomers.  Only the first components are used by the plots and
# UMAP, so compute those with RSpectra (O(N^2 k)).  The scaling matches
# FactoMineR::PCA(scale.unit = TRUE): centre and scale each column, then SVD.
partial_pca <- function(x, ncp) {
  X <- as.matrix(x)
  n <- nrow(X)

  # dist.dna(pairwise.deletion = TRUE) returns NA for pairs with no
  # overlapping non-gap sites.  FactoMineR::PCA used to impute these by the
  # mean of the variable, but a full PCA is O(N^3).  Reproduce that
  # imputation and continue with the partial SVD.
  if (anyNA(X)) {
    cm0 <- colMeans(X, na.rm = TRUE)
    cm0[!is.finite(cm0)] <- 0
    na_idx <- which(is.na(X), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      X[na_idx] <- cm0[na_idx[, "col"]]
    }
  }

  cm <- colMeans(X)
  csd <- apply(X, 2, sd)
  csd[!is.finite(csd) | csd == 0] <- 1
  Xs <- sweep(X, 2, cm, "-")
  Xs <- sweep(Xs, 2, csd, "/")

  k <- min(ncp, min(dim(Xs)) - 1)
  if (k < 1) {
    stop("Fewer than two sequences; cannot compute PCA")
  }
  if (min(dim(Xs)) < 2 * k) {
    s <- svd(Xs, nu = k, nv = k)
    # base svd returns min(n, p) singular values but only k vectors
    d <- s$d[seq_len(k)]
    v <- s$v[, seq_len(k), drop = FALSE]
  } else {
    s <- RSpectra::svds(Xs, k = k, nu = k, nv = k)
    d <- s$d
    v <- s$v
  }
  ord <- order(d, decreasing = TRUE)
  d <- d[ord]
  v <- v[, ord, drop = FALSE]

  var_coord <- v %*% diag(d, nrow = k)
  colnames(var_coord) <- paste0("Dim.", seq_len(k))
  rownames(var_coord) <- colnames(X)

  total_var <- sum(Xs^2) / (n - 1)
  eig_val <- (d^2) / (n - 1)
  pct <- eig_val / total_var * 100
  eig <- cbind(eigenvalue = eig_val,
               "percentage of variance" = pct,
               "cumulative percentage of variance" = cumsum(pct))
  list(var = list(coord = var_coord), eig = eig)
}

make_pca_plots <- function(pca_res, filename, output_path) {
  dt <- pca_res$var$coord %>% as.data.table(keep.rownames = TRUE)
  colnames(dt)[1] <- "rn"

  eigenvalues <- pca_res$eig[, 2]
  scree_data <- data.table(
    Principal_Component = seq_along(eigenvalues),
    Variance_Explained = eigenvalues
  )

  keep_n <- min(10, nrow(scree_data))
  p1 <- ggplot(scree_data[seq_len(keep_n)],
               aes(x = as.factor(Principal_Component), y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7, color = "black") +
    geom_text(aes(label = sprintf("%.2f", Variance_Explained)), vjust = -0.5) +
    labs(x = "Principal Component", y = "Variance Explained", title = "Scree Plot") +
    theme_minimal()

  ggsave(filename = paste0(output_path, filename, "_monomers_PCA_scree_plot.png"),
         plot = p1)

  dt[, ar_id := str_remove(rownames(msa), "^_R_")]
  dt[, origin := str_split(ar_id, "_", simplify = TRUE)[, 1]]

  p2 <- ggplot(dt) + geom_point(aes(x = Dim.1, y = Dim.2, color = origin), alpha = 0.4) +
    theme_bw() + theme(legend.position = "none")
  ggsave(plot = p2, filename = paste0(output_path, filename, "_monomers_PCA_plot.png"))
}

make_umap_plot <- function(umap_input, filename, output_path, msa) {
  set.seed(42)
  umap_res <- umap(umap_input)
  cat("Finished UMAP calculating: \n", opt$sequence_alignment)

  pca_r_dt <- umap_res$layout %>% as.data.table()
  pca_r_dt[, ar_id := str_remove(rownames(msa), "^_R_")]
  pca_r_dt[, origin := str_split(ar_id, "_", simplify = TRUE)[, 1]]

  p1 <- ggplot(pca_r_dt) +
    geom_point(aes(x = V1, y = V2, color = origin), alpha = 0.4) +
    theme_bw() + theme(legend.position = "none")
  ggsave(plot = p1, filename = paste0(output_path, filename, "_monomers_UMAP_plot.png"))
}

umap_on_pca <- !identical(Sys.getenv("SATXPLOR_UMAP_ON_PCA", "1"), "0")

if (opt$dim_red == "pca") {
  pca_res <- partial_pca(x, 10)
  cat("Finished PCA calculating: \n", opt$sequence_alignment)
  make_pca_plots(pca_res, filename, opt$output_path)
}

if (opt$dim_red == "umap") {
  if (umap_on_pca) {
    pca50 <- partial_pca(x, 50)
    umap_input <- pca50$var$coord
  } else {
    umap_input <- setnafill(as.data.table(x), fill = 0)
  }
  make_umap_plot(umap_input, filename, opt$output_path, msa)
}

if (opt$dim_red == "both") {
  pca50 <- partial_pca(x, 50)
  cat("Finished PCA calculating: \n", opt$sequence_alignment)
  make_pca_plots(pca50, filename, opt$output_path)

  if (umap_on_pca) {
    umap_input <- pca50$var$coord
  } else {
    umap_input <- setnafill(as.data.table(x), fill = 0)
  }
  make_umap_plot(umap_input, filename, opt$output_path, msa)
}
