suppressPackageStartupMessages({
  require(optparse)
  require(data.table)
  require(stringr)
  require(igraph)
  require(networkD3)
  require(htmlwidgets)
})

# Define command-line options
option_list <- list(
make_option(c("-m",  "--matrix"),
            type = "character",
            help = "Path to the matrix from the dimensionality reduction step"),
make_option(c("-a",  "--arrays"),
            type = "character",
            help = "Path to the matrix from the dimensionality reduction step"),
make_option(c("-k",  "--k_neigh"),
            type = "integer",
            help = "Number of nearest neighboursh"),
make_option(c("-o",  "--output_path"),
            type = "character",
            help = "Output files path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

mat_dt <- fread(opt$matrix)

first_id <- str_remove(names(mat_dt)[1], "^_R_")
sat_name <- str_match(first_id, "^.*_([^_]+)_\\d+$")[, 2]
if (is.na(sat_name)) {
  stop("Cannot parse satellite name from matrix column names")
}

ids <- names(mat_dt)
D <- as.matrix(mat_dt)
rownames(D) <- ids
colnames(D) <- ids
rm(mat_dt)
gc()

filename <- sat_name

arrays <- fread(opt$arrays)[V3 == sat_name]
if (nrow(arrays) < 2) {
  cat("Fewer than two arrays for", sat_name, "- skipping network\n")
  quit(save = "no", status = 0)
}

# Map monomer ids to array ids using the same overlap logic as before.
map_ids <- function(ids) {
  dt_synth <- data.table(id = ids)
  # Monomer ids are <contig>_<feature>_<start>.  Parse from the right so
  # contig names containing underscores are kept intact.
  dt_synth[, seqnames := str_match(id, "^(.*)_[^_]+_\\d+$")[, 2]]
  dt_synth[, start := as.double(str_extract(id, "\\d+$"))]
  dt_synth[, end := start + 1]
  setkey(dt_synth, seqnames, start, end)
  setkey(arrays, V1, V4, V5)
  foverlaps(dt_synth, arrays) %>%
    na.omit() %>%
    .[, .(id, array = paste0(seqnames, "_", V4))]
}

ids_clean <- str_remove(ids, "^_R_")
mapped <- map_ids(ids_clean)
lab <- mapped$array[match(ids_clean, mapped$id)]
keep <- !is.na(lab)
if (sum(keep) < 2) {
  cat("Fewer than two monomers could be assigned to arrays for", sat_name,
      "- skipping network\n")
  quit(save = "no", status = 0)
}
D <- D[keep, keep, drop = FALSE]
lab <- factor(lab[keep])

# Compute array x array mean distances directly.  This is numerically
# identical to melting the full N x N matrix, but never materialises N^2 rows.
if (anyNA(D)) {
  D0 <- D
  D0[is.na(D0)] <- 0
  good <- (!is.na(D)) * 1
  row_sums <- rowsum(D0, group = lab, reorder = TRUE)
  row_counts <- rowsum(good, group = lab, reorder = TRUE)
  rm(D0, good)
  gc()
  block_sums <- rowsum(t(row_sums), group = lab, reorder = TRUE)
  block_counts <- rowsum(t(row_counts), group = lab, reorder = TRUE)
  means <- block_sums / block_counts
} else {
  row_sums <- rowsum(D, group = lab, reorder = TRUE)
  block_sums <- rowsum(t(row_sums), group = lab, reorder = TRUE)
  counts <- as.numeric(table(lab))
  means <- block_sums / outer(counts, counts)
}
rm(D, row_sums); gc()

ua <- levels(lab)
m <- length(ua)
if (m < 2) {
  cat("Fewer than two arrays for", sat_name, "- skipping network\n")
  quit(save = "no", status = 0)
}
k <- min(opt$k_neigh, m - 1)

# For each var_array (row of `means`), take its k nearest name_arrays.
tmp_list <- vector("list", m)
for (i in seq_len(m)) {
  v <- means[i, ]
  v[i] <- Inf
  ord <- order(v, ua)[seq_len(k)]
  tmp_list[[i]] <- data.table(var_array = ua[i],
                              name_array = ua[ord],
                              mval = v[ord])
}
tmp <- rbindlist(tmp_list)

g <- graph_from_data_frame(tmp, directed = F)

p <- igraph_to_networkD3(g)

p$nodes$group <- str_remove(p$nodes$name, "_\\d+$")

graph <- forceNetwork(Links = p$links,  Nodes = p$nodes,  Source = "source",
                     Target = "target",  NodeID = "name",
                     Group = "group",  Value = "value",
                     zoom = TRUE,  linkDistance = 30,
                     linkWidth = 1,
                     arrows = FALSE,
                     charge = -50,
                     legend = TRUE,  opacity = 1)

saveWidget(graph, file = paste0(opt$output_path, filename, "_network.html"))
