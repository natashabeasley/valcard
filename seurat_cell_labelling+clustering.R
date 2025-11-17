# Load in FILTERED dataset
multiome_list_rm_doublets <- readRDS("/scratch/pigblast/nb443/ValCARD_seqdata.1/valcard_list_rm_doublets_min_100.rds")
# Filter out datasets with < 50 cells
multiome_list_filtered <- multiome_list_rm_doublets[sapply(multiome_list_rm_doublets, function(x) ncol(x) >= 50)]
# Merge all the objects in the list into one - processing multiome data
multiome <- merge(x = multiome_list_filtered[[1]], 
                  y = multiome_list_filtered[-1])
# Load the SCT + PCA'd heart reference
reference <- readRDS("/scratch/pigblast/nb443/heart_reference_with_sct.rds")
reference <- UpdateSeuratObject(reference)
DefaultAssay(reference) <- "SCT"
# Safety: ensure reference has PCA (and UMAP if you want to project query onto reference UMAP)
if (!"pca" %in% Reductions(reference)) {
  reference <- RunPCA(reference, verbose = FALSE)}
if (!"umap" %in% Reductions(reference)) {
  reference <- FindNeighbors(reference, dims = 1:50, verbose = FALSE)
  reference <- FindClusters(reference, resolution = 0.5, verbose = FALSE)
  reference <- RunUMAP(reference, dims = 1:50, reduction = "pca", verbose = FALSE)}
# Find anchors between reference (SCT) and query (SCT)
anchors <- FindTransferAnchors(
  reference = reference,
  query = multiome,
  normalization.method = "SCT",
  reference.reduction = "pca",  # use PCA from reference
  dims = 1:50,
  k.filter = NA,                # avoids over-filtering when reference is large
  verbose = TRUE)
DefaultAssay(reference) <- "SCT"
if (!"pca" %in% Reductions(reference)) {
  reference <- RunPCA(reference, verbose = FALSE)}
# Recompute UMAP on the reference **with** a stored model
reference <- RunUMAP(
  reference,
  reduction = "pca",
  dims = 1:50,
  return.model = TRUE,   # <-- this is the key bit
  verbose = FALSE)
multiome <- MapQuery(
  anchorset = anchors,
  query = multiome,
  reference = reference,
  refdata = list(seurat_label = "cell_type_leiden0.5"),
  reference.reduction = "pca",
  reduction.model = "umap")
# Inspect new metadata columns
grep("predicted|prediction", colnames(multiome@meta.data), value = TRUE)
# If you see "predicted.seurat_label", alias it to "seurat_label"
if ("predicted.seurat_label" %in% colnames(multiome@meta.data)) {
  multiome$seurat_label <- multiome$predicted.seurat_label
} else if ("predicted.id" %in% colnames(multiome@meta.data)) {
  # Some Seurat versions use 'predicted.id' by default
  multiome$seurat_label <- multiome$predicted.id}
# (Optional) also alias the max score for convenience
if ("prediction.score.seurat_label" %in% colnames(multiome@meta.data)) {
  multiome$prediction_score <- multiome$prediction.score.seurat_label
} else if ("prediction.score.max" %in% colnames(multiome@meta.data)) {
  multiome$prediction_score <- multiome$prediction.score.max}
# Make sure we really have the projected reduction
if (!"ref.umap" %in% Reductions(multiome)) {
  stop("No 'ref.umap' reduction found. Re-run MapQuery with reduction.model='umap', or plot another reduction (e.g., 'umap' or 'wnn.umap').")}
DimPlot(
  multiome,
  reduction = "ref.umap",
  group.by  = "seurat_label",
  label     = TRUE, repel = TRUE) + ggtitle("UMAP")
# saveRDS(multiome, file ="/scratch/pigblast/nb443/ValCARD_seqdata.1/seurat_multiome.rds")
# multiome <- readRDS("/scratch/pigblast/nb443/ValCARD_seqdata.1/seurat_multiome.rds")

meta <- multiome[[]]
# Label column (keep what you already had)
pred_col <- dplyr::case_when(
  "seurat_label" %in% names(meta) ~ "seurat_label",
  "predicted.seurat_label" %in% names(meta) ~ "predicted.seurat_label",
  "predicted.id" %in% names(meta) ~ "predicted.id",
  TRUE ~ NA_character_)
stopifnot(!is.na(pred_col))

# Score column — include your observed name
score_col <- dplyr::case_when(
  "prediction_score" %in% names(meta) ~ "prediction_score",
  "prediction.score.seurat_label" %in% names(meta) ~ "prediction.score.seurat_label",
  "prediction.score.max" %in% names(meta) ~ "prediction.score.max",
  "predicted.seurat_label.score" %in% names(meta) ~ "predicted.seurat_label.score",
  TRUE ~ NA_character_)

if (is.na(score_col)) {
  stop("No prediction score column found. Available columns:\n",
       paste(names(meta), collapse = ", "))}

# Optional: alias to a consistent name for downstream code
multiome$prediction_score <- as.numeric(multiome[[score_col]][,1])

# Build the df and plot as before
df <- meta |>
  dplyr::transmute(
    predicted_label = .data[[pred_col]],
    score = as.numeric(.data[[score_col]])
  ) |>
  dplyr::filter(!is.na(predicted_label), !is.na(score))

sum_df <- df |>
  dplyr::group_by(predicted_label) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_score = mean(score),
    se = sd(score) / sqrt(n),
    .groups = "drop"
  ) |>
  dplyr::arrange(mean_score)
p <- ggplot(sum_df, aes(x = reorder(predicted_label, mean_score), y = mean_score)) +
  geom_col(fill = "grey40") +
  geom_errorbar(aes(ymin = pmax(mean_score - se, 0), ymax = pmin(mean_score + se, 1)),
                width = 0.3) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Mean Seurat Prediction Score per Predicted Cell Type",
    x = "Predicted Cell Type",
    y = "Mean Prediction Score (± SE)"
  ) +
  theme_minimal(base_size = 13)
print(p)
# choose label column (fallback to orig.ident)
lab_col <- if ("seurat_label" %in% colnames(multiome@meta.data)) "seurat_label" else "orig.ident"
# -------- RNA-only UMAP --------
rna_assay <- if ("SCT" %in% Assays(multiome)) "SCT" else "RNA"
DefaultAssay(multiome) <- rna_assay
if (length(VariableFeatures(multiome)) == 0) {
  if (rna_assay == "SCT") {
    DefaultAssay(multiome) <- "RNA"
    multiome <- SCTransform(multiome, verbose = FALSE)
    DefaultAssay(multiome) <- "SCT"
  } else {
    multiome <- NormalizeData(multiome, verbose = FALSE)
    multiome <- FindVariableFeatures(multiome, nfeatures = 3000, verbose = FALSE)
    multiome <- ScaleData(multiome, features = VariableFeatures(multiome), verbose = FALSE)
  }
} else if (rna_assay == "RNA") {
  multiome <- ScaleData(multiome, features = VariableFeatures(multiome), verbose = FALSE)}
if (!"pca" %in% Reductions(multiome)) {
  multiome <- RunPCA(multiome, npcs = 50, verbose = FALSE)}
if (!"umap.rna" %in% Reductions(multiome)) {
  multiome <- RunUMAP(
    multiome, reduction = "pca", dims = 1:50,
    reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = FALSE)}
# 1) RNA plain (one colour)
multiome$.__one__ <- "cells"
multiome$.__one__ <- "cells"
p_rna_plain <- DimPlot(
  multiome,
  reduction = "umap.rna",
  group.by = ".__one__",
  cols = "#EC7063") +
  ggtitle("snRNA-seq") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_rna_plain)

# 2) RNA coloured by labels (no text labels on points)
p_rna_lab <- DimPlot(
  multiome,
  reduction = "umap.rna",
  group.by = lab_col,
  label = FALSE) +
  ggtitle("snRNA-seq with labels") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_rna_lab)
# -------- ATAC-only UMAP --------
atac_assay <- if ("ATAC" %in% Assays(multiome)) "ATAC" else if ("peaks" %in% Assays(multiome)) "peaks" else stop("No ATAC assay found.")
DefaultAssay(multiome) <- atac_assay
if (!"lsi" %in% Reductions(multiome)) {
  multiome <- RunTFIDF(multiome, verbose = FALSE)
  multiome <- FindTopFeatures(multiome, min.cutoff = "q0", verbose = FALSE)
  multiome <- RunSVD(multiome, verbose = FALSE)}
if (!"umap.atac" %in% Reductions(multiome)) {
  multiome <- RunUMAP(
    multiome, reduction = "lsi", dims = 2:30,
    reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE)}
# 3) ATAC plain (one colour)
p_atac_plain <- DimPlot(
  multiome,
  reduction = "umap.atac",
  group.by = ".__one__",
  cols = "#85C1E9") +
  ggtitle("snATAC-seq") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_atac_plain)
# 4) ATAC coloured by labels
p_atac_lab <- DimPlot(
  multiome,
  reduction = "umap.atac",
  group.by = lab_col,
  label = FALSE) +
  ggtitle("snATAC-seq with labels") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_atac_lab)

# -------- WNN co-embedding + overlay --------
DefaultAssay(multiome) <- rna_assay
stopifnot("pca" %in% Reductions(multiome), "lsi" %in% Reductions(multiome))
multiome <- FindMultiModalNeighbors(
  multiome, reduction.list = list("pca","lsi"),
  dims.list = list(1:50, 2:30))
if (!"umap.wnn" %in% Reductions(multiome)) {
  multiome <- RunUMAP(
    multiome, nn.name = "weighted.nn",
    reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", verbose = FALSE)}
# 5) “Combined RNA & ATAC” overlay (schematic overlay on same axes)
emb <- Embeddings(multiome, "umap.wnn") |> as.data.frame()
overlay_df <- rbind(
  transform(emb, modality = "RNA"),
  transform(emb, modality = "ATAC"))
p_overlay <- ggplot(overlay_df, aes(x = wnnUMAP_1, y = wnnUMAP_2, colour = modality)) +
  geom_point(size = 0.2, alpha = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "Combined RNA and ATAC UMAP", x = "UMAP_1", y = "UMAP_2")
print(p_overlay)
# 6) Co-embedded WNN UMAP coloured by labels
p_wnn_lab <- DimPlot(multiome, reduction = "umap.wnn",
                     group.by = lab_col, label = FALSE) +
  ggtitle("Co-embedded RNA + ATAC") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_wnn_lab)
# Save the co-embedded object
# saveRDS(multiome, file = "/scratch/pigblast/nb443/ValCARD_seqdata.1/multiome_seurat.rds")
# Check what metadata columns you have
head(multiome@meta.data)
# Count cells by sample and cell type
cell_counts <- multiome@meta.data %>%
  count(sample_id, seurat_label, name = "n_cells") %>%
  arrange(sample_id, seurat_label)
print(cell_counts)
# Replace "seurat_label" with the column in your metadata that holds cell type annotations
cell_type_counts <- multiome@meta.data %>%
  count(seurat_label, name = "n_cells") %>%
  arrange(desc(n_cells))
print(cell_type_counts)
