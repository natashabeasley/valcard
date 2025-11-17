# Peoducing QC figures for val card RNA+ATAC 
set.seed(42)
library(Seurat)
library(Signac)
library(future)
options(future.globals.maxSize = +Inf)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(DoubletFinder)
library(Matrix)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'  # make sure chromosome names match your data
genome(annotations) <- "hg38"

samples_df <- read.csv("/home/n/nb443/Documents/Val-CARD/file_locations.csv", header = TRUE)
# Rename columns to standard names
colnames(samples_df) <- c("sample_id", "path")
# Remove extra whitespace in paths
samples_df$path <- gsub(" ", "", samples_df$path)
# Remove any empty rows (optional safety check)
samples_df <- na.omit(samples_df)
samples_df <- samples_df[samples_df$sample_id != "" & samples_df$path != "", ]
tss.positions <- promoters(genes(EnsDb.Hsapiens.v86), upstream = 0, downstream = 1)
tss.positions <- keepStandardChromosomes(tss.positions, pruning.mode = "coarse")
seqlevelsStyle(tss.positions) <- "UCSC"

## ---------- helper: is this an RNA 10x dir? ----------
is_rna_path <- function(p) {
  dir.exists(file.path(p, "filtered_feature_bc_matrix")) ||
    length(list.files(p, pattern = "barcodes.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)) > 0}

## ---------- load a single RNA sample ----------
load_rna <- function(path, sample_id) {
  mdir <- if (dir.exists(file.path(path, "filtered_feature_bc_matrix"))) {
    file.path(path, "filtered_feature_bc_matrix")} else path
  
  counts <- Read10X(mdir)
  if (is.list(counts)) {
    counts <- counts[[if ("Gene Expression" %in% names(counts)) "Gene Expression" else 1]]}
  so <- CreateSeuratObject(counts, project = sample_id, min.cells = 3, min.features = 200)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so$Sample <- sample_id
  so}

## ---------- load all RNA objects (skip non-RNA paths) ----------
rna_rows <- which(vapply(samples_df$path, is_rna_path, logical(1)))
rna_objs <- lapply(rna_rows, function(i) {
  sid <- as.character(samples_df$sample_id[i])
  pth <- as.character(samples_df$path[i])
  message("Loading RNA: ", sid)
  tryCatch(load_rna(pth, sid), error = function(e) { message("  -> skipped: ", e$message); NULL })})
names(rna_objs) <- samples_df$sample_id[rna_rows]
rna_objs <- Filter(Negate(is.null), rna_objs)

stopifnot(length(rna_objs) > 0)

## ---------- gather metadata (no Seurat merge) ----------
meta_list <- lapply(names(rna_objs), function(sid) {
  md <- rna_objs[[sid]]@meta.data
  md$Sample <- sid
  md})
meta_all <- dplyr::bind_rows(meta_list)

## ---------- shared style ----------
theme_set(theme_bw(base_size = 12))
rotate_x <- theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
pal <- hue_pal()(length(unique(meta_all$Sample)))
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"),
             plot.title.position = "plot")  # optional, nicer alignment

## ============ A: Cell Number per Sample ============
cell_counts <- meta_all %>%
  count(Sample, name = "Cells") %>%
  arrange(Sample) %>%
  mutate(Sample = factor(Sample, levels = Sample))

pA <- ggplot(cell_counts, aes(Sample, Cells, fill = Sample)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = pal) +
  labs(title = "Cell Number per Sample", x = "Sample", y = "Number of Cells") +
  rotate_x + theme(legend.position = "none")
print(pA)

## ============ B: Mitochondrial Percentage per Sample ============
meta_all$Sample <- factor(meta_all$Sample,
                          levels = levels(cell_counts$Sample))  # keep same order
pB <- ggplot(meta_all, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_jitter(size = 0.15, alpha = 0.2, width = 0.2, height = 0) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Mitochondrial Percentage per Sample",
       x = "Sample", y = "Mitochondrial Content (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(pB)

## ============ C: RNA Counts per Sample ============
pC <- ggplot(meta_all, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(scale = "width", trim = FALSE) +                 # no boxplot fill
  geom_jitter(size = 0.15, alpha = 0.2, width = 0.2, height = 0) +  # lighter dots
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "RNA Counts per Sample", x = "Sample", y = "Total RNA counts") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(pC)

## ============ D: Number of Features per Sample ============
pD <- ggplot(meta_all, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(scale = "width", trim = FALSE) +                 # no white box in middle
  geom_jitter(size = 0.15, alpha = 0.2, width = 0.2, height = 0) +  # lighter dots
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Number of Features per Sample",
       x = "Sample", y = "Number of Detected Genes") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(pD)

## ============ E: Per-sample averages ============
rna_summary <- meta_all %>%
  group_by(Sample) %>%
  summarise(
    `Cell number` = dplyr::n(),
    `Average mitochondrial percentage (%)` = mean(percent.mt, na.rm = TRUE),
    `Average RNA count`     = mean(nCount_RNA,   na.rm = TRUE),
    `Average RNA features`  = mean(nFeature_RNA, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(Sample)

print(rna_summary)     # view in console
# View(rna_summary)    # uncomment in RStudio to open as a table

load_multiome <- function(sample_id, sample_path) {
  cat("Processing:", sample_id, "\n")
  # Read 10X data from filtered feature matrix folder
  data <- Read10X(data.dir = sample_path)
  # Create RNA assay
  seurat_obj <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA", project = sample_id)
  # Load gene annotations
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  # Create ATAC assay
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(
    counts = data$Peaks,
    sep = c(":", "-"),
    genome = "hg38",
    annotations = annotations)
  # Correct fragments path (pointing to outs folder)
  fragments_path <- file.path(dirname(sample_path), "atac_fragments.tsv.gz")
  # Check if fragment file exists
  if (!file.exists(fragments_path)) {
    stop(paste("Fragment file does not exist:", fragments_path))
  } else {
    cat("Using fragments file:", fragments_path, "\n")}
  # Create Fragment object for ATAC assay
  cells <- colnames(seurat_obj)
  seurat_obj[["ATAC"]]@fragments <- list(CreateFragmentObject(path = fragments_path, cells = cells))
  # Rename cells to include sample id
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_id)
  seurat_obj$sample_id <- sample_id
  # Percent mitochondrial RNA
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "RNA")
  # Calculate TSS positions if not global
  if (!exists("tss.positions")) {
    tss.positions <- promoters(genes(EnsDb.Hsapiens.v86), upstream = 0, downstream = 1)
    tss.positions <- keepStandardChromosomes(tss.positions, pruning.mode = "coarse")
    seqlevelsStyle(tss.positions) <- "UCSC"}
  # Calculate TSS enrichment
  seurat_obj <- TSSEnrichment(
    object = seurat_obj,
    tss.positions = tss.positions,
    fast = TRUE,
    assay = "ATAC",
    region_extension = 1000,
    verbose = TRUE)
  # Calculate nucleosome signal
  seurat_obj <- NucleosomeSignal(
    object = seurat_obj,
    assay = "ATAC",
    verbose = TRUE)
  # Set default assay back to RNA for filtering
  DefaultAssay(seurat_obj) <- "RNA"
  # Filter cells with QC thresholds - adjust thresholds as needed
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > 200 &
                         nFeature_RNA < 12500 &
                         percent.mt < 5 &
                         nCount_ATAC > 13 &
                         nCount_ATAC < 100000 &
                         nucleosome_signal < 1 &
                         TSS.enrichment > 1)
  return(seurat_obj)}
multiome_list <- lapply(1:nrow(samples_df), function(i) {
  load_multiome(samples_df$sample_id[i], samples_df$path[i])}) # Takes approx 2 hours
names(multiome_list) <- samples_df$sample_id
# Rename list elements by sample ID
names(multiome_list) <- samples_df$sample_id

# Save the list of Seurat objects to an .rds file - rstudio sometimes crashes at this point
# saveRDS(multiome_list, file = "/scratch/pigblast/nb443/ValCARD_seqdata.1/valcard_list_pre_doublet_removal.rds")

# Normalize each sample with SCTransform
multiome_list_sct <- lapply(X = multiome_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, verbose = TRUE, variable.features.n = 3000)})

# Check cell counts for each sample - total cells = 22,341
sapply(multiome_list_sct, ncol)

# Run PCA - default = 50 cells however some samples have less than 50 cells, reduce to 30 -> change to 100?.
multiome_list_pca <- lapply(X = multiome_list_sct, FUN = function(x) {
  # Determine safe number of PCs: can't be more than min(features, cells)
  max_pcs <- min(ncol(x), nrow(VariableFeatures(x))) - 1
  safe_pcs <- min(100, max_pcs)  # Use up to 30 PCs max, or whatever fits
  RunPCA(x, verbose = TRUE, npcs = safe_pcs)})

# Run UMAP - default n_neighbours = 30
multiome_list_umap <- lapply(X = multiome_list_pca, FUN = function(x) {
  n_cells <- ncol(x)
  neighbors <- min(30, n_cells - 1)  # Must be < number of cells
  RunUMAP(x, dims = 1:min(10, ncol(x[["pca"]])), n.neighbors = neighbors, verbose = TRUE)})

# Run DoubletFinder and remove doublets
multiome_list_rm_doublets <- lapply(X = multiome_list_umap, FUN = function(x, PCs = 1:10, pN = 0.25, doublet_rate = 0.008) {
  n_cells <- ncol(x)
  # Skip if sample has too few cells
  if (n_cells < 100) {
    message("Skipping doublet detection for sample with <100 cells: ", x@project.name)
    return(x)}
  # Step 1: Parameter sweep
  sweep.res.list <- tryCatch({
    paramSweep(x, PCs = PCs, sct = TRUE)
  }, error = function(e) {
    message("paramSweep error: ", e$message, " in ", x@project.name)
    return(NULL)})
  if (is.null(sweep.res.list)) return(x)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk <- as.numeric(as.character(bcmvn.max$pK))
  # Step 2: Estimate expected doublets
  annotations <- if (!is.null(x@meta.data$seurat_clusters)) x@meta.data$seurat_clusters else rep(1, ncol(x))
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate * n_cells)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  if (nExp_poi.adj < 2) {
    message("Too few expected doublets (", nExp_poi.adj, "), skipping: ", x@project.name)
    return(x)}
  # Step 3: Run doubletFinder
  x <- doubletFinder(x,
                     PCs = PCs,
                     pN = pN,
                     pK = optimal.pk,
                     nExp = nExp_poi.adj,
                     sct = TRUE)
  # Step 4: Subset singlets
  df_col <- grep("DF.classifications", colnames(x@meta.data), value = TRUE)
  if (length(df_col) > 0) {
    colnames(x@meta.data)[which(colnames(x@meta.data) == df_col)] <- "doublet_finder"
    x <- subset(x, subset = doublet_finder == "Singlet")}
  return(x)})

# Optionally assign names again
names(multiome_list_rm_doublets) <- names(multiome_list)

# Check cell numbers now it has been filtered and doublets removed
sapply(multiome_list_rm_doublets, ncol)
saveRDS(multiome_list_rm_doublets, file = "/scratch/pigblast/nb443/ValCARD_seqdata.1/valcard_list_rm_doublets_min_100")

# ---- 1) Collect ATAC metadata from each sample ----
stopifnot(exists("multiome_list"), length(multiome_list) > 0)
meta_atac <- dplyr::bind_rows(lapply(multiome_list, function(obj) {
  md <- obj@meta.data
  # sample label you stored earlier
  md$Sample <- if ("sample_id" %in% names(md)) md$sample_id else obj@project.name
  
  # FRiP if possible (else NA)
  if (!"FRiP" %in% names(md)) {
    if (all(c("nCount_ATAC","passed_filters") %in% names(md))) {
      md$FRiP <- md$nCount_ATAC / md$passed_filters
    } else if ("pct_reads_in_peaks" %in% names(md)) {
      md$FRiP <- md$pct_reads_in_peaks / 100
    } else {
      md$FRiP <- NA_real_}}md}))

# ---- 2) Plot styling (match RNA look) ----
theme_set(theme_bw(base_size = 12))
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))
rotate_x <- theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

# consistent sample order (alphabetical; change if you prefer by cell count)
sample_order <- meta_atac %>% count(Sample, name = "Cells") %>% arrange(Sample) %>% pull(Sample)
meta_atac$Sample <- factor(meta_atac$Sample, levels = sample_order)
pal <- hue_pal()(nlevels(meta_atac$Sample))
pt_alpha <- 0.2

violin_jitter <- function(df, y, title, ylab) {
  ggplot(df, aes(x = Sample, y = .data[[y]], fill = Sample)) +
    geom_violin(scale = "width", trim = FALSE, na.rm = TRUE) +
    geom_jitter(size = 0.15, alpha = pt_alpha, width = 0.2, height = 0, na.rm = TRUE) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = "Sample", y = ylab) +
    rotate_x + theme(legend.position = "none")}

# ---- 3) ATAC QC plots (visualize only) ----
# A) Number of ATAC features per cell
pA_atac_features <- violin_jitter(meta_atac, "nFeature_ATAC",
                                  "Number of ATAC Features per Cell", "ATAC features")
print(pA_atac_features)

# B) Fraction of Reads in Peaks (FRiP) — will be blank if FRiP not available
pB_atac_frip <- violin_jitter(meta_atac, "FRiP",
                              "Fraction of Reads in Peaks per Cell", "Fraction in Peaks")
print(pB_atac_frip)

# C) Nucleosome signal per cell
pC_atac_nucleosome <- violin_jitter(meta_atac, "nucleosome_signal",
                                    "Nucleosome Signal per Cell", "Nucleosome Signal")
print(pC_atac_nucleosome)

# D) Fragment count per cell
pD_atac_fragcount <- violin_jitter(meta_atac, "nCount_ATAC",
                                   "Fragment Count per Cell", "Fragments / Cell")
print(pD_atac_fragcount)

# E) TSS enrichment
pE_atac_tss <- violin_jitter(meta_atac, "TSS.enrichment",
                             "TSS Enrichment", "TSS Enrichment")
print(pE_atac_tss)

# ---- 4) Per-sample averages (printed) ----
atac_summary_table <- meta_atac %>%
  group_by(Sample) %>%
  summarise(
    `Average fragment count per cell` = mean(nCount_ATAC,   na.rm = TRUE),
    FRiP                               = mean(FRiP,          na.rm = TRUE),
    `Mean Nucleosome Signal per Cell`  = mean(nucleosome_signal, na.rm = TRUE),
    `Mean number of features per cell` = mean(nFeature_ATAC, na.rm = TRUE),
    `Average TSS enrichment`           = mean(TSS.enrichment, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(Sample)
print(atac_summary_table)
# View(atac_summary_table)  # uncomment in RStudio to open as a sheet
# saveRDS(multiome_list, file = "/scratch/pigblast/nb443/ValCARD_seqdata.1/merged_filtered.rds")

# Filtering RNA + ATAC seperatley
# ---- packages ----
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(dplyr)
  library(readr)})
# ---- input samples table (same as yours) ----
samples_df <- read.csv("/home/n/nb443/Documents/Val-CARD/file_locations.csv", header = TRUE)
colnames(samples_df) <- c("sample_id", "path")
samples_df$path <- gsub(" ", "", samples_df$path)
samples_df <- na.omit(samples_df)
samples_df <- samples_df[samples_df$sample_id != "" & samples_df$path != "", ]
# ---- prepare TSS positions once ----
tss.positions <- promoters(genes(EnsDb.Hsapiens.v86), upstream = 0, downstream = 1)
tss.positions <- keepStandardChromosomes(tss.positions, pruning.mode = "coarse")
seqlevelsStyle(tss.positions) <- "UCSC"
# ---- loader that DOES NOT subset; only computes QC + flags ----
load_multiome_with_qc <- function(sample_id, sample_path) {
  message("Processing: ", sample_id)
  # Read 10x multiome output (filtered_feature_bc_matrix)
  data <- Read10X(data.dir = sample_path)
  # Create base object (RNA)
  seu <- CreateSeuratObject(
    counts  = data$`Gene Expression`,
    assay   = "RNA",
    project = sample_id)
  # Annotations for ATAC
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  # Add ATAC assay
  seu[["ATAC"]] <- CreateChromatinAssay(
    counts      = data$Peaks,
    sep         = c(":", "-"),
    genome      = "hg38",
    annotations = annotations)
  # Fragments path (10x outs folder next to filtered_feature_bc_matrix)
  fragments_path <- file.path(dirname(sample_path), "atac_fragments.tsv.gz")
  if (!file.exists(fragments_path)) {
    stop(paste("Fragment file does not exist:", fragments_path))} else {
    message("Using fragments file: ", fragments_path)}
  # Attach fragments (Signac)
  cells <- colnames(seu)
  seu[["ATAC"]]@fragments <- list(CreateFragmentObject(path = fragments_path, cells = cells))
  # Rename cells to include sample id
  seu <- RenameCells(seu, add.cell.id = sample_id)
  seu$sample_id <- sample_id
  # RNA QC: %MT
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-", assay = "RNA")
  # ATAC QC: TSS, Nucleosome
  DefaultAssay(seu) <- "ATAC"
  seu <- TSSEnrichment(
    object           = seu,
    tss.positions    = tss.positions,
    fast             = TRUE,
    assay            = "ATAC",
    region_extension = 1000,
    verbose          = FALSE)
  seu <- NucleosomeSignal(object = seu, assay = "ATAC", verbose = FALSE)
  # Back to RNA as default
  DefaultAssay(seu) <- "RNA"
  # ---- add pass/fail flags (hresholds) ----
  seu$pass_rna  <- with(seu@meta.data,
                        nFeature_RNA > 200 &
                          nFeature_RNA < 12500 &
                          percent.mt   < 5)
  seu$pass_atac <- with(seu@meta.data,
                        nCount_ATAC        > 13 &
                          nCount_ATAC        < 100000 &
                          nucleosome_signal  < 1 &
                          TSS.enrichment     > 1)
  seu$pass_both <- seu$pass_rna & seu$pass_atac
  return(seu)}
# ---- run over all samples ----
obj_list <- lapply(seq_len(nrow(samples_df)), function(i) {
  load_multiome_with_qc(samples_df$sample_id[i], samples_df$path[i])})
names(obj_list) <- samples_df$sample_id
# ---- per-sample summaries ----
summary_df <- bind_rows(lapply(obj_list, function(seu) {
  md <- seu@meta.data
  tibble(
    sample_id      = unique(md$sample_id),
    cells_total    = nrow(md),
    cells_pass_rna = sum(md$pass_rna,  na.rm = TRUE),
    cells_pass_atac= sum(md$pass_atac, na.rm = TRUE),
    cells_pass_both= sum(md$pass_both, na.rm = TRUE))}))
# Useful convenience: which samples effectively “failed” one modality (no cells passing)
failed_rna_samples  <- summary_df %>% filter(cells_pass_rna == 0)  %>% pull(sample_id)
failed_atac_samples <- summary_df %>% filter(cells_pass_atac == 0) %>% pull(sample_id)
# ---- (optional) write summary to disk ----
write_csv(summary_df, "multiome_qc_summary_by_sample.csv")
# ---- (optional) get per-cell barcode lists if you want to subset later ----
cells_pass_rna_by_sample <- lapply(obj_list, function(seu) rownames(seu@meta.data)[seu$pass_rna])
cells_pass_atac_by_sample<- lapply(obj_list, function(seu) rownames(seu@meta.data)[seu$pass_atac])
cells_pass_both_by_sample<- lapply(obj_list, function(seu) rownames(seu@meta.data)[seu$pass_both])
# ---- (optional) create filtered objects now, if/when you’re ready ----
# Filtered by RNA only
obj_list_rna_only  <- lapply(obj_list, function(seu) subset(seu, cells = WhichCells(seu, expression = pass_rna)))
# Filtered by ATAC only
obj_list_atac_only <- lapply(obj_list, function(seu) subset(seu, cells = WhichCells(seu, expression = pass_atac)))
# Filtered by both
obj_list_both      <- lapply(obj_list, function(seu) subset(seu, cells = WhichCells(seu, expression = pass_both)))
library(dplyr)
library(stringr)
summary_ordered <- summary_view %>%
  mutate(sample_num = as.integer(str_extract(sample_id, "\\d+$"))) %>%
  arrange(sample_num)
print(summary_ordered, n = Inf)
