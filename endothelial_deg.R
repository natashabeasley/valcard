library(Seurat)
install.packages("Signac")
library(Signac)
library(Matrix)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
set.seed(123)

### --- For categorical analysis (B/C/D vs A) --- ###
# Load object
obj <- readRDS("/scratch/pigblast/nb443/ValCARD_seqdata.1/seurat_multiome.rds")
# Load the metadata CSV
md <- read.csv("/home/n/nb443/Documents/Val-CARD/dose_metadata.csv",
               stringsAsFactors = FALSE, check.names = FALSE)
# Make the column names R-safe (so "dose(mg/kg)" becomes "dose.mg.kg.")
colnames(md) <- make.names(colnames(md))
# Standardize sample_id on the object
obj$sample_id <- trimws(as.character(obj$sample_id))
# Map the metadata onto cells by sample_id (each cell inherits its sample’s dose info)
obj$dose_mgkg      <- md$dose.mg.kg.[ match(obj$sample_id, md$sample_id) ]
obj$dose_time_days <- md$dose_time.days.[ match(obj$sample_id, md$sample_id) ]
obj$sex <- md$sex[ match(obj$sample_id, md$sample_id) ]
# Ensure numeric
obj$dose_mgkg      <- as.numeric(obj$dose_mgkg)
obj$dose_time_days <- as.numeric(obj$dose_time_days)
print(aggregate(cbind(dose_mgkg = obj$dose_mgkg,
                      dose_time_days = obj$dose_time_days,
                      sex = obj$sex),
                by = list(sample_id = obj$sample_id), FUN = function(x) unique(na.omit(x))))
GroupA <- c("V014", "V032", "V035", "V041")
GroupB <- c("V003", "V026", "V034", "V031", "V007")
GroupC <- c("V010", "V017", "V028", "V045", "V004", "V043")
GroupD <- c("V011", "V027", "V030")
## --- DEG on Endothelial cells: Groups B/C/D vs A (NO dose covariate) --------
# keep only cardiomyocytes with group assigned (no dose filter)
DefaultAssay(obj) <- "RNA"
## map sample_id -> group on the object
grp_map <- c(
  setNames(rep("A", length(GroupA)), GroupA),
  setNames(rep("B", length(GroupB)), GroupB),
  setNames(rep("C", length(GroupC)), GroupC),
  setNames(rep("D", length(GroupD)), GroupD))
obj$group <- unname(grp_map[as.character(obj$sample_id)])
obj$group <- factor(obj$group, levels = c("A","B","C","D"))
# quick check 
print(table(obj$group, useNA = "ifany"))
cm <- subset(
  obj,
  subset = !is.na(group) &
    seurat_label %in% c(" Endothelial I", " Endothelial II"))
# normalize + join layers (Seurat v5)
cm <- NormalizeData(cm, verbose = FALSE)
cm <- JoinLayers(cm, assay = DefaultAssay(cm))
Idents(cm) <- "group"
# DEG helper (LR; no latent variables)
cm$sex <- factor(cm$sex)  
deg <- function(seu, g1, g0, covar = "sex") {
  # subset to just the groups being compared
  sub <- subset(seu, idents = c(g1, g0))
  # check covariate variation
  cov_vals <- sub[[covar]][, 1]
  has_var <- !(all(is.na(cov_vals)) || length(unique(cov_vals)) < 2)
  x <- FindMarkers(
    sub, ident.1 = g1, ident.2 = g0,
    test.use = "LR",
    latent.vars = if (has_var) covar else NULL,  
    logfc.threshold = 0.25, min.pct = 0.1,
    verbose = FALSE)
  x$FDR <- p.adjust(x$p_val, "BH")
  x$log2FC <- if ("avg_log2FC" %in% names(x)) x$avg_log2FC else x$avg_logFC
  x[order(x$FDR, -abs(x$log2FC)), , drop = FALSE]}
# balanced (same # cells per group) 
set.seed(123)
nc <- table(Idents(cm))
k  <- min(nc[names(nc) %in% c("A","B","C","D")])
pick <- unlist(lapply(c("A","B","C","D"), function(g) {
  sample(WhichCells(cm, idents = g), k)}))
cm_bal <- subset(cm, cells = pick)
deg_BA_bal <- deg(cm_bal, "B", "A")
deg_CA_bal <- deg(cm_bal, "C", "A")
deg_DA_bal <- deg(cm_bal, "D", "A")
table(Idents(cm_bal))   # should all equal eachother
# ---- save DEG results (balanced) -----------------------------------------
save_deg <- function(df, label,
                     outdir = "/rfs/CardiacSurgeryRes/nb443/VALCARD/endothelial/deg_results_balanced") {
  df_all <- tibble::rownames_to_column(as.data.frame(df), var = "gene")
  # infer group names from the label like "B_vs_A"
  g1 <- sub("_vs_.*", "", label)
  g0 <- sub(".*_vs_", "", label)
  # add direction + convenience columns
  df_all <- df_all %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        is.na(log2FC)           ~ NA_character_,
        log2FC > 0              ~ paste0("Up_in_", g1),
        log2FC < 0              ~ paste0("Down_in_", g1),
        TRUE                    ~ "No_change"),
      contrast = paste0(g1, " vs ", g0),
      sig = FDR < 0.05 & abs(log2FC) >= 0.25)
  # write files
  write.csv(df_all, file.path(outdir, paste0(label, ".all.csv")), row.names = FALSE)
  df_sig <- dplyr::filter(df_all, sig)
  write.csv(df_sig, file.path(outdir, paste0(label, ".sig.csv")), row.names = FALSE)}
save_deg(deg_BA_bal, "B_vs_A")
save_deg(deg_CA_bal, "C_vs_A")
save_deg(deg_DA_bal, "D_vs_A")

# significance rule
filt <- function(d) subset(d, FDR < 0.05 & abs(log2FC) >= 0.25)
# Counts bar plot (Up/Down)
cnt <- function(d) {
  s <- filt(d)
  data.frame(up = sum(s$log2FC > 0, na.rm=TRUE),
             down = sum(s$log2FC < 0, na.rm=TRUE),
             total = nrow(s))}
summ_bal <- rbind(
  `B vs A` = cnt(deg_BA_bal),
  `C vs A` = cnt(deg_CA_bal),
  `D vs A` = cnt(deg_DA_bal))
df_bal <- rbind(
  data.frame(contrast = rownames(summ_bal), direction = "Upregulated",   n = summ_bal[, "up"]),
  data.frame(contrast = rownames(summ_bal), direction = "Downregulated", n = summ_bal[, "down"]))
ggplot(df_bal, aes(contrast, n, fill = direction)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Number of significant DEGs (FDR < 0.05, |log2FC| ≥ 0.25)",
       title = "Significant DEGs - Endothelial cells") +
  scale_fill_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue")) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))
# Faceted volcano plots 
as_volc <- function(d) {
  d$neglog10FDR <- -log10(pmax(d$FDR, .Machine$double.xmin))
  d$dir <- ifelse(d$FDR < 0.05 & d$log2FC >  0.25, "Up",
                  ifelse(d$FDR < 0.05 & d$log2FC < -0.25, "Down", "NS"))
  d}
# rebuild v_bal keeping gene names
v_bal <- dplyr::bind_rows(
  as_volc(deg_BA_bal) %>% rownames_to_column("gene") %>% mutate(contrast = "B vs A"),
  as_volc(deg_CA_bal) %>% rownames_to_column("gene") %>% mutate(contrast = "C vs A"),
  as_volc(deg_DA_bal) %>% rownames_to_column("gene") %>% mutate(contrast = "D vs A"))
# pick top 10 Up and top 10 Down per contrast (by significance, then effect size)
top_lab <- v_bal %>%
  filter(dir %in% c("Up","Down")) %>%
  group_by(contrast, dir) %>%
  arrange(desc(neglog10FDR), desc(abs(log2FC))) %>%
  slice_head(n = 10) %>%
  ungroup()
# volcano with labels
ggplot(v_bal, aes(log2FC, neglog10FDR, color = dir)) +
  geom_point(alpha = 0.7, size = 0.9) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_wrap(~ contrast, nrow = 1) +
  geom_text_repel(
    data = top_lab,
    aes(label = gene),
    size = 2.7, box.padding = 0.3, max.overlaps = 100,
    segment.color = "grey50", min.segment.length = 0
  ) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  labs(x = "log2 fold-change", y = expression(-log[10](FDR)),
       title = "Volcano plots", color = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# pick genes to plot (example: top sig DEGs across contrasts)
pick_top <- function(d, n = 30) {
  d <- as.data.frame(d)
  d <- d[d$FDR < 0.05 & abs(d$log2FC) >= 0.25, , drop = FALSE]
  rownames(d[order(d$FDR, -abs(d$log2FC)), , drop = FALSE])[seq_len(min(n, nrow(d)))]}
genes <- unique(c(
  pick_top(deg_BA_bal, 30),
  pick_top(deg_CA_bal, 30),
  pick_top(deg_DA_bal, 30)))
genes <- genes[genes %in% rownames(cm_bal)]  # keep genes present
# scale those genes and extract the matrix
cm_bal <- ScaleData(cm_bal, features = genes, verbose = FALSE)
mat <- GetAssayData(cm_bal, slot = "scale.data")[genes, , drop = FALSE]
# order cells by group and make annotations
cells_ordered <- names(sort(Idents(cm_bal)))  # order by A, B, C, D (factor levels)
mat_ord <- mat[, cells_ordered, drop = FALSE]
group_ord <- Idents(cm_bal)[cells_ordered]
group_cols <- c(A="#4DAF4A", B="#377EB8", C="#E41A1C", D="#984EA3")
# draw heatmap
ComplexHeatmap::Heatmap(
  mat_ord, name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(Group = group_ord,
                                                     col = list(Group = group_cols)),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_split = group_ord,
  column_title = NULL,
  row_title = "Top DEGs",
  heatmap_legend_param = list(title = "Expression"))
VlnPlot(cm_bal, features = c("TFPI"), group.by = "group",
        pt.size = 0.1, cols = c("grey","red","blue","green"))

## GO analysis
org_db <- org.Hs.eg.db
# your significance rule for DEGs
filt <- function(d) subset(as.data.frame(d), FDR < 0.05 & abs(log2FC) >= 0.25)
# background/universe = all genes in the endothelial subset used for DEG
DefaultAssay(cm_bal) <- "RNA"
universe_genes <- rownames(cm_bal)
# SYMBOL -> ENTREZ
to_entrez <- function(symbols) {
  suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  ) %>% distinct(SYMBOL, .keep_all = TRUE)}
univ_entrez <- to_entrez(universe_genes)$ENTREZID
# GO:BP only
do_go_bp <- function(entrez_vec, universe_entrez = univ_entrez,
                     pcut = 0.05, qcut = 0.05) {
  enrichGO(
    gene          = entrez_vec,
    universe      = universe_entrez,
    OrgDb         = org_db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = pcut,
    qvalueCutoff  = qcut,
    readable      = TRUE)}
# Show dotplot for a GO result if it has rows
show_bp_dotplot <- function(ego, title, top = 20) {
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    print(clusterProfiler::dotplot(ego, showCategory = top) + ggtitle(title))}}
# For one contrast: run GO for Up and Down, show plots, save only significant terms to Excel
go_bp_view_and_save <- function(deg, label,
                                outdir = "/rfs/CardiacSurgeryRes/nb443/VALCARD/endothelial/go_enrichment",
                                top_show = 20) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  sig_tbl <- filt(deg)
  up_genes   <- rownames(sig_tbl[sig_tbl$log2FC > 0, , drop = FALSE])
  down_genes <- rownames(sig_tbl[sig_tbl$log2FC < 0, , drop = FALSE])
  up_entrez   <- to_entrez(up_genes)$ENTREZID
  down_entrez <- to_entrez(down_genes)$ENTREZID
  ego_up   <- if (length(up_entrez))   do_go_bp(up_entrez)   else NULL
  ego_down <- if (length(down_entrez)) do_go_bp(down_entrez) else NULL
  # View plots (Up and Down separately)
  show_bp_dotplot(ego_up,   paste0(label, " — GO:BP (Up)"),   top = top_show)
  show_bp_dotplot(ego_down, paste0(label, " — GO:BP (Down)"), top = top_show)
  #Save only significant ones to Excel (.xlsx)
  # clusterProfiler already respects p/q cutoffs; still filter by p.adjust for safety
  df_up   <- if (!is.null(ego_up))   as.data.frame(ego_up)   else NULL
  df_down <- if (!is.null(ego_down)) as.data.frame(ego_down) else NULL
  df_up   <- if (!is.null(df_up)   && nrow(df_up))   dplyr::filter(df_up, p.adjust <= 0.05)   else NULL
  df_down <- if (!is.null(df_down) && nrow(df_down)) dplyr::filter(df_down, p.adjust <= 0.05) else NULL
  # If at least one has rows, write an Excel with two sheets (missing sheets omitted)
  sheets <- list()
  if (!is.null(df_up)   && nrow(df_up))   sheets$Up_BP   <- df_up
  if (!is.null(df_down) && nrow(df_down)) sheets$Down_BP <- df_down
  if (length(sheets) > 0) {
    writexl::write_xlsx(sheets, file.path(outdir, paste0(label, "_GO_BP.xlsx")))} else {
    message(label, ": no significant GO:BP terms to save at p.adjust <= 0.05.")}
  invisible(list(up = ego_up, down = ego_down))}
# Run per condition (separate) 
go_bp_view_and_save(deg_BA_bal, "B_vs_A")
go_bp_view_and_save(deg_CA_bal, "C_vs_A")
go_bp_view_and_save(deg_DA_bal, "D_vs_A")