library(Seurat)
install.packages("Signac")
library(Signac)
library(Matrix)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
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
## --- DEG on cardiomyocytes: Groups B/C/D vs A (NO dose covariate) --------
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
    grepl("cardio", seurat_label, ignore.case = TRUE))
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
# ---- save DEG results (balanced) -----------------------------------------
dir.create("/rfs/CardiacSurgeryRes/nb443/VALCARD/cardiomyocytes/deg_results_balanced", showWarnings = FALSE)
# significance rule (reuse your existing one)
filt <- function(d) subset(d, FDR < 0.05 & abs(log2FC) >= 0.25)
save_deg <- function(df, label, outdir = "/rfs/CardiacSurgeryRes/nb443/VALCARD/cardiomyocytes/deg_results_balanced") {
  df_all <- tibble::rownames_to_column(as.data.frame(df), var = "gene")
  write.csv(df_all, file.path(outdir, paste0(label, ".all.csv")), row.names = FALSE)
  df_sig <- tibble::rownames_to_column(as.data.frame(filt(df)), var = "gene")
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
    title = "Significant DEGs") +
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

VlnPlot(cm_bal, features = c("S100A10"), group.by = "group",
                 pt.size = 0.1, cols = c("grey","red","blue","green"))
