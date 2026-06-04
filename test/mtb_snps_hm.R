library(aricode)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(tibble)
library(viridisLite)

filt = 0
out_path <- ifelse(filt, "/hlilab/gianfranco/pgtools/test/mtb_heatmap_filt.pdf", "/hlilab/gianfranco/pgtools/test/mtb_heatmap.pdf")

mtb_md_path <- "/hlilab/gianfranco/pgtools/data/mtb_asm/mtb_md.tsv"
mtb_md <- read.csv(mtb_md_path, sep="\t")

mtb_snps_path <- ifelse(filt, "/hlilab/gianfranco/pgtools/results/mtb_asm_spec_filt.tsv", "/hlilab/gianfranco/pgtools/results/mtb_asm_spec.tsv")
mtb_snps <- read.csv(mtb_snps_path, sep="\t", check.names = FALSE)
colnames(mtb_snps)[4:ncol(mtb_snps)] <- sub("^.*asm/(.*)\\.fa\\.gz$", "\\1", colnames(mtb_snps)[4:ncol(mtb_snps)])

# Remove reference
mtb_snps <- mtb_snps %>% dplyr::select(-H37Rv)

# result <- mtb_snps %>%
#   group_by(across(1:2)) %>%
#   filter(
#     sum(pangenome, na.rm = TRUE) < 170,
#     all(pangenome >= 5)
#   ) %>%
#   ungroup()

mtb_mat <- as.matrix(mtb_snps[, 4:ncol(mtb_snps)])
storage.mode(mtb_mat) <- "numeric"

# Transform into binary matrix
mtb_mat_bin <- ifelse(mtb_mat != 0, 1, 0)

# Compute pairwise distances between samples (JI)
jd <- dist(t(mtb_mat_bin), method = "binary")

# Hierarchical clustering
hc <- hclust(jd, method = "complete")
mtb_clusters <- data.frame(ClusterID = cutree(hc, h = 0.1)) %>%
  rownames_to_column("SampleID")

# Set same order for metadata and lineages
mtb_md <- mtb_md[match(mtb_clusters$SampleID, mtb_md$SampleID), ]
NMI(mtb_clusters$ClusterID, mtb_md$PrimaryLineage)

value_colors <- colorRamp2(
  seq(min(1 - as.matrix(jd)), 1, length.out = 256),
  rocket(256, direction = -1))
cluster_colors <- c(
  "lineage1" = "#E63946",
  "lineage2" = "#2A9D8F",
  "lineage3" = "#E9C46A",
  "lineage4" = "#457B9D",
  "lineage5" = "#F4A261",
  "lineage6" = "#6A4C93",
  "lineage7" = "#2DC653",
  "lineage8" = "#FF6B6B"
)

pdf(out_path, width = 10, height = 8)

left_ha <- rowAnnotation(
  # Lineage = factor(mtb_clusters$ClusterID),
  Lineage = factor(mtb_md$PrimaryLineage),
  col = list(Lineage = cluster_colors),
  annotation_name_gp = gpar(fontsize = 0)
)

ht <- Heatmap(
  1 - as.matrix(jd),
  col = value_colors,
  left_annotation  = left_ha,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_dend    = FALSE,
  show_row_names   = FALSE,
  show_column_names = FALSE,
  name = "JS"
)

draw(ht)
dev.off()

