# -[]百分比
# -[]保存DEG結果
# -[]儲存交集、差集結果

# -[]根據樣本種類對cell type做DEG

# ── 套件 ─────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(VennDiagram)
library(gridExtra)

# ── 1. 找顯著 marker (每個 Seurat 各跑一次) ────────────
get_markers <- function(seu, id_col,
                        logfc = 0.25, p_adj = 0.05, min_pct = 0.1){
  Idents(seu) <- id_col
  markers <- FindAllMarkers(seu,
                            only.pos     = TRUE,
                            min.pct      = min_pct,
                            logfc.threshold = logfc,
                            test.use     = "wilcox") %>%    # 可改 t, LR, etc.
    filter(p_val_adj < p_adj) %>%                          # 依需求調整
    select(gene, cluster)
  return(markers)
}

markers1 <- get_markers(seurat_obj1, "Cell_Type_Compare_Merged")
markers2 <- get_markers(seurat_obj2, "Cell_Type_Compare_Merged")

# ── 2. 組 marker gene list ────────────────────────────
gene_list_from_markers <- function(markers_df){
  split(markers_df$gene, markers_df$cluster)
}

gene_sets1 <- gene_list_from_markers(markers1)
gene_sets2 <- gene_list_from_markers(markers2)

# 只保留兩個物件都有的 cell type
common_ct <- intersect(names(gene_sets1), names(gene_sets2))
gene_sets1 <- gene_sets1[common_ct]
gene_sets2 <- gene_sets2[common_ct]

if(length(common_ct) == 0){
  stop("兩個 Seurat 物件沒有任何共同的 Cell_Type_Compare_Merged。")
}

# ── 3. 畫 Venn & 收集到 list ──────────────────────────
venn_plots <- lapply(common_ct, function(ct){
  genes_a <- unique(gene_sets1[[ct]])
  genes_b <- unique(gene_sets2[[ct]])
  
  venn.plot <- venn.diagram(
    x = list(`Obj1` = genes_a, `Obj2` = genes_b),
    filename = NULL,
    category.names = c("Obj1", "Obj2"),
    main = ct,
    fill = c("#1f78b4", "#33a02c"),
    alpha = 0.6,
    cex   = 1.2,
    cat.cex = 1,
    margin = 0.05
  )
  grid::grid.grabExpr(grid::grid.draw(venn.plot))
})

# ── 4A. 直接顯示 (RStudio/Plots) ───────────────────────
# gridExtra::grid.arrange(grobs = venn_plots, ncol = 2)

# ── 套件 ─────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(VennDiagram)
library(gridExtra)

# ── 1. 找顯著 marker (每個 Seurat 各跑一次) ────────────
get_markers <- function(seu, id_col,
                        logfc = 0.25, p_adj = 0.05, min_pct = 0.1){
  Idents(seu) <- id_col
  markers <- FindAllMarkers(seu,
                            only.pos     = TRUE,
                            min.pct      = min_pct,
                            logfc.threshold = logfc,
                            test.use     = "wilcox") %>%    # 可改 t, LR, etc.
    filter(p_val_adj < p_adj) %>%                          # 依需求調整
    select(gene, cluster)
  return(markers)
}

markers1 <- get_markers(seurat_obj1, "Cell_Type_Compare_Merged")
markers2 <- get_markers(seurat_obj2, "Cell_Type_Compare_Merged")

# ── 2. 組 marker gene list ────────────────────────────
gene_list_from_markers <- function(markers_df){
  split(markers_df$gene, markers_df$cluster)
}

gene_sets1 <- gene_list_from_markers(markers1)
gene_sets2 <- gene_list_from_markers(markers2)

# 只保留兩個物件都有的 cell type
common_ct <- intersect(names(gene_sets1), names(gene_sets2))
gene_sets1 <- gene_sets1[common_ct]
gene_sets2 <- gene_sets2[common_ct]

if(length(common_ct) == 0){
  stop("兩個 Seurat 物件沒有任何共同的 Cell_Type_Compare_Merged。")
}

# ── 3. 畫 Venn & 收集到 list ──────────────────────────
venn_plots <- lapply(common_ct, function(ct){
  genes_a <- unique(gene_sets1[[ct]])
  genes_b <- unique(gene_sets2[[ct]])
  
  venn.plot <- venn.diagram(
    x = list(`Obj1` = genes_a, `Obj2` = genes_b),
    filename = NULL,
    category.names = c("Obj1", "Obj2"),
    main = ct,
    fill = c("#1f78b4", "#33a02c"),
    alpha = 0.6,
    cex   = 1.2,
    cat.cex = 1,
    margin = 0.05
  )
  grid::grid.grabExpr(grid::grid.draw(venn.plot))
})

# ── 4A. 直接顯示 (RStudio/Plots) ───────────────────────
gridExtra::grid.arrange(grobs = venn_plots, ncol = 5)

# ── 4B. 或存成 PDF 檔 ─────────────────────────────────
pdf("Venn_marker_overlap.pdf", width = 17, height = 10)
gridExtra::grid.arrange(grobs = venn_plots, ncol = 5)
dev.off()


