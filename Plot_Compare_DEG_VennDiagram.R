# -[T]百分比
# -[]保存DEG結果
# -[]儲存交集、差集結果
# -[]TopN gene 模式
# -[]Object名稱設定

# -[]Fig legend

# -[]根據樣本種類對cell type做DEG



# ---------- 套件 ----------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(VennDiagram)
  library(gridExtra)
})

# ---------- 可調參數 ----------
Name_Seurat_obj1 <- "Keloid_Charlene"      # 物件 1 名稱
Name_Seurat_obj2 <- "Keloid_Jojie"         # 物件 2 名稱
Color_Seurat_obj1  <- "#1f78b4"            # 物件 1 顏色
Color_Seurat_obj2  <- "#33a02c"            # 物件 2 顏色
outdir    <- "Venn_output"    # 所有結果輸出的資料夾
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 找顯著 marker ----------
get_markers <- function(seu, id_col,
                        logfc = 0.25, p_adj = 0.05, min_pct = 0.10) {
  
  Idents(seu) <- id_col
  markers <- FindAllMarkers(
    seu,
    only.pos       = TRUE,
    min.pct        = min_pct,
    logfc.threshold= logfc,
    test.use       = "wilcox"
  ) %>% 
    filter(p_val_adj < p_adj) %>% 
    select(gene, cluster)
  
  return(markers)
}

markers1 <- get_markers(seurat_obj1, "Cell_Type_Compare_Merged")
markers2 <- get_markers(seurat_obj2, "Cell_Type_Compare_Merged")

# → 1-A. 直接存檔（完整 DEG 結果）
write.csv(markers1,
          file = file.path(outdir, paste0(Name_Seurat_obj1, "_DEG_markers.csv")),
          row.names = FALSE)
write.csv(markers2,
          file = file.path(outdir, paste0(Name_Seurat_obj2, "_DEG_markers.csv")),
          row.names = FALSE)

# ---------- 2. 組 marker gene list ----------
gene_list_from_markers <- function(markers_df){
  split(markers_df$gene, markers_df$cluster)
}

gene_sets1 <- gene_list_from_markers(markers1)
gene_sets2 <- gene_list_from_markers(markers2)

# 只保留兩物件都存在的 Cell Type
common_ct <- intersect(names(gene_sets1), names(gene_sets2))
if (length(common_ct) == 0)
  stop("兩個 Seurat 物件沒有任何共同的 Cell_Type_Compare_Merged。")

gene_sets1 <- gene_sets1[common_ct]
gene_sets2 <- gene_sets2[common_ct]

# ---------- 3. 畫 VennDiagram + 輸出基因清單 ----------
# --------- 在 3. 畫 VennDiagram + 輸出基因清單 的迴圈內 ----------

venn_plots <- lapply(common_ct, function(ct){
  
  # --- A. 取基因集合 ---
  genes_a  <- unique(gene_sets1[[ct]])
  genes_b  <- unique(gene_sets2[[ct]])
  inter_g  <- intersect(genes_a, genes_b)
  only_a_g <- setdiff(genes_a, genes_b)
  only_b_g <- setdiff(genes_b, genes_a)
  
  # === 先把 ct 變成安全檔名 ===
  safe_ct <- gsub("[\\s/\\\\:*?\"<>|]", "_", ct)   # 全轉成底線
  
  # --- B. 存文字檔 ---
  writeLines(only_a_g,
             file.path(outdir, sprintf("%s_%s_only.txt", safe_ct, Name_Seurat_obj1)))
  writeLines(inter_g,
             file.path(outdir, sprintf("%s_Intersection.txt", safe_ct)))
  writeLines(only_b_g,
             file.path(outdir, sprintf("%s_%s_only.txt", safe_ct, Name_Seurat_obj2)))
  
  # --- C. 畫 Venn + 存 PNG ---
  total_pop <- length(union(genes_a, genes_b))
  venn.plot <- venn.diagram(
    x              = setNames(list(genes_a, genes_b),
                              c(Name_Seurat_obj1, Name_Seurat_obj2)),
    filename       = NULL,
    category.names = c(Name_Seurat_obj1, Name_Seurat_obj2),
    fill           = c(Color_Seurat_obj1, Color_Seurat_obj2),
    alpha          = 0.60,
    cex            = 1.2,
    cat.cex        = 1,
    margin         = 0.05,
    main           = ct,
    print.mode     = c("raw", "percent"),
    total.population = total_pop
  )
  
  png(file.path(outdir, sprintf("%s_Venn.png", safe_ct)),
      width = 1200, height = 800, res = 120)
  grid::grid.draw(venn.plot)
  dev.off()
  
  # --- D. 回傳 grob 供排版 ---
  grid::grid.grabExpr(grid::grid.draw(venn.plot))
})


gridExtra::grid.arrange(grobs = venn_plots,
                        ncol  = min(5, length(venn_plots)))

# ---------- 4. 輸出所有 Venn 圖至單一 PDF ----------
pdf(file.path(outdir, "All_CellTypes_Venn.pdf"),
    width = 12, height = 8)
gridExtra::grid.arrange(grobs = venn_plots,
                        ncol  = min(5, length(venn_plots)))
dev.off()

message("✅  完成！所有結果已輸出至資料夾： ", normalizePath(outdir))
