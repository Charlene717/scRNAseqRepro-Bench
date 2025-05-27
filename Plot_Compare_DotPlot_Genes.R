# ── 載入套件（同前） ───────────────────────────────────────
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. 建立 dotplot 數據 (單一 Seurat) ─────────────────────
prep_dot_data <- function(seu, genes,
                          celltype_col,
                          add_col   = NULL,   # ★ 新增：額外欄位
                          combine   = TRUE,   # ★ 新增：是否組合兩欄
                          obj_name){
  
  # ---- 1.1 基因檢查 ----
  present_genes <- genes[genes %in% rownames(seu)]
  missing_genes <- setdiff(genes, present_genes)
  if(length(missing_genes) > 0){
    message(obj_name, "：缺少基因 ", paste(missing_genes, collapse = ", "))
  }
  
  # ---- 1.2 建立分群變數 ----
  meta_ct <- seu[[celltype_col]][,1]
  
  if(!is.null(add_col)){
    if(!add_col %in% colnames(seu@meta.data)){
      stop("欄位 '", add_col, "' 不存在於 meta.data")
    }
    meta_add <- seu[[add_col]][,1]
    if(combine){
      group_var <- paste(meta_ct, meta_add, sep = " | ")
    } else {
      group_var <- interaction(meta_ct, meta_add, sep = " | ")
    }
  } else {
    group_var <- meta_ct
  }
  
  # ---- 1.3 計算平均與百分比 ----
  expr_mat <- FetchData(seu, vars = present_genes, slot = "data")
  dot_df <- expr_mat |>
    as_tibble() |>
    mutate(Group = group_var) |>
    pivot_longer(-Group, names_to = "Gene", values_to = "Expr") |>
    group_by(Group, Gene) |>
    summarise(avg_exp = mean(Expr),
              pct_exp = mean(Expr > 0) * 100,
              .groups = "drop") |>
    mutate(Object = obj_name)
  
  return(dot_df)
}

# ── 2. 主函式：bubble_compare ───────────────────────────────
bubble_compare <- function(seu1, seu2,
                           genes,
                           celltype_col = "Cell_Type_Compare_Merged",
                           add_col      = NULL,      # ★ 新增
                           combine      = TRUE,      # ★ 新增
                           name1        = "Seurat obj1",
                           name2        = "Seurat obj2",
                           palette      = c("#fee8c8","#fdbb84","#e34a33")){
  
  df1 <- prep_dot_data(seu1, genes, celltype_col, add_col, combine, name1)
  df2 <- prep_dot_data(seu2, genes, celltype_col, add_col, combine, name2)
  
  plot_df <- bind_rows(df1, df2) |>
    mutate(Gene = factor(Gene, levels = genes))
  
  ggplot(plot_df,
         aes(x = Group, y = Gene,
             size = pct_exp, colour = avg_exp)) +
    geom_point(alpha = .8) +
    scale_size(range = c(1, 10), name = "% Expressing") +
    scale_colour_gradientn(colours = palette, name = "Avg Expr") +
    coord_flip() +
    facet_wrap(~Object, ncol = 1) +
    labs(title = "Gene expression comparison",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text  = element_text(size = 13, face = "bold"))
}

bubble_compare(seurat_obj1, seurat_obj2,
               genes = c("PHLDA3","CSRP1","NAV1","DAM177B","NEDD4"),
               name1 = Name_Seurat_obj1,
               name2 = Name_Seurat_obj2)


##############################################################################
# 欲保留的細胞類型
# cell_keep <- c("Schwann cells", "Fibroblasts", "Keratinocytes")
cell_keep <- c("Schwann cells")


seurat_subset1 <- seurat_obj1[, seurat_obj1$Cell_Type_Compare_Merged %in% cell_keep]

seurat_subset2 <- seurat_obj2[, seurat_obj2$Cell_Type_Compare_Merged %in% cell_keep]


bubble_compare(seurat_subset1, seurat_subset2,
               genes        = c("PHLDA3","CSRP1","NAV1","DAM177B","NEDD4"),
               celltype_col = "Cell_Type_Compare_Merged",
               add_col      = "skin_group",   # <─ 新增第二欄
               combine      = TRUE,           # <─ 預設 TRUE，貼在一起
               name1        = Name_Seurat_obj1,
               name2        = Name_Seurat_obj2)

bubble_compare(
  seurat_obj1,
  seurat_obj2,
  genes        = c(
    "PHLDA3",
    "NRXN1",
    "CCN2",
    "MPZ",
    "PTN",
    "S100B",   # 許旺細胞經典 Ca2+ 結合蛋白
    "NGFR",    # p75^NTR，未成熟/修復型 Schwann cell 標誌
    "NES",     # Nestin，去分化/修復型上調
    "IGFBP5",  # 促纖維化分泌因子
    "IGFBP3",  # 調節生長因子結合蛋白
    "CCN3",    # NOV，細胞通訊網絡因子3
    "TNFAIP6", # TSG-6，調節炎症與基質
    "COL1A1",  # I 型膠原
    "COL3A1"   # III 型膠原
  ),
  celltype_col = "Cell_Type_Compare_Merged",
  # add_col      = "skin_group",   # <─ 新增第二欄
  combine      = FALSE,           # <─ 預設 TRUE，貼在一起
  name1        = Name_Seurat_obj1,
  name2        = Name_Seurat_obj2
)

