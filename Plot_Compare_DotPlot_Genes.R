# ── 載入套件 ────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ── 自訂函式：為單一 Seurat 物件建立 dotplot 資料 ───────────
prep_dot_data <- function(seu, genes, celltype_col, obj_name){
  # 1. 只取目前物件有的基因
  present_genes <- genes[genes %in% rownames(seu)]
  missing_genes <- setdiff(genes, present_genes)
  if(length(missing_genes) > 0){
    message(obj_name, "：缺少基因 ", paste(missing_genes, collapse = ", "))
  }
  
  # 2. 擷取表現資料
  expr_mat <- FetchData(seu, vars = present_genes, slot = "data")
  meta      <- seu@meta.data[[celltype_col]]
  
  # 3. 計算平均表現量與百分比
  dot_df <- expr_mat |>
    as_tibble() |>
    mutate(CellType = meta) |>
    pivot_longer(-CellType, names_to = "Gene", values_to = "Expr") |>
    group_by(CellType, Gene) |>
    summarise(avg_exp  = mean(Expr),
              pct_exp  = mean(Expr > 0) * 100,
              .groups = "drop") |>
    mutate(Object = obj_name)
  
  return(dot_df)
}

# ── 主要函式：畫氣泡圖 ──────────────────────────────────────
bubble_compare <- function(seu1, seu2,
                           genes   = c("PHLDA3","CSRP1","NAV1","DAM177B","NEDD4"),
                           celltype_col = "Cell_Type_Compare_Merged",
                           name1   = "Seurat obj1",
                           name2   = "Seurat obj2",
                           palette = c("#fee8c8","#fdbb84","#e34a33")){
  
  # 1. 準備資料
  df1 <- prep_dot_data(seu1, genes, celltype_col, name1)
  df2 <- prep_dot_data(seu2, genes, celltype_col, name2)
  plot_df <- bind_rows(df1, df2) |>
    mutate(Gene = factor(Gene, levels = genes))        # 固定基因順序
  
  # 2. 畫圖
  ggplot(plot_df,
         aes(x = CellType, y = Gene,
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
          strip.text = element_text(size = 13, face = "bold"))
}

# ── 範例呼叫 ────────────────────────────────────────────────
bubble_compare(seurat_obj1, seurat_obj2,
               genes = c("PHLDA3","CSRP1","NAV1","DAM177B","NEDD4"),
               name1 = Name_Seurat_obj1,
               name2 = Name_Seurat_obj2)
