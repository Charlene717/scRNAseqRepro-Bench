# 載入所需套件
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 讀取兩個 Seurat 物件 --------------------------------------------------

# 路徑一 (2025051906YDG.rds)
seurat_obj1 <- readRDS("C:/Charlene/Code_GitHub_BioInport2025/KGD_CellTypeAnnot_Skin/Export_2025051906YDG_Keloid_Charlene_2Step_RefComb_IntHarmony/2025051906YDG.rds")
# seurat_obj1 <- readRDS("C:/Charlene/Code_GitHub_BioInport2025/KGD_Lab_Code_scRNA-seq_Ch/Export_2025051707YXL_Keloid_For_Reproduction/2025051707YXL_seuratObject_Sample_seurat.rds")

seurat_obj1$Cell_Type_Compare <- seurat_obj1$CellType_SeuratClusters_ChatGPTDR

# 路徑二 (TNtype.combined-37.rds)
seurat_obj2 <- readRDS("C:/Charlene/Dataset_KGD_Lab/#_Seurat_Object/Keloid_Jojie/TNtype.combined-37.rds")
seurat_obj2$Cell_Type_KGD <- Idents(seurat_obj2)
seurat_obj2$Cell_Type_Compare <- Idents(seurat_obj2)

# 2. 建立合併子群的規則 ------------------------------------------------------
# 下面以 case_when() 為例，示範如何將多個子群合併為更大的類別。
# 請根據實際資料的子群名稱進行對應調整。

seurat_obj1$Cell_Type_Compare_Merged <- dplyr::case_when(
  seurat_obj1$Cell_Type_Compare %in% c("Pro-inflammatory fibroblast_1", 
                                       "Pro-inflammatory fibroblast_2") ~ "Pro-inflammatory fibroblast",
  seurat_obj1$Cell_Type_Compare %in% c("Secretory reticular fibroblast_1",
                                       "Secretory reticular fibroblast_2") ~ "Secretory reticular fibroblast",
  seurat_obj1$Cell_Type_Compare %in% c("Fibroblasts_1", "Fibroblasts_2", "Fibroblasts_3", "Fibroblasts_4") ~ "Fibroblasts",
  seurat_obj1$Cell_Type_Compare %in% c("Spinous keratinocyte_1", "Spinous keratinocyte_2",
                                       "Spinous keratinocyte_3", "Spinous keratinocyte_4", "Spinous keratinocyte_5") ~ "Spinous keratinocyte",
  # 依你的需求持續擴充或修改... 
  TRUE ~ seurat_obj1$Cell_Type_Compare  # 若不在上述對應中，則保留原標註
)

# 對 seurat_obj2 做同樣的處理
seurat_obj2$Cell_Type_Compare_Merged <- dplyr::case_when(
  seurat_obj2$Cell_Type_Compare %in% c("Pro-inflammatory fibroblast_1", 
                                       "Pro-inflammatory fibroblast_2") ~ "Pro-inflammatory fibroblast",
  seurat_obj2$Cell_Type_Compare %in% c("Secretory reticular fibroblast_1",
                                       "Secretory reticular fibroblast_2") ~ "Secretory reticular fibroblast",
  # 依需求增加合併規則...
  TRUE ~ seurat_obj2$Cell_Type_Compare
)

# 3. 檢視新的合併欄位 --------------------------------------------------------
# 查看新的合併後子群分布
unique(seurat_obj1$Cell_Type_Compare_Merged)
unique(seurat_obj2$Cell_Type_Compare_Merged)

# 也可以各自統計一下
table1_merged <- seurat_obj1@meta.data %>%
  count(Cell_Type_Compare_Merged, name = "CellCount")
table2_merged <- seurat_obj2@meta.data %>%
  count(Cell_Type_Compare_Merged, name = "CellCount")

print(table1_merged)
print(table2_merged)

# 4. 後續可視化 --------------------------------------------------------------
# 可以用合併後的欄位做分組繪圖 (group.by)
VlnPlot(seurat_obj1, features = c("nFeature_RNA", "nCount_RNA"), group.by = "Cell_Type_Compare_Merged", ncol = 2) +
  ggtitle("seurat_obj1 QC (Merged)")
VlnPlot(seurat_obj2, features = c("nFeature_RNA", "nCount_RNA"), group.by = "Cell_Type_Compare_Merged", ncol = 2) +
  ggtitle("seurat_obj2 QC (Merged)")

DimPlot(seurat_obj1, reduction = "umap", group.by = "Cell_Type_Compare_Merged") +
  ggtitle("UMAP - seurat_obj1 (Merged)")
DimPlot(seurat_obj2, reduction = "umap", group.by = "Cell_Type_Compare_Merged") +
  ggtitle("UMAP - seurat_obj2 (Merged)")
