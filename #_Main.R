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


# 2. 基本資訊檢查 (QC) -------------------------------------------------------

# 查看物件概況
# 例如：細胞數量、基因數量、metadata 欄位等
seurat_obj1
seurat_obj2

# 查看 metadata 欄位
colnames(seurat_obj1@meta.data)
colnames(seurat_obj2@meta.data)

# 可以檢查常見的 QC 指標 (nFeature_RNA、nCount_RNA、percent.mt 等)
# 如果物件中尚未計算 percent.mt，可在讀檔後自行計算
# 這裡示範直接繪製小提琴圖比較分佈，若物件中有相同的欄位名
VlnPlot(seurat_obj1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) +
  ggtitle("seurat_obj1 QC")
VlnPlot(seurat_obj2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) +
  ggtitle("seurat_obj2 QC")


# 3. 比較細胞標註（Cell Type）相關資訊 --------------------------------------

# 若兩個物件都有類似的欄位存放 Cell Type，例如：
# seurat_obj1$Cell_Type_Compare
# seurat_obj2$Cell_Type_Compare
# 可使用 unique() 查看有哪些類型：

unique(seurat_obj1$Cell_Type_Compare)
unique(seurat_obj2$Cell_Type_Compare)

# 若欄位名稱不同，請替換成對應的 metadata 欄位

# 也可以加上表格統計 (dplyr::count) 看每種 Cell Type 有多少細胞
table1 <- seurat_obj1@meta.data %>%
  count(Cell_Type_Compare, name = "CellCount")

table2 <- seurat_obj2@meta.data %>%
  count(Cell_Type_Compare, name = "CellCount")

print(table1)
print(table2)


# 4. 後續可視化與進階比較 (選擇性) ------------------------------------------

# (1) UMAP / t-SNE 可視化
# 如果物件中已經計算好 UMAP/t-SNE，可直接繪製：
DimPlot(seurat_obj1, reduction = "umap", group.by = "Cell_Type_Compare") +
  ggtitle("UMAP - seurat_obj1")

DimPlot(seurat_obj2, reduction = "umap", group.by = "Cell_Type_Compare") +
  ggtitle("UMAP - seurat_obj2")

# (2) Marker 基因的表現比較
# 可以選擇幾個 Marker 基因，查看兩個物件中表現情況：
markers_to_check <- c("COL1A1", "ACTA2", "KRT14")  # 範例基因，實際依你需要修改

FeaturePlot(seurat_obj1, features = markers_to_check, reduction = "umap", ncol = 2)
FeaturePlot(seurat_obj2, features = markers_to_check, reduction = "umap", ncol = 2)

# (3) 如果需要進行整合分析或合併比較，請參考 Seurat integration 或 merge 流程
# 例如：
# combined_obj <- merge(seurat_obj1, seurat_obj2)
# 之後可依分析需求執行標準流程 (NormalizeData, FindVariableFeatures, ScaleData, RunPCA, RunUMAP, etc.)

