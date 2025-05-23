##### Presetting ######
rm(list = ls()) # Clean all variables ##* Comment out if running the entire script

## Speed up (Old Version)
memory.limit(150000)

# ## Speed up
# if(!require('future')) install.packages('future'); library(future)
# ## https://github.com/immunogenomics/presto
# if(!require('presto')) devtools::install_github("immunogenomics/presto"); library(presto) # Speeds up FindAllMarkers
# plan(multicore, workers = 20)
# options(future.globals.maxSize = 2048*100 * 1024^2) # Set memory limit to ~204.8 GB

#### Load Packages ####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}


#### Load Data ####
# 路徑一 (2025051906YDG.rds)
seurat_obj1 <- readRDS("C:/Charlene/Code_GitHub_BioInport2025/KGD_CellTypeAnnot_Skin/Export_2025051906YDG_Keloid_Charlene_2Step_RefComb_IntHarmony/2025051906YDG.rds")
# seurat_obj1 <- readRDS("C:/Charlene/Code_GitHub_BioInport2025/KGD_Lab_Code_scRNA-seq_Ch/Export_2025051707YXL_Keloid_For_Reproduction/2025051707YXL_seuratObject_Sample_seurat.rds")

seurat_obj1$Cell_Type_Compare <- seurat_obj1$CellType_SeuratClusters_ChatGPTDR

# 路徑二 (TNtype.combined-37.rds)
seurat_obj2 <- readRDS("C:/Charlene/Dataset_KGD_Lab/#_Seurat_Object/Keloid_Jojie/TNtype.combined-37.rds")
seurat_obj2$Cell_Type_KGD <- Idents(seurat_obj2)
seurat_obj2$Cell_Type_Compare <- Idents(seurat_obj2)


#### Rename cell type ####
# 請先載入所需套件
library(stringr)

merge_cell_type <- function(vec){
  vec |>
    str_remove("([_-]?\\d+)$") |>              # 1. 去掉尾端連接符＋數字
    str_trim() |>                              # 2. 去空白
    str_to_lower() |>                          # 3. 全轉小寫
    str_replace("([a-z])cells$", "\\1 cells") |># 4. 若結尾為 ...cells 且前面緊貼字母，補空格
    str_to_sentence() |>                       # 5. 首字母大寫
    (\(x) ifelse(str_ends(x, regex("s$", ignore_case = TRUE)),
                 x,
                 paste0(x, "s")))()            # 6. 若結尾不是 s 再補 s
}

seurat_obj1$Cell_Type_Compare_Merged <- merge_cell_type(seurat_obj1$Cell_Type_Compare)
seurat_obj2$Cell_Type_Compare_Merged <- merge_cell_type(seurat_obj2$Cell_Type_Compare)



# ── 專門統一 Keratinocytes / Fibroblasts ────────────────────
library(stringr)
std_kera_fibro <- function(vec){
  vec |>
    str_replace(regex(".*keratinocytes?$", TRUE), "Keratinocytes") |>
    str_replace(regex(".*fibroblasts?$",  TRUE), "Fibroblasts")
}

seurat_obj1$Cell_Type_Compare_Merged <- std_kera_fibro(seurat_obj1$Cell_Type_Compare_Merged)
seurat_obj2$Cell_Type_Compare_Merged <- std_kera_fibro(seurat_obj2$Cell_Type_Compare_Merged)


#### Compare ####
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



#### Visualization ####



# 可以用合併後的欄位做分組繪圖 (group.by)
VlnPlot(seurat_obj1, features = c("nFeature_RNA", "nCount_RNA"), group.by = "Cell_Type_Compare_Merged", ncol = 2) +
  ggtitle("seurat_obj1 QC (Merged)")
VlnPlot(seurat_obj2, features = c("nFeature_RNA", "nCount_RNA"), group.by = "Cell_Type_Compare_Merged", ncol = 2) +
  ggtitle("seurat_obj2 QC (Merged)")

DimPlot(seurat_obj1, reduction = "umap", group.by = "Cell_Type_Compare_Merged") +
  ggtitle("UMAP - seurat_obj1 (Merged)")
DimPlot(seurat_obj2, reduction = "umap", group.by = "Cell_Type_Compare_Merged") +
  ggtitle("UMAP - seurat_obj2 (Merged)")
