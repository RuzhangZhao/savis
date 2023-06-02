# savis


Single-cell RNAseq Adaptive Visualization 

```R
library(Seurat)
library(savis)
library(mixhvg)
npcs=30
obj<-CreateSeuratObject(expr,verbose=F)
obj$cell_label<-cell_label
obj<-NormalizeData(obj,verbose=F)
obj<-FindVariableFeaturesMix(obj,verbose = F)
obj<-ScaleData(obj,verbose=F)
suppressWarnings(obj<-RunPCA(obj,npcs=npcs,verbose=F))
obj<-RunUMAP(obj,dims = 1:npcs,verbose=F)
obj<-RunPreSAVIS(obj,verbose=F)
obj<-RunSAVIS(obj,verbose=F)
adaDimPlot(obj,group.by = cell_label,color.mode = 0)

```

