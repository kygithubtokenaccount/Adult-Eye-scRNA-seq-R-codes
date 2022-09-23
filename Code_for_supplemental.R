# R code for Supplemental Figures

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(sctransform)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

# Supp Fig 1
# load monocle3 CDS object of adult 1 day, 3 days and 7 days old male adult eyes as "a17m"
cell_type_color <- c("R1-6" = "#f263df",
                     "R7" = "#639dfc",
                     "R8" = "#1bbec3",
                     "Cone" = "#17b93f",
                     "1_pigments" = "#b69e1c",
                     "2,3_pigments" = "#f6756f")
plot_cells(a17m, color_cells_by="assigned_cell_type", label_cell_groups = F, show_trajectory_graph = F, cell_size = 1)+scale_color_manual(values = cell_type_color) + NoAxes() + NoLegend() #S1A
plot_cells(a17m, color_cells_by = "orig.ident", label_cell_groups = F, show_trajectory_graph = F, cell_size = 1)+NoAxes()+NoLegend() #S1B


# load 1 day, 3 days, 7 days old male adult eyes Seurat analyzed data sets as a1d, a3d and a7d, respectively

# Supp Fig 2
DimPlot(a1d, label = F, pt.size = 2, cols = c("grey","grey","grey","red","grey","grey")) + NoAxes() + NoLegend() #S2 A
DimPlot(a3d, label = F, pt.size = 2, cols = c("grey","grey","grey","red","grey","grey")) + NoAxes() + NoLegend() #S2 A
FeaturePlot(a1d, features = "sens", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 B
FeaturePlot(a3d, features = "sens", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 B
FeaturePlot(a1d, features = "Rh5", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 C
FeaturePlot(a3d, features = "Rh5", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 C
FeaturePlot(a1d, features = "Rh6", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 D
FeaturePlot(a3d, features = "Rh6", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 D
FeaturePlot(a1d, features = "CG2082", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 E
FeaturePlot(a3d, features = "CG2082", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S2 E

# Supp Fig 3
DimPlot(a1d, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","red","grey")) + NoAxes() + NoLegend() #S3 A
DimPlot(a3d, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","red","grey")) + NoAxes() + NoLegend() #S3 A
FeaturePlot(a1d, features = "pros", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 B
FeaturePlot(a3d, features = "pros", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 B
FeaturePlot(a1d, features = "Rh3", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 C
FeaturePlot(a3d, features = "Rh3", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 C
FeaturePlot(a1d, features = "Rh4", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 D
FeaturePlot(a3d, features = "Rh4", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 D
FeaturePlot(a1d, features = "Cep135", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 E
FeaturePlot(a3d, features = "Cep135", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S3 E

# Supp Fig 4
# load 1 day, 3 days and 7 days old male adult eyes with full Rh, Dorsal 3rd R7 and DRA annotated as a1drhc, a3drhc, and a7drhc, respectively
DimPlot(a1drhc, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","red","grey")) + NoAxes() + NoLegend() #S4 A
DimPlot(a3drhc, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","red","grey")) + NoAxes() + NoLegend() #S4 A
DimPlot(a7drhc, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","red","grey")) + NoAxes() + NoLegend() #S4 A
FeaturePlot(a1drhc, features = "Rh3", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 B
FeaturePlot(a3drhc, features = "Rh3", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 B
FeaturePlot(a7drhc, features = "Rh3", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 B
FeaturePlot(a1drhc, features = "hth", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 C
FeaturePlot(a3drhc, features = "hth", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 C
FeaturePlot(a7drhc, features = "hth", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 C
FeaturePlot(a1drhc, features = "Skeletor", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 D
FeaturePlot(a3drhc, features = "Skeletor", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 D
FeaturePlot(a7drhc, features = "Skeletor", pt.size = 2, order = F) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S4 D

# Supp Fig 5
DimPlot(a1d, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","grey","red")) + NoAxes() + NoLegend() #S5 A
DimPlot(a3d, label = F, pt.size = 2, cols = c("grey","grey","grey","grey","grey","red")) + NoAxes() + NoLegend() #S5 A
FeaturePlot(a1d, features = "ninaE", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S5 B
FeaturePlot(a3d, features = "ninaE", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S5 B
FeaturePlot(a1d, features = "oc", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S5 C
FeaturePlot(a3d, features = "oc", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S5 C
FeaturePlot(a1d, features = "Zasp66", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S5 D
FeaturePlot(a3d, features = "Zasp66", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S5 D
VlnPlot(a1d, features = "ninaE") + NoLegend()
VlnPlot(a3d, features = "ninaE") + NoLegend()

# Supp Fig 6
DimPlot(a1d, label = F, pt.size = 3, cols = c("grey","grey","red","grey","grey","grey")) + NoAxes() + NoLegend() #S6 A
DimPlot(a3d, label = F, pt.size = 3, cols = c("grey","grey","red","grey","grey","grey")) + NoAxes() + NoLegend() #S6 A
FeaturePlot(a1d, features = "ct", pt.size = 3) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S6 B
FeaturePlot(a3d, features = "ct", pt.size = 3) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S6 B
FeaturePlot(a1d, features = "Crys", pt.size = 3) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S6 C
FeaturePlot(a3d, features = "Crys", pt.size = 3) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S6 C
FeaturePlot(a1d, features = "CG5597", pt.size = 3) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S6 D
FeaturePlot(a3d, features = "CG5597", pt.size = 3) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S6 D

# Supp Fig 7
DimPlot(a1d, label = F, pt.size = 2, cols = c("red","red","grey","grey","grey","grey")) + NoAxes() + NoLegend() #S7A
DimPlot(a3d, label = F, pt.size = 2, cols = c("red","red","grey","grey","grey","grey")) + NoAxes() + NoLegend() #S7A
FeaturePlot(a1d, features = "w", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 B
FeaturePlot(a3d, features = "w", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 B
FeaturePlot(a1d, features = "Pdh", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 C
FeaturePlot(a3d, features = "Pdh", pt.size = 2) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 C
VlnPlot(a1d, features = "Pdh") + NoLegend() #S7 D
VlnPlot(a3d, features = "Pdh") + NoLegend() #S7 E
FeaturePlot(a1d, features = "santa-maria", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 F
FeaturePlot(a3d, features = "santa-maria", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 F
FeaturePlot(a1d, features = "wrapper", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 G
FeaturePlot(a3d, features = "wrapper", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S7 G

# Supp Fig 8
# load in the following seurat objects
# a1dr56 = 1 day old data with Rh5 and Rh6 removed
# a3dr56 = 1 day old data with Rh5 and Rh6 removed
# a7dr56 = 1 day old data with Rh5 and Rh6 removed
FeaturePlot(a1dr56, features = "Rh3", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S8 A
FeaturePlot(a3dr56, features = "Rh3", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S8 A
FeaturePlot(a7dr56, features = "Rh3", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S8 A
FeaturePlot(a1dr56, features = "Rh4", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S8 B
FeaturePlot(a3dr56, features = "Rh4", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S8 B
FeaturePlot(a7dr56, features = "Rh4", pt.size = 2, order = T) + NoAxes() + NoLegend() & theme(plot.title= element_blank()) #S8 B
