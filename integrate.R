
library(dplyr)
library(harmony)
library(Seurat)
library(ggplot2)


# Initialize the Seurat object with the merged BPA nes

load('data/gestation/mouse/BPA_GO_merge_mouse_ges.rda')
load('data/lactation/mouse/BPA_GO_merge_mouse_lac.rda')
load('data/BPA_GO_comb.rda')
meta_data = meta_comb


dup = meta_data[duplicated(meta_data$cell_id),]


#meta_filter = meta_data[which(meta_data$dataset %in% c('L1','L2','G1','G2')),]

nes_sparse = as(nes_comb, "sparseMatrix")
nes_sparse = as(nes_merge, "sparseMatrix")
#nes_sparse = as(nes_comb, "sparseMatrix")
seurat <- CreateSeuratObject(counts = nes_sparse)
seurat<- FindVariableFeatures(seurat, selection.method = "vst")
all.bpa <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.bpa)

#meta_data = meta_comb
# Add meta data 
seurat[['cell.id']] = meta_data$cell_id
seurat[['dataset']] = meta_data$dataset

# Run PCA
seurat <- RunPCA(seurat)

ElbowPlot(seurat)



### Integration
# Extract the PCs
pc = as.matrix(attr(seurat[['pca']], "cell.embeddings"))

# Use first 15 PCs to do integration
V <- pc[,1:10]
d = cbind(meta_data,V)

p1 = ggplot() + geom_point(data = d, aes(x=PC_1, y =PC_2, col=dataset)) + theme(aspect.ratio=1)
p1

harmony_embeddings <- harmony::HarmonyMatrix(
  V, meta_data, 'dataset', do_pca = FALSE, verbose=FALSE
)

d2 = cbind(meta_data,harmony_embeddings)

p2 = ggplot() + geom_point(data = d2, aes(x=PC_1, y =PC_2, col=dataset)) + theme(aspect.ratio=1)
p2


save(harmony_embeddings, file = "data/harmony_comb.rda")



# Clustering
load('data/gestation/mouse/harmony_CP_Mm_ges.rda')

# Replace the cell embeddings as harmony embeddings
#cell_emb = attr(seurat[['pca']], "cell.embeddings")
attr(seurat[['pca']], "cell.embeddings") <- harmony_embeddings

#attr(seurat@reductions[["harmony"]], "harmony.embeddings")<- harmony_embeddings


seurat <- FindNeighbors(seurat,dims = 1:10, k.param = 20)
seurat <- FindClusters(seurat, resolution = 0.4)

head(Idents(seurat), 5)

seurat <- RunUMAP(seurat, dims = 1:10)
#DimPlot(seurat, reduction = "umap")

umap = seurat@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(meta_data) %>% cbind(seurat@meta.data$seurat_clusters)
colnames(umap)[colnames(umap) == 'seurat@meta.data$seurat_clusters'] = 'cluster'

#umap = drop_na(umap)



cluster = as.character(seurat@meta.data[["seurat_clusters"]])
colourCount = length(unique(cluster))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

p1_umap = ggplot(umap) + geom_point( aes(x = UMAP_1, y = UMAP_2, col = cluster),size=0.5) + theme_classic()+ theme(aspect.ratio=1 ) + guides(color = guide_legend(override.aes = list(size = 5), title = 'Cluster')) + scale_color_manual(values = getPalette(colourCount)) 
p1_umap


p2_umap = ggplot(umap, aes(x = UMAP_1, y = UMAP_2, col = dataset)) + geom_point(size=0.5)+ theme_classic() + theme(aspect.ratio=1) + guides(color = guide_legend(override.aes = list(size = 5))) 
p2_umap

p3_umap = ggplot(umap, aes(x = UMAP_1, y = UMAP_2, col = SuperCluster)) + geom_point(size=0.5) + theme_classic()+ theme(aspect.ratio=1 ) + guides(color = guide_legend(override.aes = list(size = 5), title = 'Reference_sc')) 
p3_umap

p4_umap = ggplot(umap, aes(x = UMAP_1, y = UMAP_2, col = CellType)) + geom_point(size=0.5) + theme_classic()+ theme(aspect.ratio=1 ) + guides(color = guide_legend(override.aes = list(size = 5), title = 'Reference_sn')) 
p4_umap



umap_hsd = umap[which(umap$SuperCluster=='Hsd'),]
umap_hsd$Colors = ifelse(umap_hsd$SuperCluster=='Hsd','#000000','#D3D3D3')

p3_umap = ggplot(umap_hsd, aes(x = UMAP_1, y = UMAP_2, col = Colors)) + geom_point(size=0.5) + theme(aspect.ratio=1 ) + guides(color = guide_legend(override.aes = list(size = 5))) 
p3_umap  




load("data/gestation/mouse/seurat_CP_Mm_ges.rda")
meta_data = seurat@meta.data

save(seurat,umap, file = "data/gestation/mouse/seurat_CP_Mm_ges.rda")



attr(seurat,'meta.data') <- meta_data
