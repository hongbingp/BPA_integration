library(Seurat)

load('data/lactation/mouse/matrix_lac_mouse.rda')

# Initialize the Seurat object with the raw

seurat_L1 <- CreateSeuratObject(counts = m_sc1)
seurat_L2 <- CreateSeuratObject(counts = m_sc2)
seurat_Lsn <- CreateSeuratObject(counts = m_sn1)

seurat_L1[["percent.mt"]] <- PercentageFeatureSet(seurat_L1, pattern = "^Mt")
#head(seurat_L1@meta.data, 30)
seurat_L2[["percent.mt"]] <- PercentageFeatureSet(seurat_L2, pattern = "^Mt")
seurat_Lsn[["percent.mt"]] <- PercentageFeatureSet(seurat_Lsn, pattern = "^Mt")



load('data/gestation/mouse/matrix.rda')

seurat_G1 <- CreateSeuratObject(counts = m_sc1)
seurat_G2 <- CreateSeuratObject(counts = m_sc2)
seurat_Gsn <- CreateSeuratObject(counts = m_sn1)

seurat_G1[["percent.mt"]] <- PercentageFeatureSet(seurat_G1, pattern = "^Mt")
#head(seurat_G1@meta.data, 30)
seurat_G2[["percent.mt"]] <- PercentageFeatureSet(seurat_G2, pattern = "^Mt")
seurat_Gsn[["percent.mt"]] <- PercentageFeatureSet(seurat_Gsn, pattern = "^Mt")


# Visualize QC metrics as a violin plot
VlnPlot(seurat_L1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(seurat_Lsn, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_Lsn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(seurat_Gsn, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_Gsn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



# QC

pmt = seurat_L1@meta.data[["percent.mt"]]
mad(seurat_L1@meta.data[["percent.mt"]],constant = 1)

seurat_L1 <- subset(seurat_L1, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 5)
seurat_L2 <- subset(seurat_L2, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 5)
seurat_Lsn <- subset(seurat_Lsn, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 0.7)

seurat_G1 <- subset(seurat_G1, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 5)
seurat_G2 <- subset(seurat_G2, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 5)
seurat_Gsn <- subset(seurat_Gsn, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 0.7)

# Get the matrix after quality control
m_L1_qc = attr(seurat_L1[['RNA']], 'counts')
m_L2_qc = attr(seurat_L2[['RNA']], 'counts')
m_Lsn_qc = attr(seurat_Lsn[['RNA']], 'counts')

m_G1_qc = attr(seurat_G1[['RNA']], 'counts')
m_G2_qc = attr(seurat_G2[['RNA']], 'counts')
m_Gsn_qc = attr(seurat_Gsn[['RNA']], 'counts')

save(m_L1_qc, m_L2_qc, m_Lsn_qc,m_G1_qc, m_G2_qc, m_Gsn_qc, file = "data/matrix_Mm_qc.rda")
