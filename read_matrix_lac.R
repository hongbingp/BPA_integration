.libPaths("/home/chaochen/R/x86_64-pc-linux-gnu-library/4.1")

setwd('/media/jiahui/hdd2/seq_integration/data')

library(Seurat)

m_sc1 <- ReadMtx(
  mtx = "lactation/mouse/sc1/GSM2834502_L_1_matrix.mtx", features = "lactation/mouse/sc1/GSM2834502_L_1_genes.tsv",
  cells = "lactation/mouse/sc1/GSM2834502_L_1_barcodes.tsv"
)
m_sc2 <- ReadMtx(
  mtx = "lactation/mouse/sc2/GSM2834503_L_2_matrix.mtx", features = "lactation/mouse/sc2/GSM2834503_L_2_genes.tsv",
  cells = "lactation/mouse/sc2/GSM2834503_L_2_barcodes.tsv"
)
m_sn1 <- ReadMtx(
  mtx = "lactation/mouse/sn1/matrix.mtx", features = "lactation/mouse/sn1/features.tsv",
  cells = "lactation/mouse/sn1/barcodes.tsv"
)


save(m_sc1, m_sc2, m_sn1, file = "E:/seq_integration/data/lactation/mouse/matrix_lac_mouse.rda")



m_lmc1 <- ReadMtx(
  mtx = "lactation/human/E-MTAB-9841/LMC1_matrix.mtx", features = "lactation/human/E-MTAB-9841/LMC1_features.tsv",
  cells = "lactation/human/E-MTAB-9841/LMC1_barcodes.tsv"
)

m_lmc2 <- ReadMtx(
  mtx = "lactation/human/E-MTAB-9841/LMC2_matrix.mtx", features = "lactation/human/E-MTAB-9841/LMC2_features.tsv",
  cells = "lactation/human/E-MTAB-9841/LMC2_barcodes.tsv"
)

m_lmc3 <- ReadMtx(
  mtx = "lactation/human/E-MTAB-9841/LMC3_matrix.mtx", features = "lactation/human/E-MTAB-9841/LMC3_features.tsv",
  cells = "lactation/human/E-MTAB-9841/LMC3_barcodes.tsv"
)

m_lmc4 <- ReadMtx(
  mtx = "lactation/human/E-MTAB-9841/LMC4_matrix.mtx", features = "lactation/human/E-MTAB-9841/LMC4_features.tsv",
  cells = "lactation/human/E-MTAB-9841/LMC4_barcodes.tsv"
)


save(m_lmc2, m_lmc3, file = "lactation/human/matrix_lac_human23.rda")
save(m_lmc1,m_lmc2, m_lmc3, m_lmc4,file = "lactation/human/matrix_lac_human.rda")
