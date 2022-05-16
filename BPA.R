setwd('/media/jiahui/hdd2/seq_integration/')

library(viper)
library(ggplot2)


# CP
load("CP-Mm-MSigDB-regulon.rda")
load("data/matrix_Mm_qc.rda")
# Filter the BPA
len <- unlist(lapply(gset, function(x) length(x$tfmode)))

len_df = data.frame(len)
len_df_cut = data.frame(len_df[(len < 250),])
colnames(len_df_cut) = 'len'
p1 = ggplot(len_df_cut,aes(x=len)) + geom_histogram(binwidth=10)
p1

gset <- gset[len >= 50 & len <= 100]

# G1
m1 = as.matrix(m_G1_qc)
#m <- m_sc1[apply(m_sc1, 1, sd) > 0, ]
nes_G1 <- aREA(m1, gset)$nes
save(nes_G1, file = "data/gestation/mouse/BPA_CP_G1.rda")

# G2
m2 = as.matrix(m_G2_qc)
nes_G2 <- aREA(m2, gset)$nes
save(nes_G2, file = "data/gestation/mouse/BPA_CP_G2.rda")

# Gsn
m3 = as.matrix(m_Gsn_qc)
nes_Gsn <- aREA(m3, gset)$nes
save(nes_Gsn, file = "data/gestation/mouse/BPA_CP_Gsn.rda")




# CP human
load("CP-Hs-MSigDB-regulon.rda")
load("data/lactation/human/sampledmatrix_lac_human.rda")
# Filter the BPA
len <- unlist(lapply(gset, function(x) length(x$tfmode)))

len_df = data.frame(len)
len_df_cut = data.frame(len_df[(len < 250),])
colnames(len_df_cut) = 'len'
p1 = ggplot(len_df_cut,aes(x=len)) + geom_histogram(binwidth=10)
p1

gset <- gset[len >= 50 & len <= 100]

# lmc1
m1 = as.matrix(m_lmc1_qc)
#m <- m_lmc1[apply(m_lmc1, 1, sd) > 0, ]
nes1 <- aREA(m1, gset)$nes
save(nes1, file = "data/lactation/human/BPA_CP_lmc1.rda")

# lmc2
m2 = as.matrix(m_lmc2_qc)
nes2 <- aREA(m2, gset)$nes
save(nes2, file = "data/lactation/human/BPA_CP_lmc2.rda")

# lmc3
m3 = as.matrix(m_lmc3_qc)
nes3 <- aREA(m3, gset)$nes
save(nes3, file = "data/lactation/human/BPA_CP_lmc3.rda")

# lmc4
m4 = as.matrix(m_lmc4_qc)
nes4 <- aREA(m4, gset)$nes
save(nes4, file = "data/lactation/human/BPA_CP_lmc4.rda")









# GO
load("GO-BP-Mm-MSigDB-regulon.rda")
load("data/matrix_Mm_qc.rda")
#m <- m_sc1[apply(m_sc1, 1, sd) > 0, ]
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]

# G1
m1 = as.matrix(m_G1_qc)
#m <- m_sc1[apply(m_sc1, 1, sd) > 0, ]
nes_G1 <- aREA(m1, gset)$nes
save(nes_G1, file = "data/gestation/mouse/BPA_GO_G1.rda")

# G2
m2 = as.matrix(m_G2_qc)
nes_G2 <- aREA(m2, gset)$nes
save(nes_G2, file = "data/gestation/mouse/BPA_GO_G2.rda")

# Gsn
m3 = as.matrix(m_Gsn_qc)
nes_Gsn <- aREA(m3, gset)$nes
save(nes_Gsn, file = "data/gestation/mouse/BPA_GO_Gsn.rda")


# L1
m1 = as.matrix(m_L1_qc)
#m <- m_sc1[apply(m_sc1, 1, sd) > 0, ]
nes_L1 <- aREA(m1, gset)$nes
save(nes_L1, file = "data/lactation/mouse/BPA_GO_L1.rda")

# L2
m2 = as.matrix(m_L2_qc)
nes_L2 <- aREA(m2, gset)$nes
save(nes_L2, file = "data/lactation/mouse/BPA_GO_L2.rda")

# Lsn
m3 = as.matrix(m_Lsn_qc)
nes_Lsn <- aREA(m3, gset)$nes
save(nes_Lsn, file = "data/lactation/mouse/BPA_GO_Lsn.rda")






# GO human
load("GO-BP-Hs-MSigDB-regulon.rda")
load("data/lactation/human/sampledmatrix_lac_human.rda")
# Filter the BPA
len <- unlist(lapply(gset, function(x) length(x$tfmode)))

len_df = data.frame(len)
len_df_cut = data.frame(len_df[(len < 250),])
colnames(len_df_cut) = 'len'
p1 = ggplot(len_df_cut,aes(x=len)) + geom_histogram(binwidth=10)
p1

gset <- gset[len >= 50 & len <= 100]

# lmc1
m1 = as.matrix(m_lmc1_qc)
#m <- m_lmc1[apply(m_lmc1, 1, sd) > 0, ]
nes1 <- aREA(m1, gset)$nes
save(nes1, file = "data/lactation/human/BPA_GO_lmc1.rda")

# lmc2
m2 = as.matrix(m_lmc2_qc)
nes2 <- aREA(m2, gset)$nes
save(nes2, file = "data/lactation/human/BPA_GO_lmc2.rda")

# lmc3
m3 = as.matrix(m_lmc3_qc)
nes3 <- aREA(m3, gset)$nes
save(nes3, file = "data/lactation/human/BPA_GO_lmc3.rda")

# lmc4
m4 = as.matrix(m_lmc4_qc)
nes4 <- aREA(m4, gset)$nes
save(nes4, file = "data/lactation/human/BPA_GO_lmc4.rda")






