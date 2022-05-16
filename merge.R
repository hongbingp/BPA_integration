
# GO

load("data/lactation/mouse/BPA_GO_L1.rda")
load("data/lactation/mouse/BPA_GO_L2.rda")
load("data/lactation/mouse/BPA_GO_Lsn.rda")

load("data/gestation/mouse/BPA_GO_G1.rda")
load("data/gestation/mouse/BPA_GO_G2.rda")
load("data/gestation/mouse/BPA_GO_Gsn.rda")

int = intersect(rownames(nes_G1),rownames(nes_Gsn))
nes_G1 = nes_G1[int,]
nes_G2 = nes_G2[int,]
nes_Gsn = nes_Gsn[int,]

#diffcn = setdiff(rownames(nes1),rownames(nes3_ges))

#nes3t_del = nes3t[,-which(colnames(nes3t)==diffcn)]

# Create meta_data
meta1_ges = data.frame(cell_id = colnames(nes_G1),dataset='G1')
meta2_ges = data.frame(cell_id = colnames(nes_G2),dataset='G2')
meta3_ges = data.frame(cell_id = colnames(nes_Gsn),dataset='G_sn')
meta1_lac = data.frame(cell_id = colnames(nes_L1),dataset='L1')
meta2_lac = data.frame(cell_id = colnames(nes_L2),dataset='L2')
meta3_lac = data.frame(cell_id = colnames(nes_Lsn),dataset='L_sn')

# Merge
meta_data = rbind(meta1_ges,meta2_ges,meta3_ges,meta1_lac,meta2_lac,meta3_lac)
nes_merge = cbind(nes_G1,nes_G2,nes_Gsn,nes_L1,nes_L2,nes_Lsn)
save(nes_merge,meta_data,file='data/BPA_GO_merge_mouse.rda')

meta_data = rbind(meta1_ges,meta2_ges,meta1_lac,meta2_lac)
nes_merge = cbind(nes_G1,nes_G2,nes_L1,nes_L2)
save(nes_merge,meta_data,file='data/BPA_GO_merge_mouse_nosn.rda')

meta_data = rbind(meta1_ges,meta2_ges,meta3_ges)
nes_merge = cbind(nes_G1,nes_G2,nes_Gsn)
save(nes_merge,meta_data,file='data/gestation/mouse/BPA_CP_merge_mouse_ges.rda')

meta_data = rbind(meta1_lac,meta2_lac,meta3_lac)
nes_merge = cbind(nes_L1,nes_L2,nes_Lsn)
save(nes_merge,meta_data,file='data/lactation/mouse/BPA_CP_merge_mouse_lac.rda')







# CP

load("data/lactation/mouse/BPA_GO_L1.rda")
load("data/lactation/mouse/BPA_GO_L2.rda")
load("data/lactation/mouse/BPA_GO_Lsn.rda")

load("data/gestation/mouse/BPA_CP_G1.rda")
load("data/gestation/mouse/BPA_CP_G2.rda")
load("data/gestation/mouse/BPA_CP_Gsn.rda")

uni = union(rownames(nes_G1),rownames(nes_Lsn))
#diffcn = setdiff(rownames(nes1),rownames(nes3_ges))

#nes3t_del = nes3t[,-which(colnames(nes3t)==diffcn)]

# Create meta_data
meta1_ges = data.frame(cell_id = colnames(nes_G1),dataset='G1')
meta2_ges = data.frame(cell_id = colnames(nes_G2),dataset='G2')
meta3_ges = data.frame(cell_id = colnames(nes_Gsn),dataset='G_sn')
meta1_lac = data.frame(cell_id = colnames(nes_L1),dataset='L1')
meta2_lac = data.frame(cell_id = colnames(nes_L2),dataset='L2')
meta3_lac = data.frame(cell_id = colnames(nes_Lsn),dataset='L_sn')

# Merge
meta_data = rbind(meta1_ges,meta2_ges,meta3_ges,meta1_lac,meta2_lac,meta3_lac)
nes_merge = cbind(nes_G1,nes_G2,nes_Gsn,nes_L1,nes_L2,nes_Lsn)
save(nes_merge,meta_data,file='data/BPA_GO_merge_mouse.rda')

meta_data = rbind(meta1_ges,meta2_ges,meta1_lac,meta2_lac)
nes_merge = cbind(nes_G1,nes_G2,nes_L1,nes_L2)
save(nes_merge,meta_data,file='data/BPA_GO_merge_mouse_nosn.rda')

meta_data = rbind(meta1_ges,meta2_ges,meta3_ges)
nes_merge = cbind(nes_G1,nes_G2,nes_Gsn)
save(nes_merge,meta_data,file='data/gestation/mouse/BPA_GO_merge_mouse_ges.rda')

meta_data = rbind(meta1_lac,meta2_lac,meta3_lac)
nes_merge = cbind(nes_L1,nes_L2,nes_Lsn)
save(nes_merge,meta_data,file='data/lactation/mouse/BPA_GO_merge_mouse_lac.rda')


load("data/lactation/mouse/BPA_CP_sc1_lac.rda")
load("data/lactation/mouse/BPA_CP_sc2_lac.rda")
load("data/lactation/mouse/BPA_CP_sn1_lac.rda")

load("data/gestation/mouse/BPA_CP_sc1_ges.rda")
load("data/gestation/mouse/BPA_CP_sc2_ges.rda")
load("data/gestation/mouse/BPA_CP_sn1_ges.rda")

uni = union(rownames(nes1),rownames(nes3_ges))
diffcn = setdiff(rownames(nes1),rownames(nes3_ges))

#nes3t_del = nes3t[,-which(colnames(nes3t)==diffcn)]

rownames(nes3) == rownames(nes1)

# Create meta_data
meta1_ges = data.frame(cell_id = colnames(nes1_ges),dataset='G1')
meta2_ges = data.frame(cell_id = colnames(nes2_ges),dataset='G2')
meta3_ges = data.frame(cell_id = colnames(nes3_ges),dataset='G_sn')
meta1_lac = data.frame(cell_id = colnames(nes1),dataset='L1')
meta2_lac = data.frame(cell_id = colnames(nes2),dataset='L2')
meta3_lac = data.frame(cell_id = colnames(nes3),dataset='L_sn')

# Merge
meta_data = rbind(meta1_ges,meta2_ges,meta3_ges,meta1_lac,meta2_lac,meta3_lac)
nes_merge = cbind(nes1_ges,nes2_ges,nes3_ges,nes1,nes2,nes3)
save(nes_merge,meta_data,file='data/BPA_CP_merge_mouse.rda')

meta_data = rbind(meta1_ges,meta2_ges,meta1_lac,meta2_lac)
nes_merge = cbind(nes1_ges,nes2_ges,nes1,nes2)
save(nes_merge,meta_data,file='data/BPA_CP_merge_mouse_nosn.rda')



nes1t = t(nes1)
nes2t = t(nes2)
nes3t = t(nes3)

diffcn = setdiff(colnames(nes3t),colnames(nes1t))

#nes3t_del = nes3t[,-which(colnames(nes3t)==diffcn)]

colnames(nes3t) == colnames(nes1t)

# Create meta_data
meta1 = data.frame(cell_id = rownames(nes1t),dataset='sc1')
meta2 = data.frame(cell_id = rownames(nes2t),dataset='sc2')
meta3 = data.frame(cell_id = rownames(nes3t),dataset='sn1')
meta_data = rbind(meta1,meta2,meta3)

# Merge
nes_merge = rbind(nes1t,nes2t,nes3t)
nes_merge_t = t(nes_merge)
save(nes_merge_t,meta_data,file='data/lactation/mouse/BPA_GO_merge_lac.rda')



