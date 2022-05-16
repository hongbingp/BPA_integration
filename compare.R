annotation <- read.csv('data/GSE106273_ClusterAnnotation.csv')
annotation <- annotation[which(annotation$Condition %in% c('G','L')), ]
annotation$SubCluster <- gsub('10.1','Fbb',annotation$SubCluster)
annotation$SubCluster <- gsub('10.3','Endo',annotation$SubCluster)
annotation$SubCluster <- gsub('6-1','Imn',annotation$SubCluster)
annotation$SubCluster <- gsub('6-2','Imn',annotation$SubCluster)

annotation$SuperCluster[annotation$SubClusterNumbers %in% c('C16','C17')] <- 'Imn'
annotation$SuperCluster[annotation$SubClusterNumbers == 'C18'] <- 'Fbb'
annotation$SuperCluster[annotation$SubClusterNumbers == 'C19'] <- 'Endo'

write.csv(annotation,'data/annotation_gl.csv')




annotation <- read.csv('data/annotation_gl.csv')
meta = seurat@meta.data
meta = separate(meta,cell.id,c('barcode1','barcode2'),'-')

data = left_join(meta,annotation,by=c('barcode1'='bcs','dataset'='SampleID'))
#sc = data2[which(is.na(data2$SubCluster) & data2$dataset %in% c('sc1','sc2')),]

meta_data = data


############
meta_data[is.na(meta_data$SubCluster),]$SubCluster = 'NA'
meta_data[is.na(meta_data$SuperCluster),]$SuperCluster = 'NA'


fal = data[which(data$keep==FALSE),]


sn <- read.table('data/p13mg.meta.txt')
sn$bcs = gsub('P13_','',rownames(sn))

umap$cell.id = rownames(umap)

umap = left_join(umap,sn,by=c('cell.id'='bcs'))

