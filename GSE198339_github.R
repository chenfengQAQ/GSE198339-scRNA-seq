library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(monocle)
library(CellChat)
library(circlize)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(clustree)
library(data.table)
library(multtest)
library(metap)
library(wordcloud)
library(ggpubr)
library(ggsci)
library(future)
library(data.table)

options(scipen = 20000,future.globals.maxSize=1000*1024^2)

setwd('D:/GEO_test/GSE198339/')

#cell clustering --- control1
control1_h5<-Read10X_h5('GSM5945252_Participant_1_raw_gene_bc_matrices_h5.h5')
control1<-CreateSeuratObject(counts = control1_h5,project = 'control1',
                             min.cells = 3,min.features = 200)
control1[['percent.mt']]<-PercentageFeatureSet(control1,pattern = '^MT-')
VlnPlot(control1,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
control1<-subset(control1,subset=nFeature_RNA>200 & nFeature_RNA<2000&
                   percent.mt<25)
control1<-NormalizeData(control1)
control1<-FindVariableFeatures(control1)
control1<-ScaleData(control1)
control1<-RunPCA(control1)

saveRDS(control1,'D:/GEO_test/GSE198339/control1/control1_runpca.rds')

DimPlot(control1,reduction = 'pca',pt.size = 1)
DimHeatmap(control1,cells = 500,dims = 1:30,balanced = TRUE)
ElbowPlot(control1,ndims = 50)


#The cumulative contribution of the principal component is greater than 90%
#The PC contributes less than 5% to the variance
#The difference between two consecutive PCs is less than 0.1%
pct <- control1[["pca"]]@stdev / sum(control1[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1],
           sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
                decreasing = T)[1] + 1)
ElbowPlot(control1,ndims = 50)$data %>% ggplot()+
  geom_point(aes(x=dims,y=stdev))+
  geom_vline(xintercept = pcs,color = 'darkred')+
  theme_bw()+labs(title = 'Elbow plot : quantitative approach')

plot_df<-data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df,aes(cumu,pct,label = rank,color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

#different resolutions
par(mfrow=c(1,1))
ptm = Sys.time()
for(i in seq(0.1,1,by=0.1)){
  control1<-FindNeighbors(control1,dims = 1:12)
  control1<-FindClusters(control1,resolution = i)
  control1<-RunUMAP(control1,dims = 1:12)
  picture1<-DimPlot(control1,reduction = 'umap',label = TRUE,pt.size = 1.5)
  ggsave(filename = paste0('dim12_res',i,'.png'),plot = picture1,
         path = 'D:/GEO_test/GSE198339/control1/diff_resolution_control1/')
}
Sys.time() - ptm
saveRDS(control1,'D:/GEO_test/GSE198339/control1/control1_dim12_re0.1-1.rds')

#Clusters phylogenetic tree
clustree(control1@meta.data,prefix = 'RNA_snn_res.')

#Cell cycle
control1<-CellCycleScoring(object = control1,s.features = cc.genes$s.genes,
                           g2m.features = cc.genes$g2m.genes)
DimPlot(control1,reduction = 'umap',group.by = 'Phase',label = TRUE,
        pt.size = 1.5)

#Determined resolution
control1<-readRDS('D:/GEO_test/GSE198339/control1/control1_runpca.rds')
control1<-FindNeighbors(control1,dims = 1:12)
control1<-FindClusters(control1,resolution = 0.5)
control1<-RunUMAP(control1,dims = 1:12)
DimPlot(control1,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(control1,'D:/GEO_test/GSE198339/control1/control1_dim11_re0.5_unnamed.rds')

#Different visual plots of T marker genes
VlnPlot(control1,features = c('CD3D','CD3G','CD3E'),pt.size = 0)
FeaturePlot(control1,features = c('CD3D','CD3G','CD3E'))
RidgePlot(control1,features = c('CD3D','CD3G','CD3E'))
DotPlot(control1,features = c('CD3D','CD3G','CD3E'))

#Differential genes between cell clusters
demarkers<-FindAllMarkers(control1,only.pos = TRUE)
fwrite(demarkers,'D:/GEO_test/GSE198339/control1/control1_clusters_demarkers.csv',
       row.names = TRUE)

#Cell clusters annotation
DotPlot(control1,features = c('CD3D','CD3G','CD3E',   #0  3  5
                              'MS4A1','CD19','CD79A',   #1  7
                              'NKG7','GNLY',   #4
                              'CST3','LYZ',  #6
                              'CD68','CD163','CD14'
))+RotatedAxis()

new.cluster.ids <- c('T','B','Unknown','T','NK','T','DC','B')
names(new.cluster.ids) <- levels(control1)
control1 <- RenameIdents(control1, new.cluster.ids)
DimPlot(control1, reduction = "umap", label = TRUE, pt.size = 1.5,repel = TRUE)
saveRDS(control1,'D:/GEO_test/GSE198339/control1/control1_dim11_re0.5_named.rds')


#cell clustering --- AS1
setwd('D:/GEO_test/GSE198339/')
AS1_h5<-Read10X_h5('GSM5945256_Participant_5_raw_gene_bc_matrices_h5.h5')
AS1<-CreateSeuratObject(counts = AS1_h5,project = 'AS1',
                        min.cells = 3,min.features = 200)
AS1[['percent.mt']]<-PercentageFeatureSet(AS1,pattern = '^MT-')
VlnPlot(AS1,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
AS1<-subset(AS1,subset=nFeature_RNA>200 & nFeature_RNA<2000&
              percent.mt<25)
AS1<-NormalizeData(AS1)
AS1<-FindVariableFeatures(AS1)
AS1<-ScaleData(AS1)
AS1<-RunPCA(AS1)
# saveRDS(AS1,'D:/GEO_test/GSE198339/AS1/AS1_runpca.rds')
AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_runpca.rds')
DimPlot(AS1,reduction = 'pca',pt.size = 1)
DimHeatmap(AS1,cells = 500,dims = 1:30,balanced = TRUE)
ElbowPlot(AS1,ndims = 50)

#The cumulative contribution of the principal component is greater than 90%
#The PC contributes less than 5% to the variance
#The difference between two consecutive PCs is less than 0.1%
pct <- AS1[["pca"]]@stdev / sum(AS1[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1],
           sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
                decreasing = T)[1] + 1)
ElbowPlot(AS1,ndims = 50)$data %>% ggplot()+
  geom_point(aes(x=dims,y=stdev))+
  geom_vline(xintercept = pcs,color = 'darkred')+
  theme_bw()+labs(title = 'Elbow plot : quantitative approach')

plot_df<-data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df,aes(cumu,pct,label = rank,color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

#different resolutions
par(mfrow=c(1,1))
ptm = Sys.time()
for(i in seq(0.1,1,by=0.1)){
  AS1<-FindNeighbors(AS1,dims = 1:12)
  AS1<-FindClusters(AS1,resolution = i)
  AS1<-RunUMAP(AS1,dims = 1:12)
  picture1<-DimPlot(AS1,reduction = 'umap',label = TRUE,pt.size = 1.5)
  ggsave(filename = paste0('dim12_res',i,'.png'),plot = picture1,
         path = 'D:/GEO_test/GSE198339/AS1/diff_resolution_AS1/')
}
Sys.time() - ptm
saveRDS(AS1,'D:/GEO_test/GSE198339/AS1/AS1_dim12_re0.1-1.rds')

#Clusters phylogenetic tree
clustree(AS1@meta.data,prefix = 'RNA_snn_res.')

#Cell cycle
AS1<-CellCycleScoring(object = AS1,s.features = cc.genes$s.genes,
                      g2m.features = cc.genes$g2m.genes)
DimPlot(AS1,reduction = 'umap',group.by = 'Phase',label = TRUE,
        pt.size = 1.5)

#Determined resolution
AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_runpca.rds')
AS1<-FindNeighbors(AS1,dims = 1:12)
AS1<-FindClusters(AS1,resolution = 0.5)
AS1<-RunUMAP(AS1,dims = 1:12)
DimPlot(AS1,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(AS1,'D:/GEO_test/GSE198339/AS1/AS1_dim12_re0.5_unnamed.rds')

#Different visual plots of T marker genes
VlnPlot(AS1,features = c('CD3D','CD3G','CD3E'),pt.size = 0)
FeaturePlot(AS1,features = c('CD3D','CD3G','CD3E'))
RidgePlot(AS1,features = c('CD3D','CD3G','CD3E'))
DotPlot(AS1,features = c('CD3D','CD3G','CD3E'))

#Differential genes between cell clusters
demarkers<-FindAllMarkers(AS1,only.pos = TRUE)
fwrite(demarkers,'D:/GEO_test/GSE198339/AS1/AS1_clusters_demarkers.csv',
       row.names = TRUE)

#Cell clusters annotation
DotPlot(AS1,features = c('CD3D','CD3G','CD3E',   #0  3  5
                         'MS4A1','CD19','CD79A',   #1  7
                         'NKG7','GNLY',   #4
                         'CST3','LYZ',  #6
                         'CD68','CD163','CD14'
))+RotatedAxis()

new.cluster.ids <- c('T','B','Unknown','T','NK','T','DC','B')
names(new.cluster.ids) <- levels(AS1)
AS1 <- RenameIdents(AS1, new.cluster.ids)
DimPlot(AS1, reduction = "umap", label = TRUE, pt.size = 1.5,repel = TRUE)
saveRDS(AS1,'D:/GEO_test/GSE198339/AS1/AS1_dim12_re0.5_named.rds')


#DoubletFinder --- control1
control1<-readRDS('D:/GEO_test/GSE198339/control1/control1_dim11_re0.5_unnamed.rds')
ptm = Sys.time()
sweep.res.list_control1 <- paramSweep(control1, PCs = 1:11, sct = FALSE)
sweep.stats_control1 <- summarizeSweep(sweep.res.list_control1, GT = FALSE)
bcmvn_control1 <- find.pK(sweep.stats_control1)
best_pk <- as.numeric(as.vector(bcmvn_control1$pK[which.max(bcmvn_control1$BCmetric)]))
best_pk
homotypic.prop <- modelHomotypic(control1$seurat_clusters)
homotypic.prop
DoubletRate<-0.008*(nrow(control1@meta.data)/1000)
nExp_poi <- round(DoubletRate * nrow(control1@meta.data))
nExp_poi
nExp_poi_adj<-round(nExp_poi*(1-homotypic.prop))
nExp_poi_adj
control1<-doubletFinder(control1,PCs = 1:11,pN = 0.25,pK = best_pk,nExp = nExp_poi,
                        reuse.pANN = FALSE,sct = FALSE)
control1<-doubletFinder(control1,PCs = 1:11,pN = 0.25,pK = best_pk,nExp = nExp_poi_adj,
                        reuse.pANN = FALSE,sct = FALSE)
DimPlot(control1,reduction = 'umap',group.by = 'DF.classifications_0.25_0.18_24')
Sys.time() - ptm
control1<-subset(control1,subset = DF.classifications_0.25_0.18_24 == 'Singlet')
saveRDS(control1,'D:/GEO_test/GSE198339/control1/control1_dim11_res0.5_unnamed_Singlet.rds')


#DoubletFinder --- AS1
AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_dim12_re0.5_unnamed.rds')
ptm = Sys.time()
sweep.res.list_AS1 <- paramSweep(AS1, PCs = 1:12, sct = FALSE)
sweep.stats_AS1 <- summarizeSweep(sweep.res.list_AS1, GT = FALSE)
bcmvn_AS1 <- find.pK(sweep.stats_AS1)
best_pk <- as.numeric(as.vector(bcmvn_AS1$pK[which.max(bcmvn_AS1$BCmetric)]))
best_pk
homotypic.prop <- modelHomotypic(AS1$seurat_clusters)
homotypic.prop
DoubletRate<-0.008*(nrow(AS1@meta.data)/1000)
nExp_poi <- round(DoubletRate * nrow(AS1@meta.data))
nExp_poi
nExp_poi_adj<-round(nExp_poi*(1-homotypic.prop))
nExp_poi_adj
AS1<-doubletFinder(AS1,PCs = 1:12,pN = 0.25,pK = best_pk,nExp = nExp_poi,
                   reuse.pANN = FALSE,sct = FALSE)
AS1<-doubletFinder(AS1,PCs = 1:12,pN = 0.25,pK = best_pk,nExp = nExp_poi_adj,
                   reuse.pANN = FALSE,sct = FALSE)
DimPlot(AS1,reduction = 'umap',group.by = 'DF.classifications_0.25_0.05_16')
Sys.time() - ptm
AS1<-subset(AS1,subset = DF.classifications_0.25_0.05_16 == 'Singlet')
saveRDS(AS1,'D:/GEO_test/GSE198339/AS1/AS1_dim12_res0.5_unnamed_Singlet.rds')


#merge_data --- control1 and AS1
control1<-readRDS('D:/GEO_test/GSE198339/control1/control1_dim11_res0.5_unnamed_Singlet.rds')
AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_dim12_res0.5_unnamed_Singlet.rds')

control1_again<-CreateSeuratObject(counts = control1@assays$RNA,
                                   meta.data = control1@meta.data[,1:3],
                                   project = 'control1')
AS1_again<-CreateSeuratObject(counts = AS1@assays$RNA,
                              meta.data = AS1@meta.data[,1:3],
                              project = 'AS1')

#cell clustering 
data.merged<-merge(control1_again,AS1_again)
data.merged[['percent.mt']]<-PercentageFeatureSet(data.merged,pattern = '^MT-')
VlnPlot(data.merged,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
data.merged<-NormalizeData(data.merged)
data.merged<-FindVariableFeatures(data.merged)
data.merged<-ScaleData(data.merged)
data.merged<-RunPCA(data.merged)

saveRDS(data.merged,'D:/GEO_test/GSE198339/control1_AS1_merge/data.merged_runpca.rds')

DimPlot(data.merged,reduction = 'pca',pt.size = 1)
DimHeatmap(data.merged,cells = 500,dims = 1:30,balanced = TRUE)
ElbowPlot(data.merged,ndims = 50)

#The cumulative contribution of the principal component is greater than 90%
#The PC contributes less than 5% to the variance
#The difference between two consecutive PCs is less than 0.1%
pct <- data.merged[["pca"]]@stdev / sum(data.merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1],
           sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
                decreasing = T)[1] + 1)
ElbowPlot(data.merged,ndims = 50)$data %>% ggplot()+
  geom_point(aes(x=dims,y=stdev))+
  geom_vline(xintercept = pcs,color = 'darkred')+
  theme_bw()+labs(title = 'Elbow plot : quantitative approach')

plot_df<-data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df,aes(cumu,pct,label = rank,color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

#different resolutions
par(mfrow=c(1,1))
ptm = Sys.time()
for(i in seq(0.1,1,by=0.1)){
  data.merged<-FindNeighbors(data.merged,dims = 1:20)
  data.merged<-FindClusters(data.merged,resolution = i)
  data.merged<-RunUMAP(data.merged,dims = 1:20)
  picture1<-DimPlot(data.merged,reduction = 'umap',label = TRUE,pt.size = 1.5)
  ggsave(filename = paste0('dim20_res',i,'.png'),plot = picture1,
         path = 'D:/GEO_test/GSE198339/control1_AS1_merge/diff_resolution_merge/')
}
Sys.time() - ptm
saveRDS(data.merged,'D:/GEO_test/GSE198339/control1_AS1_merge/data.merge_dim20_res0.1-1.rds')

#Clusters phylogenetic tree
clustree(data.merged@meta.data,prefix = 'RNA_snn_res.')

#Determined resolution
data.merged<-readRDS('D:/GEO_test/GSE198339/control1_AS1_merge/data.merged_runpca.rds')
data.merged<-FindNeighbors(data.merged,dims = 1:20)
data.merged<-FindClusters(data.merged,resolution = 0.2)
data.merged<-RunUMAP(data.merged,dims = 1:20)
DimPlot(data.merged,reduction = 'umap',pt.size = 1.5,label = TRUE,repel = TRUE,
        group.by = 'orig.ident')
saveRDS(data.merged,'D:/GEO_test/GSE198339/control1_AS1_merge/data.merge_dim20_res0.2.rds')

#Differential genes between cell clusters from control1 and AS1
orig.ident_df1<-data.frame(table(data.merged@meta.data$orig.ident))
data.merged[['status']]<-rep(c('control','AS'),c(orig.ident_df1$Freq[2],
                                                 orig.ident_df1$Freq[1]))
data.merged[['RNA']]<-JoinLayers(data.merged[['RNA']])

Idents(data.merged)<-'seurat_clusters'
DefaultAssay(data.merged)<-'RNA'

data_middle_2<-data.frame()
ptm = Sys.time()
for(i in levels(data.merged)){
  data.merged.markers<-FindConservedMarkers(data.merged,ident.1 = i,
                                            grouping.var = 'status',
                                            verbose = FALSE,only.pos=TRUE)
  data.merged.markers<-data.frame(gene_name = row.names(data.merged.markers),
                                  data.merged.markers)
  row.names(data.merged.markers)<-NULL
  seurat_clusters<-data.frame(seurat_clusters = rep(i,length(row.names(data.merged.markers))))
  data_middle_1<-cbind.data.frame(data.merged.markers,seurat_clusters)
  data_middle_2<-rbind.data.frame(data_middle_2,data_middle_1)
}
Sys.time() - ptm
write.csv(data_middle_2,'D:/GEO_test/GSE198339/control1_AS1_merge/data.merged.markers.csv')

DimPlot(data.merged,reduction = 'umap',pt.size = 1.5,label = TRUE,repel = TRUE)

#Cell clusters annotation
DotPlot(data.merged,features = c('CD3D','CD3G','CD3E',   #0  3
                                 'MS4A1','CD19','CD79A',  #1  6
                                 'NKG7','GNLY',  #4
                                 'CST3','LYZ',  #5
                                 'CD68','CD163','CD14'
))+RotatedAxis()

new.cluster.ids <- c('T','B','Unknown','T','NK','DC','B')
names(new.cluster.ids) <- levels(data.merged)
data.merged <- RenameIdents(data.merged, new.cluster.ids)
DimPlot(data.merged, reduction = "umap", label = TRUE, pt.size = 1.5,repel = TRUE)
saveRDS(data.merged,'D:/GEO_test/GSE198339/control1_AS1_merge/data.merge_dim20_res0.2_named.rds')

#Cell proportion
cell_proportion<-as.data.frame(prop.table(
  table(Idents(data.merged),data.merged@meta.data$status)
  ,margin = 2))
names(cell_proportion)<-c('cluster','status','proportion')

ggplot(cell_proportion,aes(x=status,y=proportion,fill=cluster))+
  geom_bar(position = 'stack',stat="identity")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(color="black"))+
  scale_fill_manual(values = rainbow(25))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))+
  geom_text(aes(label=proportion),position = position_stack(vjust = 0.5),
            color='black')

#Differential cell clusters genes between control1 and AS1
data.merged@meta.data$celltype<-Idents(data.merged)
data.merged@meta.data$celltype.status<-paste(data.merged@meta.data$celltype,
                                             data.merged@meta.data$status,
                                             sep = '_')

Idents(data.merged)<-'celltype.status'
demarkers3<-FindMarkers(data.merged,ident.1 = 'T_AS',
                        ident.2 = 'T_control',verbose = FALSE)

write.csv(demarkers,'D:/GEO_test/GSE198339/control1_AS1_merge/T demarkers_group.csv')


#Subcluster analysis --- T
#sub_cell clustering 
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/T')
Idents(data.merged)<-'celltype'
sub_cell<-subset(data.merged,idents=c('T'))
sub_cell<-NormalizeData(sub_cell)
sub_cell<-FindVariableFeatures(sub_cell)
sub_cell<-ScaleData(sub_cell)
sub_cell<-RunPCA(sub_cell)

saveRDS(sub_cell,"D:/GEO_test/GSE198339/control1_AS1_merge/T/T_runpca.rds")

DimPlot(sub_cell,reduction = 'pca',pt.size = 1)
DimHeatmap(sub_cell,cells = 500,dims = 1:30,balanced = TRUE)
ElbowPlot(sub_cell,ndims = 50)

#The cumulative contribution of the principal component is greater than 90%
#The PC contributes less than 5% to the variance
#The difference between two consecutive PCs is less than 0.1%
pct <- sub_cell[["pca"]]@stdev / sum(sub_cell[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1],
           sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
                decreasing = T)[1] + 1)
ElbowPlot(sub_cell,ndims = 50)$data %>% ggplot()+
  geom_point(aes(x=dims,y=stdev))+
  geom_vline(xintercept = pcs,color = 'darkred')+
  theme_bw()+labs(title = 'Elbow plot : quantitative approach')
plot_df<-data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df,aes(cumu,pct,label = rank,color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

#different resolutions
par(mfrow=c(1,1))
ptm = Sys.time()
for(i in seq(0.1,1,by=0.1)){
  sub_cell<-FindNeighbors(sub_cell,dims = 1:10)
  sub_cell<-FindClusters(sub_cell,resolution = i)
  sub_cell<-RunUMAP(sub_cell,dims = 1:10)
  picture1<-DimPlot(sub_cell,reduction = 'umap',label = TRUE,pt.size = 1.5)
  ggsave(filename = paste0('dim10_res',i,'.png'),plot = picture1,
         path = 'D:/GEO_test/GSE198339/control1_AS1_merge/T/diff_resolution_T')
}
Sys.time() - ptm
saveRDS(sub_cell,"D:/GEO_test/GSE198339/control1_AS1_merge/T/T_dim10_res0.1-1.rds")

#Clusters phylogenetic tree
clustree(sub_cell@meta.data,prefix = 'RNA_snn_res.')

#Cell cycle
sub_cell<-CellCycleScoring(object = sub_cell,s.features = cc.genes$s.genes,
                           g2m.features = cc.genes$g2m.genes)
DimPlot(sub_cell,reduction = 'umap',group.by = 'Phase',label = TRUE,
        pt.size = 1.5)

#Determined resolution
sub_cell<-readRDS("D:/GEO_test/GSE198339/control1_AS1_merge/T/T_runpca.rds")
sub_cell<-FindNeighbors(sub_cell,dims = 1:10)
sub_cell<-FindClusters(sub_cell,resolution = 0.1)
sub_cell<-RunUMAP(sub_cell,dims = 1:10)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(sub_cell,"D:/GEO_test/GSE198339/control1_AS1_merge/T/T_dim10_res0.1_unnamed.rds")

#Differential genes between sub_cell clusters from control1 and AS1
Idents(sub_cell)<-'seurat_clusters'
DefaultAssay(sub_cell)<-'RNA'

data_middle_4<-data.frame()
ptm = Sys.time()
for(i in levels(sub_cell)){
  sub_cell.markers<-FindConservedMarkers(sub_cell,ident.1 = i,
                                         grouping.var = 'status',
                                         verbose = FALSE,only.pos=TRUE)
  
  sub_cell.markers<-data.frame(gene_name = row.names(sub_cell.markers),
                               sub_cell.markers)
  row.names(sub_cell.markers)<-NULL
  seurat_clusters<-data.frame(seurat_clusters = rep(i,length(row.names(sub_cell.markers))))
  data_middle_3<-cbind.data.frame(sub_cell.markers,seurat_clusters)
  data_middle_4<-rbind.data.frame(data_middle_4,data_middle_3)
}
Sys.time() - ptm

#Exclude ribosomal genes
markers.noRio<-data_middle_4[!grepl('^RP[SL]',data_middle_4$gene_name,
                                    ignore.case = FALSE),]

write.csv(markers.noRio,'D:/GEO_test/GSE198339/control1_AS1_merge/T/sub_cell.markers.noRio.csv')

#Sub_cell clusters annotation
DotPlot(sub_cell,features = c('CCR7','SELL',  #0
                              'CCL5','GZMA','GZMK','CD8A' #1
))+RotatedAxis()

new.cluster.ids<-c('navie CD4+T','effector memory CD8+T')
names(new.cluster.ids) <- levels(sub_cell)
sub_cell <- RenameIdents(sub_cell, new.cluster.ids)
sub_cell@meta.data$sub_celltype<-Idents(sub_cell)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE,
        split.by = 'orig.ident')
saveRDS(sub_cell,"D:/GEO_test/GSE198339/control1_AS1_merge/T/T_dim10_res0.1_named.rds")

#Sub_cell proportion
Idents(sub_cell)<-'sub_celltype'
cell_proportion<-as.data.frame(prop.table(
  table(Idents(sub_cell),sub_cell@meta.data$status)
  ,margin = 2))
names(cell_proportion)<-c('cluster','status','proportion')

ggplot(cell_proportion,aes(x=status,y=proportion,fill=cluster))+
  geom_bar(position = 'stack',stat="identity")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(color="black"))+
  scale_fill_manual(values = rainbow(25))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))+
  geom_text(aes(label=proportion),position = position_stack(vjust = 0.5),
            color='black')

#Differential sub_cell clusters genes between control1 and AS1
sub_cell@meta.data$sub_celltype.status<-paste(sub_cell@meta.data$sub_celltype,
                                              sub_cell@meta.data$status,
                                              sep = '_')
Idents(sub_cell)<-'sub_celltype.status'
demarkers4<-FindMarkers(sub_cell,ident.1 = 'navie CD4+T_AS',
                        ident.2 = 'navie CD4+T_control',verbose = FALSE)

demarkers4.noRio<-demarkers4[!grepl('^RP[SL]',row.names(demarkers4),
                                    ignore.case = FALSE),]
write.csv(demarkers4.noRio,'D:/GEO_test/GSE198339/control1_AS1_merge/T/navie CD4+T.markers.noRio.csv')


#cellchat --- control1 and AS1
#load conda environment and python packages
Sys.setenv(PATH= paste("D:/Miniconda3-2024/envs/py39/Library/bin",Sys.getenv()["PATH"],sep=";"))
Sys.setenv(RETICULATE_PYTHON = "D:/Miniconda3-2024/envs/py39/python.exe")
library(reticulate)
use_condaenv("py39")
py_config()
import("numpy")
import("umap")
np<-import("numpy")
print(np$version$full_version)

data.merged<-readRDS('D:/GEO_test/GSE198339/control1_AS1_merge/data.merge_dim20_res0.2_named.rds')
sub_cell_control1<-subset(data.merged,subset = orig.ident =='control1')
sub_cell_AS1<-subset(data.merged,subset = orig.ident =='AS1')

#cellchat --- control1
cellchat_control1<-createCellChat(object = sub_cell_control1,group.by = 'ident',
                                  assay = 'RNA')
cellchatdb<-CellChatDB.human
cellchatdb.use<-subsetDB(cellchatdb,search = 'Secreted Signaling',
                         key = 'annotation')
cellchat_control1@DB<-cellchatdb.use
cellchat_control1<-subsetData(cellchat_control1)

#Parallel computing
future::plan("multisession", workers = 4)

cellchat_control1<-identifyOverExpressedGenes(cellchat_control1)
cellchat_control1<-identifyOverExpressedInteractions(cellchat_control1)
cellchat_control1<-projectData(cellchat_control1,PPI.human)
cellchat_control1<-computeCommunProb(cellchat_control1,type = 'triMean')
cellchat_control1<-filterCommunication(cellchat_control1,min.cells = 10)
cellchat_control1<-computeCommunProbPathway(cellchat_control1)
cellchat_control1<-aggregateNet(cellchat_control1)

saveRDS(cellchat_control1,'D:/GEO_test/GSE198339/control1_AS1_merge/cellchat/cellchat_control1.rds')

cellchat_control1<-netAnalysis_computeCentrality(cellchat_control1)

#cellchat --- AS1
cellchat_AS1<-createCellChat(object = sub_cell_AS1,group.by = 'ident',
                             assay = 'RNA')
cellchat_AS1@DB<-cellchatdb.use
cellchat_AS1<-subsetData(cellchat_AS1)
cellchat_AS1<-identifyOverExpressedGenes(cellchat_AS1)
cellchat_AS1<-identifyOverExpressedInteractions(cellchat_AS1)
cellchat_AS1<-projectData(cellchat_AS1,PPI.human)
cellchat_AS1<-computeCommunProb(cellchat_AS1,type = 'triMean')
cellchat_AS1<-filterCommunication(cellchat_AS1,min.cells = 10)
cellchat_AS1<-computeCommunProbPathway(cellchat_AS1)
cellchat_AS1<-aggregateNet(cellchat_AS1)

saveRDS(cellchat_AS1,'D:/GEO_test/GSE198339/control1_AS1_merge/cellchat/cellchat_AS1.rds')

cellchat_AS1<-netAnalysis_computeCentrality(cellchat_AS1)

#merge 
object.list<-list(control1 = cellchat_control1,AS1 = cellchat_AS1)
saveRDS(object.list,'D:/GEO_test/GSE198339/control1_AS1_merge/cellchat/object.list.rds')

cellchat_merge<-mergeCellChat(object.list,add.names = names(object.list))
saveRDS(cellchat_merge,'D:/GEO_test/GSE198339/control1_AS1_merge/cellchat/cellchat_merge.rds')

#cellchat visual
#Infer the number and intensity of signal effects of two groups from the global perspective
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/cellchat/')
gg1<-compareInteractions(cellchat_merge,group = c(1,2),show.legend = FALSE)
gg2<-compareInteractions(cellchat_merge,group = c(1,2),show.legend = FALSE,
                         measure = 'weight')
gg1+gg2


#Network diagram of signal numbers and strength
par(mfrow = c(1,2),xpd=TRUE)
netVisual_diffInteraction(cellchat_merge,weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merge,weight.scale = TRUE,measure = 'weight')

#Heatmap of signal numbers and strength
gg1<-netVisual_heatmap(cellchat_merge)
gg2<-netVisual_heatmap(cellchat_merge,measure = 'weight')
gg1+gg2

#Network diagram of signal numbers between control1 and AS1
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", 
                                       names(object.list)[i]))
}


#Changes in incoming and outgoing signaling pattern of different groups of cell groups
num.link <- sapply(object.list, 
                   function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


#Specific signal changes in B cells  between control1 and AS1
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "T", 
                                            signaling.exclude = "CXCL")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "B",
                                            signaling.exclude = c("CXCL"))
patchwork::wrap_plots(plots = list(gg1,gg2))


#Functional similarity of signaling pathways between the two groups
cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge,type = "functional")
cellchat_merge <- netEmbedding(cellchat_merge,type = "functional")
cellchat_merge <- netClustering(cellchat_merge,type = "functional")
netVisual_embeddingPairwise(cellchat_merge,type = "functional",label.size = 3.5)

#Structural similarity of signaling pathways between the two groups
cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge,type = "structural")
cellchat_merge <- netEmbedding(cellchat_merge,type = "structural")
cellchat_merge <- netClustering(cellchat_merge,type = "structural")
netVisual_embeddingPairwise(cellchat_merge,type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat_merge,type = "structural", nCol = 2)

#Euclidean distance of signal pathways between two groups
rankSimilarity(cellchat_merge, type = "functional")
rankSimilarity(cellchat_merge, type = "structural")

#Conserved signaling pathways in both groups
gg1 <- rankNet(cellchat_merge, mode = "comparison",measure = "weight", 
               sources.use = NULL,targets.use = NULL,stacked = T,do.stat = TRUE)
gg2 <- rankNet(cellchat_merge, mode = "comparison",measure = "weight", 
               sources.use = NULL,targets.use = NULL,stacked = F,do.stat = TRUE)
gg1 + gg2


#Heatmap of outgoing signaling pattern between two groups
i = 1
pathway.union <- union(object.list[[1]]@netP$pathways, 
                       object.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]],pattern = "outgoing",
                                        signaling = pathway.union,
                                        title = names(object.list)[1],
                                        width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]],
                                        pattern = "outgoing",
                                        signaling = pathway.union,
                                        title = names(object.list)[2],
                                        width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Heatmap of incoming signaling pattern between two groups
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]],pattern = "incoming",
                                        signaling = pathway.union,
                                        title = names(object.list)[1],
                                        width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]],
                                        pattern = "incoming",
                                        signaling = pathway.union,
                                        title = names(object.list)[2],
                                        width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Heatmap of overall signaling pattern between two groups
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]],pattern = "all",
                                        signaling = pathway.union,
                                        title = names(object.list)[1],
                                        width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]],
                                        pattern = "all",
                                        signaling = pathway.union,
                                        title = names(object.list)[2],
                                        width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Probability of communication between ligands of T cell groups and receptors of other cell groups
# table(data.merged@meta.data$celltype)
# table(data.merged@meta.data$seurat_clusters)
netVisual_bubble(cellchat_merge,sources.use = 1,targets.use = c(2:5),
                 comparison = c(1, 2), angle.x = 45)

netVisual_bubble(cellchat_merge,sources.use = 2,targets.use = c(1,3:5),
                 comparison = c(1, 2), angle.x = 45)


#Ligands of T cell groups and receptors of other cell groups --- upregulate and downregulate
gg1 <- netVisual_bubble(cellchat_merge, sources.use = 1,targets.use = c(2:5),
                        comparison = c(1, 2),max.dataset = 2,
                        title.name = "Increased signaling in AS",
                        angle.x = 45,remove.isolate = FALSE)
gg2 <- netVisual_bubble(cellchat_merge,sources.use = 1,targets.use = c(2:5), 
                        comparison = c(1, 2),max.dataset = 1,
                        title.name = "Decreased signaling in AS",
                        angle.x = 45,remove.isolate = FALSE)
gg1 + gg2


#Upregulate and downregulate ligand-receptors by DEA
pos.dataset = "AS1"
features.name = paste0(pos.dataset, ".merged")

cellchat_merge <- identifyOverExpressedGenes(cellchat_merge,
                                             group.dataset = "datasets",
                                             pos.dataset = pos.dataset,
                                             features.name = features.name,
                                             only.pos = FALSE, thresh.pc = 0.1,
                                             thresh.fc = 0.05,thresh.p = 0.05,
                                             group.DE.combined = TRUE) 
net <- netMappingDEG(cellchat_merge, features.name = features.name, 
                     variable.all = TRUE)
net.up <- subsetCommunication(cellchat_merge,net = net,datasets = "AS1",
                              ligand.logFC = 0.05,receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat_merge,net = net,datasets = "control1",
                                ligand.logFC = -0.05,receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merge)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merge)

#Visualization of upregulate and downregulate signaling ligand-receptors
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_merge,pairLR.use = pairLR.use.up,sources.use = 1,
                        targets.use = c(2:11),comparison = c(1, 2),angle.x = 90,
                        remove.isolate = FALSE,
                        title.name = paste0("Up-regulated signaling in ",
                                            names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_merge, pairLR.use = pairLR.use.down,sources.use = 1,
                        targets.use = c(2:11),comparison = c(1, 2), angle.x = 90,
                        remove.isolate = FALSE,
                        title.name = paste0("Down-regulated signaling in ",
                                            names(object.list)[2]))
gg1 + gg2
gg2

#Visualization of upregulate and downregulate signaling ligand-receptors --- Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]],sources.use = 1,targets.use = c(2:11),
                     slot.name = 'net',net = net.up,lab.cex = 0.8,small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",
                                         names(object.list)[2]),show.legend = TRUE)
netVisual_chord_gene(object.list[[1]],sources.use = 1,targets.use = c(2:11),
                     slot.name = 'net',net = net.down,lab.cex = 0.8,small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ",
                                         names(object.list)[2]))

#Visualization of upregulate and downregulate signaling ligand-receptors --- Word cloud
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)


#Comparing signal changes inMIF signaling pathways
#Network diagram of signal strength
pathways.show <- c("MIF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), 
                           attribute = pathways.show) 
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "circle",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]),
                      label.edge = TRUE)
}

#Heatmap of signal pattern
pathways.show <- c("MIF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]],signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#Chord diagram of signal pattern
pathways.show <- c("MIF") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,
                      layout = "chord",
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}

#Expression levels of ligand-receptors genes from MIF
cellchat_merge@meta$datasets = factor(cellchat_merge@meta$datasets,
                                      levels = c("control1", "AS1"))
plotGeneExpression(cellchat_merge,signaling = "MIF",split.by = "datasets",
                   colors.ggplot = T,type = "violin")

saveRDS(cellchat_merge,'D:/GEO_test/GSE198339/control1_AS1_merge/cellchat/cellchat_merge_again.rds')


#monocle
#Singlet cell clusters annotation
control1<-readRDS('D:/GEO_test/GSE198339/control1/control1_dim11_res0.5_unnamed_Singlet.rds')

DotPlot(control1,features = c('CD3D','CD3G','CD3E',   #0,2,4
                              'MS4A1','CD19','CD79A',   #1 7
                              'NKG7','GNLY',   #5
                              'CST3','LYZ',  #6
                              'CD68','CD163','CD14' #8
))+RotatedAxis()
new.cluster.ids <- c('T','B','T','Unknown','T','NC','DC','B','Mono and Macro')
names(new.cluster.ids) <- levels(control1)
control1 <- RenameIdents(control1, new.cluster.ids)
control1@meta.data$celltype<-Idents(control1)
DimPlot(control1, reduction = "umap", label = TRUE, pt.size = 1.5,repel = TRUE)
saveRDS(control1,'D:/GEO_test/GSE198339/control1/control1_dim11_res0.5_named_Singlet.rds')

AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_dim12_res0.5_unnamed_Singlet.rds')
DotPlot(AS1,features = c('CD3D','CD3G','CD3E',   #0  3  5
                         'MS4A1','CD19','CD79A',   #1  7
                         'NKG7','GNLY',   #4
                         'CST3','LYZ',  #6
                         'CD68','CD163','CD14'
))+RotatedAxis()
new.cluster.ids <- c('T','B','Unknown','T','NK','T','DC','B')
names(new.cluster.ids) <- levels(AS1)
AS1 <- RenameIdents(AS1, new.cluster.ids)
AS1@meta.data$celltype<-Idents(AS1)
DimPlot(AS1, reduction = "umap", label = TRUE, pt.size = 1.5,repel = TRUE)
saveRDS(AS1,'D:/GEO_test/GSE198339/AS1/AS1_dim12_res0.5_named_Singlet.rds')

#T cell grouop monocle2 --- control1
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/monocle/control1_T/')
Idents(control1)<-'celltype'
sub_cell<-subset(control1,idents=c('T'))
sub_cell<-NormalizeData(sub_cell)
sub_cell<-FindVariableFeatures(sub_cell)
sub_cell<-ScaleData(sub_cell)
sub_cell<-RunPCA(sub_cell)

saveRDS(sub_cell,"T_runpca_control1.rds")

DimPlot(sub_cell,reduction = 'pca',pt.size = 1)
DimHeatmap(sub_cell,cells = 500,dims = 1:30,balanced = TRUE)
ElbowPlot(sub_cell,ndims = 50)

pct <- sub_cell[["pca"]]@stdev / sum(sub_cell[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1],
           sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
                decreasing = T)[1] + 1)
ElbowPlot(sub_cell,ndims = 50)$data %>% ggplot()+
  geom_point(aes(x=dims,y=stdev))+
  geom_vline(xintercept = pcs,color = 'darkred')+
  theme_bw()+labs(title = 'Elbow plot : quantitative approach')
plot_df<-data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df,aes(cumu,pct,label = rank,color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

par(mfrow=c(1,1))
ptm = Sys.time()
for(i in seq(0.1,1,by=0.1)){
  sub_cell<-FindNeighbors(sub_cell,dims = 1:18)
  sub_cell<-FindClusters(sub_cell,resolution = i)
  sub_cell<-RunUMAP(sub_cell,dims = 1:18)
  picture1<-DimPlot(sub_cell,reduction = 'umap',label = TRUE,pt.size = 1.5)
  ggsave(filename = paste0('dim18_res',i,'.png'),plot = picture1,
         path = 'D:/GEO_test/GSE198339/control1_AS1_merge/monocle/control1_T')
}
Sys.time() - ptm
saveRDS(sub_cell,"T_dim18_res0.1-1.rds")

clustree(sub_cell@meta.data,prefix = 'RNA_snn_res.')

sub_cell<-CellCycleScoring(object = sub_cell,s.features = cc.genes$s.genes,
                           g2m.features = cc.genes$g2m.genes)
DimPlot(sub_cell,reduction = 'umap',group.by = 'Phase',label = TRUE,
        pt.size = 1.5)

sub_cell<-readRDS('T_runpca_control1.rds')
sub_cell<-FindNeighbors(sub_cell,dims = 1:18)
sub_cell<-FindClusters(sub_cell,resolution = 0.2)
sub_cell<-RunUMAP(sub_cell,dims = 1:18)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(sub_cell,"T_dim18_res0.2_unnamed.rds")


Idents(sub_cell)<-'seurat_clusters'
DefaultAssay(sub_cell)<-'RNA'
markers<-FindAllMarkers(sub_cell,only.pos = TRUE)

markers.noRio<-markers[!grepl('^RP[SL]',markers$gene,ignore.case = FALSE),]

write.csv(markers.noRio,'sub_cell.markers.noRio.csv')

DotPlot(sub_cell,features = c('CCR7','SELL','LTB',
                              'GZMA','GZMK','CCL5'))+RotatedAxis()

new.cluster.ids<-c('navie CD4+T','effector memory CD8+T')
names(new.cluster.ids) <- levels(sub_cell)
sub_cell <- RenameIdents(sub_cell, new.cluster.ids)
sub_cell@meta.data$sub_celltype<-Idents(sub_cell)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(sub_cell,'T_dim18_res0.2_named.rds')

#Build cds object
expr_matrix<-as(as.matrix(sub_cell@assays$RNA$counts),'sparseMatrix')
p_data<-sub_cell@meta.data
f_data<-data.frame(gene_short_name = row.names(sub_cell),
                   row.names = row.names(sub_cell))
pd<-new('AnnotatedDataFrame',p_data)
fd<-new('AnnotatedDataFrame',f_data)
cds<-newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                    lowerDetectionLimit = 0.5,
                    expressionFamily = negbinomial.size())

cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
cds<-detectGenes(cds,min_expr = 0.1)
expressed_gene<-row.names(subset(fData(cds),num_cells_expressed >= 10))
plan('sequential')
diff<-differentialGeneTest(cds[expressed_gene,],
                           fullModelFormulaStr = '~sub_celltype',
                           cores = 1)
deg<-subset(diff,qval< 0.01)
deg<-deg[order(deg$qval,decreasing = FALSE),]

save(diff,deg,file = 'train_monocle_genes.Rdata')
# load('train_monocle_genes.Rdata')

write.csv(deg,'train.monocle.deg.csv')

ordergene<-row.names(deg)[1:500]
cds<-setOrderingFilter(cds,ordergene)
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,method = 'DDRTree')

cds<-orderCells(cds)

#Visualization --- trajectory
plot_cell_trajectory(cds,color_by = 'Pseudotime',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'State',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
  facet_wrap('~sub_celltype',nrow = 1)
p1<-plot_cell_trajectory(cds,x = 1,y = 2,color_by = 'sub_celltype')+
  theme(legend.position = 'none',panel.border = element_blank())
p2<-plot_complex_cell_trajectory(cds,x = 1,y = 2,color_by = 'sub_celltype')+
  theme(legend.title = element_blank())
p1|p2

plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
  scale_color_npg()
plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
  scale_color_nejm()
# plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
#   scale_color_manual(values = )

#Visualization --- time density
data1<-pData(cds)
ggplot(data1,aes(x=Pseudotime,colour = sub_celltype,fill = sub_celltype))+
  geom_density(bw = 0.5,size=1,alpha=0.5)+theme_classic2()

#Visualization --- genes from different states
s.genes<-c('CCR7','SELL','LTB',
           'GZMA','GZMK','CCL5')
p1<-plot_genes_in_pseudotime(cds[s.genes,],color_by = 'sub_celltype')
p2<-plot_genes_violin(cds[s.genes,],color_by = 'sub_celltype')
p1|p2

#Visualization --- trajectory of CCR7 gene
colnames(pData(cds))
pData(cds)$CCR7 = log2(exprs(cds)['CCR7',]+1)
plot_cell_trajectory(cds,color_by = 'CCR7')

saveRDS(cds,'cds_final_all_deg.rds')

#Pseudotime differential genes
Time_diff<-differentialGeneTest(cds[ordergene,],cores = 1,
                                fullModelFormulaStr = '~sm.ns(Pseudotime)')
write.csv(Time_diff,'Time_diff_500deg.csv')
Time_gene<-Time_diff %>% pull(gene_short_name) %>% as.character()
p<-plot_pseudotime_heatmap(cds[Time_gene,],num_clusters = 4,show_rownames = TRUE,
                           return_heatmap = TRUE)

#Different cluster genes
p$tree_row
cluster_gene<-cutree(p$tree_row,k = 4)
clustering<-data.frame(cluster_gene)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-'Gene_clusters'
write.csv(clustering,'Time_clusters_genes_500deg.csv')

hp.gene<-p$tree_row$labels[p$tree_row$order]
Time_diff_sig<-Time_diff[hp.gene,c('gene_short_name','pval','qval')]
Time_diff_sig<-subset(Time_diff_sig,qval<0.1)
write.csv(Time_diff_sig,'Time_diff_sig_500deg.csv')


#Trajectory branch
plot_cell_trajectory(cds,color_by = 'State',size=1,show_backbone = TRUE)
BEAM_res<-BEAM(cds[ordergene,],branch_point = 4,cores = 2,
               progenitor_method = 'duplicate')
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
BEAM_res<-BEAM_res[,c('gene_short_name','pval','qval')]
write.csv(BEAM_res,'BEAM_res.csv')


plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval<1e-4)),],
                            branch_point = 4,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)
#Visualization --- top100 genes
BEAM_genes<-top_n(BEAM_res,n=100,desc(qval))%>%pull(gene_short_name)%>%as.character()
p<-plot_genes_branched_heatmap(cds[BEAM_genes,],branch_point = 4,show_rownames = TRUE,
                               return_heatmap = TRUE,num_clusters = 4)
hp.gene<-p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig<-BEAM_res[hp.gene,c('gene_short_name','pval','qval')]
write.csv(BEAM_sig,'BEAM_sig.csv')


#T cell grouop monocle2 --- AS1
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/monocle/AS1_T/')
Idents(AS1)<-'celltype'
sub_cell<-subset(AS1,idents=c('T'))
sub_cell<-NormalizeData(sub_cell)
sub_cell<-FindVariableFeatures(sub_cell)
sub_cell<-ScaleData(sub_cell)
sub_cell<-RunPCA(sub_cell)

saveRDS(sub_cell,"T_runpca_AS1.rds")

DimPlot(sub_cell,reduction = 'pca',pt.size = 1)
DimHeatmap(sub_cell,cells = 500,dims = 1:30,balanced = TRUE)
ElbowPlot(sub_cell,ndims = 50)

pct <- sub_cell[["pca"]]@stdev / sum(sub_cell[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1],
           sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
                decreasing = T)[1] + 1)
ElbowPlot(sub_cell,ndims = 50)$data %>% ggplot()+
  geom_point(aes(x=dims,y=stdev))+
  geom_vline(xintercept = pcs,color = 'darkred')+
  theme_bw()+labs(title = 'Elbow plot : quantitative approach')
plot_df<-data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df,aes(cumu,pct,label = rank,color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

par(mfrow=c(1,1))
ptm = Sys.time()
for(i in seq(0.1,1,by=0.1)){
  sub_cell<-FindNeighbors(sub_cell,dims = 1:12)
  sub_cell<-FindClusters(sub_cell,resolution = i)
  sub_cell<-RunUMAP(sub_cell,dims = 1:12)
  picture1<-DimPlot(sub_cell,reduction = 'umap',label = TRUE,pt.size = 1.5)
  ggsave(filename = paste0('dim18_res',i,'.png'),plot = picture1,
         path = 'D:/GEO_test/GSE198339/control1_AS1_merge/monocle/AS1_T')
}
Sys.time() - ptm
saveRDS(sub_cell,"T_dim12_res0.1-1.rds")

clustree(sub_cell@meta.data,prefix = 'RNA_snn_res.')

sub_cell<-CellCycleScoring(object = sub_cell,s.features = cc.genes$s.genes,
                           g2m.features = cc.genes$g2m.genes)
DimPlot(sub_cell,reduction = 'umap',group.by = 'Phase',label = TRUE,
        pt.size = 1.5)

sub_cell<-readRDS('T_runpca_AS1.rds')
sub_cell<-FindNeighbors(sub_cell,dims = 1:18)
sub_cell<-FindClusters(sub_cell,resolution = 0.5)
sub_cell<-RunUMAP(sub_cell,dims = 1:18)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(sub_cell,"T_dim12_res0.5_unnamed.rds")


Idents(sub_cell)<-'seurat_clusters'
DefaultAssay(sub_cell)<-'RNA'
markers<-FindAllMarkers(sub_cell,only.pos = TRUE)

markers.noRio<-markers[!grepl('^RP[SL]',markers$gene,ignore.case = FALSE),]

write.csv(markers.noRio,'sub_cell.markers.noRio.csv')

DotPlot(sub_cell,features = c('CCR7','SELL','LTB',
                              'GZMA','GZMK','CCL5'))+RotatedAxis()

new.cluster.ids<-c('navie CD4+T','effector memory CD8+T')
names(new.cluster.ids) <- levels(sub_cell)
sub_cell <- RenameIdents(sub_cell, new.cluster.ids)
sub_cell@meta.data$sub_celltype<-Idents(sub_cell)
DimPlot(sub_cell,reduction = 'umap',pt.size = 1.5,label = TRUE)
saveRDS(sub_cell,'T_dim18_res0.2_named.rds')

#Build cds object
expr_matrix<-as(as.matrix(sub_cell@assays$RNA$counts),'sparseMatrix')
p_data<-sub_cell@meta.data
f_data<-data.frame(gene_short_name = row.names(sub_cell),
                   row.names = row.names(sub_cell))
pd<-new('AnnotatedDataFrame',p_data)
fd<-new('AnnotatedDataFrame',f_data)
cds<-newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                    lowerDetectionLimit = 0.5,
                    expressionFamily = negbinomial.size())

cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
cds<-detectGenes(cds,min_expr = 0.1)
expressed_gene<-row.names(subset(fData(cds),num_cells_expressed >= 10))
plan('sequential')
diff<-differentialGeneTest(cds[expressed_gene,],
                           fullModelFormulaStr = '~sub_celltype',
                           cores = 1)
deg<-subset(diff,qval< 0.01)
deg<-deg[order(deg$qval,decreasing = FALSE),]

save(diff,deg,file = 'train_monocle_genes.Rdata')
# load('train_monocle_genes.Rdata')

write.csv(deg,'train.monocle.deg.csv')
# deg<-read.csv('train.monocle.deg.csv',row.names = 1,header = TRUE)

ordergene<-row.names(deg)
cds<-setOrderingFilter(cds,ordergene)
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,method = 'DDRTree')

cds<-orderCells(cds)

#Visualization --- trajectory
plot_cell_trajectory(cds,color_by = 'Pseudotime',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'State',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
  facet_wrap('~sub_celltype',nrow = 1)
p1<-plot_cell_trajectory(cds,x = 1,y = 2,color_by = 'sub_celltype')+
  theme(legend.position = 'none',panel.border = element_blank())
p2<-plot_complex_cell_trajectory(cds,x = 1,y = 2,color_by = 'sub_celltype')+
  theme(legend.title = element_blank())
p1|p2

plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
  scale_color_npg()
plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
  scale_color_nejm()
# plot_cell_trajectory(cds,color_by = 'sub_celltype',size=1,show_backbone = TRUE)+
#   scale_color_manual(values = )

#Visualization --- time density
data1<-pData(cds)
ggplot(data1,aes(x=Pseudotime,colour = sub_celltype,fill = sub_celltype))+
  geom_density(bw = 0.5,size=1,alpha=0.5)+theme_classic2()

#Visualization --- genes from different states
s.genes<-c('CCR7','SELL','LTB',
           'GZMA','GZMK','CCL5')
p1<-plot_genes_in_pseudotime(cds[s.genes,],color_by = 'sub_celltype')
p2<-plot_genes_violin(cds[s.genes,],color_by = 'sub_celltype')
p1|p2

#Visualization --- trajectory of CCR7 gene
colnames(pData(cds))
pData(cds)$CCR7 = log2(exprs(cds)['CCR7',]+1)
plot_cell_trajectory(cds,color_by = 'CCR7')

saveRDS(cds,'cds_final_all_deg.rds')

#Pseudotime differential genes
Time_diff<-differentialGeneTest(cds[ordergene,],cores = 1,
                                fullModelFormulaStr = '~sm.ns(Pseudotime)')
write.csv(Time_diff,'Time_diff_500deg.csv')
Time_gene<-Time_diff %>% pull(gene_short_name) %>% as.character()
p<-plot_pseudotime_heatmap(cds[Time_gene,],num_clusters = 4,show_rownames = TRUE,
                           return_heatmap = TRUE)

#Different cluster genes
p$tree_row
cluster_gene<-cutree(p$tree_row,k = 4)
clustering<-data.frame(cluster_gene)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-'Gene_clusters'
write.csv(clustering,'Time_clusters_genes_500deg.csv')

hp.gene<-p$tree_row$labels[p$tree_row$order]
Time_diff_sig<-Time_diff[hp.gene,c('gene_short_name','pval','qval')]
Time_diff_sig<-subset(Time_diff_sig,qval<0.1)
write.csv(Time_diff_sig,'Time_diff_sig_500deg.csv')

#Trajectory branch
plot_cell_trajectory(cds,color_by = 'State',size=1,show_backbone = TRUE)
BEAM_res<-BEAM(cds[ordergene,],branch_point = 1,cores = 2,
               progenitor_method = 'duplicate')
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
BEAM_res<-BEAM_res[,c('gene_short_name','pval','qval')]
write.csv(BEAM_res,'BEAM_res.csv')

#Visualization --- all genes
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval<1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)

#Visualization --- top100 genes
BEAM_genes<-top_n(BEAM_res,n=100,desc(qval))%>%pull(gene_short_name)%>%as.character()
p<-plot_genes_branched_heatmap(cds[BEAM_genes,],branch_point = 1,show_rownames = TRUE,
                               return_heatmap = TRUE,num_clusters = 4)
hp.gene<-p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig<-BEAM_res[hp.gene,c('gene_short_name','pval','qval')]
write.csv(BEAM_sig,'BEAM_sig.csv')




#SCENIC
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(Matrix)
library(SeuratData)
library(pheatmap)
library(scales)
library(ggrepel)

options(scipen = 20000)
windowsFonts(A = windowsFont('Times New Roman'),
             B = windowsFont('Arial'))


#Convert Seurat object to Anndata object
library(SeuratDisk)
library(sceasy)
use_python('D:/Miniconda3-2024/envs/sceasy/python.exe',required = TRUE)
library(reticulate)
use_condaenv('sceasy')
py_config()
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/scenic')
control1<-readRDS('D:/GEO_test/GSE198339/control1/control1_dim11_res0.5_unnamed_Singlet.rds')
control1[['RNA']]<-as(control1[['RNA']],'Assay')
convertFormat(control1,from = 'seurat',to = 'anndata',outFile = 'control1.h5ad')

AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_dim12_res0.5_unnamed_Singlet.rds')
AS1[['RNA']]<-as(AS1[['RNA']],'Assay')
convertFormat(AS1,from = 'seurat',to = 'anndata',outFile = 'AS1.h5ad')

#make XXX.loom  --- in jupyter notebook
# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import json
import zlib
import base64

import umap

wdir = "D:/GEO_test/GSE198339/control1_AS1_merge/scenic"
os.chdir(wdir)

adata_control1 = sc.read_h5ad('control1.h5ad')
adata_AS1 = sc.read_h5ad('AS1.h5ad')

#control1.loom
row_attrs = {
  "Gene": np.array(adata_control1.var.index) ,
}
col_attrs = {
  "CellID": np.array(adata_control1.obs.index) ,
  "nGene": np.array( np.sum(adata_control1.X.transpose()>0 , axis=0)).flatten() ,
  "nUMI": np.array( np.sum(adata_control1.X.transpose() , axis=0)).flatten() ,
}
lp.create('control1_pyscenic.loom',adata_control1.X.transpose(), row_attrs, col_attrs)

#AS1.loom
row_attrs = {
  "Gene": np.array(adata_AS1.var.index) ,
}
col_attrs = {
  "CellID": np.array(adata_AS1.obs.index) ,
  "nGene": np.array( np.sum(adata_AS1.X.transpose()>0 , axis=0)).flatten() ,
  "nUMI": np.array( np.sum(adata_AS1.X.transpose() , axis=0)).flatten() ,
}
lp.create('AS1_pyscenic.loom',adata_AS1.X.transpose(), row_attrs, col_attrs)

#pyscenic --- in Windows command line
###First step
pyscenic grn --num_workers 6 --output control1_adj.tsv --method grnboost2 control1_pyscenic.loom tf_list_human.txt

pyscenic grn --num_workers 6 --output AS1_adj.tsv --method grnboost2 AS1_pyscenic.loom tf_list_human.txt



###Second step
pyscenic ctx control1_adj.tsv D:/SCENIC_database/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname D:/SCENIC_database/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname control1_pyscenic.loom --mode "dask_multiprocessing" --output control1_500bp_reg.csv --num_workers 6 --mask_dropouts

pyscenic ctx AS1_adj.tsv D:/SCENIC_database/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname D:/SCENIC_database/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname AS1_pyscenic.loom --mode "dask_multiprocessing" --output AS1_500bp_reg.csv --num_workers 6 --mask_dropouts


###Third step
pyscenic aucell control1_pyscenic.loom control1_500bp_reg.csv --output control1_500bp_scenic.loom --num_workers 6

pyscenic aucell AS1_pyscenic.loom AS1_500bp_reg.csv --output AS1_500bp_scenic.loom --num_workers 6


#Intergrate output --- in jupyter notebook
#collect SCENIC AUCell output --- AS1
f_pyscenic_output = 'D:/GEO_test/GSE198339/control1_AS1_merge/scenic/AS1_500bp_scenic.loom'
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx.to_csv('AS1_auc_500bp.csv',sep = '\t')
lf.close()

#UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation',).fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_umap_AS1_500bp.txt", sep='\t')
#tSNE
tsne = TSNE( n_jobs=5 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_tsne_AS1_500bp.txt", sep='\t')

#scenic output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'scenic_umap_AS1_500bp.txt', sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv( 'scenic_tsne_AS1_500bp.txt', sep='\t', header=0, index_col=0 )

auc_mtx.columns = auc_mtx.columns.str.replace('(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
#regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
  tmp = x.get('regulon').replace("(","_(")
x.update( {'regulon': tmp} )

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
  pd.DataFrame(adata_AS1.obsm['X_umap'],index=adata_AS1.obs.index)[0] ,
  pd.DataFrame(adata_AS1.obsm['X_pca'],index=adata_AS1.obs.index)[0] ,
  dr_tsne['X'] ,
  dr_umap['X']
], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
  pd.DataFrame(adata_AS1.obsm['X_umap'],index=adata_AS1.obs.index)[1] ,
  pd.DataFrame(adata_AS1.obsm['X_pca'],index=adata_AS1.obs.index)[1] ,
  dr_tsne['Y'] ,
  dr_umap['Y']
], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3','4']

metaJson = {}

metaJson['embeddings'] = [
  {
    "id": 1,
    "name": f"Scanpy UMAP  (highly variable genes)"
  },
  {
    "id": 2,
    "name": "Scanpy PC1/PC2"
  },
  {
    "id": 3,
    "name": "SCENIC AUC t-SNE"
  },
  {
    "id": 4,
    "name": "SCENIC AUC UMAP"
  },
]

metaJson["clusterings"] = [{
  "id": 0,
  "group": "Scanpy",
  "name": "Scanpy louvain default resolution",
  "clusters": [],
}]

metaJson["metrics"] = [
  {
    "name": "nUMI"
  }, {
    "name": "nGene"
  }, {
    "name": "Percent_mito"
  }
]

metaJson["annotations"] = [
  {
    "name": "Louvain_clusters_Scanpy",
    "values": list(set( adata_AS1.obs['seurat_clusters'].astype(np.str) ))
  },
  #{
  #    "name": "Genotype",
  #    "values": list(set(adata.obs['Genotype'].values))
  #},
  #{
  #    "name": "Timepoint",
  #    "values": list(set(adata.obs['Timepoint'].values))
  #},
  #{
  #    "name": "Sample",
  #    "values": list(set(adata.obs['Sample'].values))
  #}
]

#SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in adata_AS1.obs['seurat_clusters']])) + 1):
  clustDict = {}
clustDict['id'] = i
clustDict['description'] = f'Unannotated Cluster {i + 1}'
metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()
clusterings["0"] = adata_AS1.obs['seurat_clusters'].values.astype(np.int64)

def dfToNamedMatrix(df):
  arr_ip = [tuple(i) for i in df.values]
dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
arr = np.array(arr_ip, dtype=dtyp)
return arr

col_attrs = {
  "CellID": np.array(adata_AS1.obs.index),
  "nUMI": np.array(adata_AS1.obs['nCount_RNA'].values),
  "nGene": np.array(adata_AS1.obs['nFeature_RNA'].values),
  "seurat_clusters_Scanpy": np.array( adata_AS1.obs['seurat_clusters'].values ),
  #"Genotype": np.array(adata.obs['Genotype'].values),
  #"Timepoint": np.array(adata.obs['Timepoint'].values),
  #"Sample": np.array(adata.obs['Sample'].values),
  "Percent_mito": np.array(adata_AS1.obs['percent.mt'].values),
  "Embeddings_X": dfToNamedMatrix(Embeddings_X),
  "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
  "RegulonsAUC": dfToNamedMatrix(auc_mtx),
  "Clusterings": dfToNamedMatrix(clusterings),
  "ClusterID": np.array(adata_AS1.obs['seurat_clusters'].values)
}

row_attrs = {
  "Gene": lf.ra.Gene,
  "Regulons": regulons,
}

attrs = {
  "title": "sampleTitle",
  "MetaData": json.dumps(metaJson),
  "Genome": 'hg38',
  "SCopeTreeL1": "",
  "SCopeTreeL2": "",
  "SCopeTreeL3": ""
}

#compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

f_final_loom = 'AS1_500bp_scenic_integrated-output.loom'
lp.create(
  filename = f_final_loom ,
  layers=lf[:,:],
  row_attrs=row_attrs, 
  col_attrs=col_attrs, 
  # file_attrs=attrs
)

lf.close()

#collect SCENIC AUCell output --- AS1
f_pyscenic_output = 'D:/GEO_test/GSE198339/control1_AS1_merge/scenic/control1_500bp_scenic.loom'
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx.to_csv('control1_auc_500bp.csv',sep = '\t')
lf.close()

#UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation',).fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_umap_control1_500bp.txt", sep='\t')
#tSNE
tsne = TSNE( n_jobs=5 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_tsne_control1_500bp.txt", sep='\t')

#scenic output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'scenic_umap_control1_500bp.txt', sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv( 'scenic_tsne_control1_500bp.txt', sep='\t', header=0, index_col=0 )

auc_mtx.columns = auc_mtx.columns.str.replace('(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
#regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
  tmp = x.get('regulon').replace("(","_(")
x.update( {'regulon': tmp} )

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
  pd.DataFrame(adata_control1.obsm['X_umap'],index=adata_control1.obs.index)[0] ,
  pd.DataFrame(adata_control1.obsm['X_pca'],index=adata_control1.obs.index)[0] ,
  dr_tsne['X'] ,
  dr_umap['X']
], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
  pd.DataFrame(adata_control1.obsm['X_umap'],index=adata_control1.obs.index)[1] ,
  pd.DataFrame(adata_control1.obsm['X_pca'],index=adata_control1.obs.index)[1] ,
  dr_tsne['Y'] ,
  dr_umap['Y']
], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3','4']

metaJson = {}

metaJson['embeddings'] = [
  {
    "id": 1,
    "name": f"Scanpy UMAP  (highly variable genes)"
  },
  {
    "id": 2,
    "name": "Scanpy PC1/PC2"
  },
  {
    "id": 3,
    "name": "SCENIC AUC t-SNE"
  },
  {
    "id": 4,
    "name": "SCENIC AUC UMAP"
  },
]

metaJson["clusterings"] = [{
  "id": 0,
  "group": "Scanpy",
  "name": "Scanpy louvain default resolution",
  "clusters": [],
}]

metaJson["metrics"] = [
  {
    "name": "nUMI"
  }, {
    "name": "nGene"
  }, {
    "name": "Percent_mito"
  }
]

metaJson["annotations"] = [
  {
    "name": "Louvain_clusters_Scanpy",
    "values": list(set( adata_control1.obs['seurat_clusters'].astype(np.str) ))
  },
  #{
  #    "name": "Genotype",
  #    "values": list(set(adata.obs['Genotype'].values))
  #},
  #{
  #    "name": "Timepoint",
  #    "values": list(set(adata.obs['Timepoint'].values))
  #},
  #{
  #    "name": "Sample",
  #    "values": list(set(adata.obs['Sample'].values))
  #}
]

#SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in adata_control1.obs['seurat_clusters']])) + 1):
  clustDict = {}
clustDict['id'] = i
clustDict['description'] = f'Unannotated Cluster {i + 1}'
metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()
clusterings["0"] = adata_control1.obs['seurat_clusters'].values.astype(np.int64)

def dfToNamedMatrix(df):
  arr_ip = [tuple(i) for i in df.values]
dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
arr = np.array(arr_ip, dtype=dtyp)
return arr

col_attrs = {
  "CellID": np.array(adata_control1.obs.index),
  "nUMI": np.array(adata_control1.obs['nCount_RNA'].values),
  "nGene": np.array(adata_control1.obs['nFeature_RNA'].values),
  "seurat_clusters_Scanpy": np.array( adata_control1.obs['seurat_clusters'].values ),
  #"Genotype": np.array(adata.obs['Genotype'].values),
  #"Timepoint": np.array(adata.obs['Timepoint'].values),
  #"Sample": np.array(adata.obs['Sample'].values),
  "Percent_mito": np.array(adata_control1.obs['percent.mt'].values),
  "Embeddings_X": dfToNamedMatrix(Embeddings_X),
  "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
  "RegulonsAUC": dfToNamedMatrix(auc_mtx),
  "Clusterings": dfToNamedMatrix(clusterings),
  "ClusterID": np.array(adata_control1.obs['seurat_clusters'].values)
}

row_attrs = {
  "Gene": lf.ra.Gene,
  "Regulons": regulons,
}

attrs = {
  "title": "sampleTitle",
  "MetaData": json.dumps(metaJson),
  "Genome": 'hg38',
  "SCopeTreeL1": "",
  "SCopeTreeL2": "",
  "SCopeTreeL3": ""
}

#compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

f_final_loom = 'control1_500bp_scenic_integrated-output.loom'
lp.create(
  filename = f_final_loom ,
  layers=lf[:,:],
  row_attrs=row_attrs, 
  col_attrs=col_attrs, 
  file_attrs=attrs
)

lf.close()

#scenic转录因子可视化---control1
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/scenic/control1')
control1<-readRDS('D:/GEO_test/GSE198339/control1/control1_dim11_res0.5_named_Singlet.rds')

loom_control1<-open_loom('D:/GEO_test/GSE198339/control1_AS1_merge/scenic/control1_500bp_scenic_integrated-output.loom') 

regulons_incidMat_control1<-get_regulons(loom_control1, column.attr.name="Regulons")
regulons_incidMat_control1[1:4,1:4] 
regulons_control1<-regulonsToGeneLists(regulons_incidMat_control1)

regulons_control1_df<-data.frame()
for(i in 1:length(regulons_control1)){
  data_test1<-data.frame(regulons_control1[i])
  names(data_test1)<-gsub("\\.", "-", names(data_test1))
  name_middle<-paste0(sub("_(.*)$", "", names(data_test1)),' (',
                      length(row.names(data_test1)),'g)')
  data_test1$regulon_name<-rep(name_middle,length(row.names(data_test1)))
  names(data_test1)[1]<-'regulon_genes'
  regulons_control1_df<-rbind.data.frame(regulons_control1_df,data_test1)
}
fwrite(x = regulons_control1_df,file = 'regulons_500bp_control1.csv',row.names = FALSE)

regulonAUC_control1<-get_regulons_AUC(loom_control1,column.attr.name='RegulonsAUC')
regulonAucThresholds_control1<-get_regulon_thresholds(loom_control1)

save(regulons_incidMat_control1,regulons_control1,regulonAUC_control1,
     regulonAucThresholds_control1,
     regulons_control1_df,
     file = 'regulon_AUC_thresholds_500bp_control1.Rdata')
# load('regulon_AUC_thresholds_500bp_control1.Rdata')

embeddings<-get_embeddings(loom_control1)  
close_loom(loom_control1)

rownames(regulonAUC_control1)
colnames(regulonAUC_control1)
names(regulons_control1)

control1$sub_celltype<-control1$celltype

Idents(control1)<-control1$sub_celltype

sub_regulonAUC_control1<-regulonAUC_control1[,match(colnames(control1),colnames(regulonAUC_control1))]
row.names(sub_regulonAUC_control1)<-unique(regulons_control1_df$regulon_name)
row.names(sub_regulonAUC_control1@assays@data@listData$AUC)<-unique(regulons_control1_df$regulon_name)
dim(sub_regulonAUC_control1)

identical(colnames(sub_regulonAUC_control1),colnames(control1))

cellClusters_control1<-data.frame(row.names = colnames(control1), 
                                  seurat_clusters = as.character(control1$seurat_clusters))
cellTypes_control1<-data.frame(row.names = colnames(control1), 
                               celltype = control1$sub_celltype)
# head(cellTypes_control1)
# head(cellClusters_control1)
# sub_regulonAUC_control1[1:4,1:4]

save(sub_regulonAUC_control1,cellTypes_control1,cellClusters_control1,control1,
     file = 'for_rss_and_visual_500bp_control1.Rdata')

# load('for_rss_and_visual_500bp_control1.Rdata')

regulonsToPlot<-c('EOMES (79g)')
regulonsToPlot
control1@meta.data<-data.frame(control1@meta.data,
                               sub_regulonAUC_control1@assays@data$AUC[regulonsToPlot,])
colnames(control1@meta.data)[length(colnames(control1@meta.data))]<-regulonsToPlot
Idents(control1)<-control1$sub_celltype
table(Idents(control1))

DotPlot(control1, features = unique(regulonsToPlot))+ RotatedAxis()
RidgePlot(control1, features = regulonsToPlot , ncol = 1)
VlnPlot(control1, features = regulonsToPlot,pt.size = 0 ) 
FeaturePlot(control1, features = regulonsToPlot)

#ras_mean
cellsPerGroup<-split(rownames(cellTypes_control1), 
                     cellTypes_control1[,'celltype']) 
sub_regulonAUC_control1<-sub_regulonAUC_control1[onlyNonDuplicatedExtended(rownames(sub_regulonAUC_control1)),] 
dim(sub_regulonAUC_control1)

regulonActivity_byGroup<-sapply(cellsPerGroup,
                                function(cells) 
                                  rowMeans(getAUC(sub_regulonAUC_control1)[,cells]))

regulonActivity_byGroup_Scaled<-t(scale(t(regulonActivity_byGroup),
                                        center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

png('ras_control1_500bp.png',width = 1400,height = 5200,res = 300)
pheatmap(regulonActivity_byGroup_Scaled,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue","white", "firebrick3"))(50),
         cluster_cols = TRUE,cluster_rows = TRUE,
         border=FALSE,
         fontsize = 5,cellwidth = 20,cellheight = 5,
         annotation_names_col = TRUE,annotation_legend = TRUE,
         fontfamily = 'A')
dev.off()

save(cellsPerGroup,regulonActivity_byGroup,regulonActivity_byGroup_Scaled,
     file = 'ras_mean_500bp_control1.Rdata')

# load('ras_mean_500bp_control1.Rdata')

#rss
rss_control1<-calcRSS(AUC=getAUC(sub_regulonAUC_control1), 
                      cellAnnotation=cellTypes_control1[colnames(sub_regulonAUC_control1),'celltype'])
rss_control1=na.omit(rss_control1) 

rss_control1_top5<-data.frame()
for(i in 1:length(colnames(rss_control1))){
  rss_control1_test1<-data.frame(rss_control1[,i])
  rss_control1_test1$regulons<-row.names(rss_control1)
  rss_control1_test1$celltype<-rep(colnames(rss_control1)[i],times=length(row.names(rss_control1_test1)))
  rss_control1_test2<-rss_control1_test1[order(rss_control1_test1[,1],decreasing = TRUE),]
  names(rss_control1_test2)[1]<-'rss_values'
  rss_control1_top5<-rbind.data.frame(rss_control1_top5,rss_control1_test2[1:5,])
}
fwrite(x = rss_control1_top5,file = 'rss_control1_500bp_top5.csv',row.names = FALSE)

save(rss_control1,rss_control1_top5,file = 'rss_500bp_control1.Rdata')
# load('rss_500bp_control1.Rdata')

#binary
regulonAucThresholds_control1_df<-data.frame(regulonAucThresholds_control1)
head(regulonAucThresholds_control1_df)

regulonAucThresholds_control1_df$thresholds<-as.numeric(row.names(regulonAucThresholds_control1_df))
row.names(regulonAucThresholds_control1_df)<-as.numeric(1:length(row.names(regulonAucThresholds_control1_df)))
names(regulonAucThresholds_control1_df)[1]<-'regulons'
regulonAucThresholds_control1_df$regulons<-unique(regulons_control1_df$regulon_name)
sub_regulonAUC_control1_df<-as.data.frame(sub_regulonAUC_control1@assays@data@listData$AUC)

head(regulonAucThresholds_control1_df)
sub_regulonAUC_control1_df[1:5,1:5]

if(identical(row.names(sub_regulonAUC_control1_df),
             regulonAucThresholds_control1_df$regulons)){
  print("Consistent of content and order")
}else{
  print("Inconsistent of content and order")
}

regulonsAUC_control1_binary_1<-cbind.data.frame(sub_regulonAUC_control1_df,
                                                regulonAucThresholds_control1_df[,2])
names(regulonsAUC_control1_binary_1)[length(regulonsAUC_control1_binary_1)]<-'thresholds'
regulonsAUC_control1_binary_1_row<-regulonsAUC_control1_binary_1

regulonsAUC_control1_binary_1[1:5,1:5]
regulonsAUC_control1_binary_1_row[1:5,1:5]
regulonsAUC_control1_binary_1_row$thresholds[1:5]

for (i in 1:(ncol(regulonsAUC_control1_binary_1)-1)) {
  regulonsAUC_control1_binary_1[,i]<-ifelse(regulonsAUC_control1_binary_1[,i] < regulonsAUC_control1_binary_1$thresholds,0,1)
}
regulonsAUC_control1_binary_1<-regulonsAUC_control1_binary_1[,-length(colnames(regulonsAUC_control1_binary_1))]
regulonsAUC_control1_binary_1[1:5,1:5]

fwrite(x = regulonsAUC_control1_binary_1,row.names = TRUE,
       file = 'regulonsAUC_control1_binary_500bp.csv')

regulonsAUC_control1_binary_1_t<-data.frame(t(regulonsAUC_control1_binary_1))
regulonsAUC_control1_binary_1_t[1:5,1:5]
colnames(regulonsAUC_control1_binary_1_t)<-row.names(regulonsAUC_control1_binary_1)

pheatmap(regulonsAUC_control1_binary_1_t[1:10,1:10],
         show_rownames = TRUE,show_colnames = TRUE,
         color = c('white','black'),
         cluster_cols = TRUE,cluster_rows = TRUE,
         border=FALSE,
         fontsize = 10,cellwidth = 10,cellheight = 10,
         annotation_names_col = TRUE,annotation_legend = TRUE,
         fontfamily = 'A')

#UMI group
cellTypes_control1_again<-cellTypes_control1%>%arrange(celltype)
head(cellTypes_control1_again)
regulonsAUC_control1_binary_2<-regulonsAUC_control1_binary_1_t
regulonsAUC_control1_binary_2[1:5,1:5]

regulonsAUC_control1_binary_3<-regulonsAUC_control1_binary_2[row.names(cellTypes_control1_again),]
regulonsAUC_control1_binary_3[1:5,1:5]
if(identical(row.names(regulonsAUC_control1_binary_3),
             row.names(cellTypes_control1_again))){
  print("Consistent of content and order")
}else{
  print("Inconsistent of content and order")
}

regulonsAUC_control1_binary_4<-as.data.frame(t(regulonsAUC_control1_binary_3))
regulonsAUC_control1_binary_4[1:5,1:5]
regulonsAUC_control1_binary_group<-data.frame(table(cellTypes_control1_again$celltype))
head(regulonsAUC_control1_binary_group)
regulonsAUC_control1_binary_group_again<-data.frame()
for(i in 1:length(row.names(regulonsAUC_control1_binary_group))){
  data_middle_1<-data.frame(celltype=factor(rep(regulonsAUC_control1_binary_group$Var1[i],
                                                each=regulonsAUC_control1_binary_group$Freq[i])))
  regulonsAUC_control1_binary_group_again<-rbind.data.frame(regulonsAUC_control1_binary_group_again,
                                                            data_middle_1)
}

if(identical(regulonsAUC_control1_binary_group_again$celltype,
             cellTypes_control1_again$celltype)){
  print("Consistent of content and order")
}else{
  print("Inconsistent of content and order")
}

row.names(regulonsAUC_control1_binary_group_again)<-colnames(regulonsAUC_control1_binary_4)
regulonsAUC_control1_binary_group_again$celltype<-as.factor(regulonsAUC_control1_binary_group_again$celltype)
levels(regulonsAUC_control1_binary_group_again$celltype)

png('auc_binary_control1_500bp.png',width = 2000,height = 5000,res = 300)
pheatmap(regulonsAUC_control1_binary_4,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = c('white','black'),
         cluster_cols = TRUE,border=FALSE,
         annotation_col = regulonsAUC_control1_binary_group_again,
         cluster_rows = TRUE,fontsize = 5,cellwidth = 0.1,cellheight = 5,
         annotation_names_col = TRUE,annotation_legend = TRUE)
dev.off()

save(regulonsAUC_control1_binary_1,regulonsAUC_control1_binary_3,
     regulonsAUC_control1_binary_4,regulonsAUC_control1_binary_group,
     regulonsAUC_control1_binary_group_again,
     file = 'binary_visual_control1_500bp.Rdata')

# load('binary_visual_control1_500bp.Rdata')


#scenic转录因子可视化---AS1
setwd('D:/GEO_test/GSE198339/control1_AS1_merge/scenic/AS1')
AS1<-readRDS('D:/GEO_test/GSE198339/AS1/AS1_dim12_res0.5_named_Singlet.rds')

loom_AS1<-open_loom('D:/GEO_test/GSE198339/control1_AS1_merge/scenic/AS1_500bp_scenic_integrated-output.loom') 

regulons_incidMat_AS1<-get_regulons(loom_AS1, column.attr.name="Regulons")
regulons_incidMat_AS1[1:4,1:4] 
regulons_AS1<-regulonsToGeneLists(regulons_incidMat_AS1)

regulons_AS1_df<-data.frame()
for(i in 1:length(regulons_AS1)){
  data_test1<-data.frame(regulons_AS1[i])
  names(data_test1)<-gsub("\\.", "-", names(data_test1))
  name_middle<-paste0(sub("_(.*)$", "", names(data_test1)),' (',
                      length(row.names(data_test1)),'g)')
  data_test1$regulon_name<-rep(name_middle,length(row.names(data_test1)))
  names(data_test1)[1]<-'regulon_genes'
  regulons_AS1_df<-rbind.data.frame(regulons_AS1_df,data_test1)
}
fwrite(x = regulons_AS1_df,file = 'regulons_500bp_AS1.csv',row.names = FALSE)

regulonAUC_AS1<-get_regulons_AUC(loom_AS1,column.attr.name='RegulonsAUC')
regulonAucThresholds_AS1<-get_regulon_thresholds(loom_AS1)

save(regulons_incidMat_AS1,regulons_AS1,regulonAUC_AS1,
     regulonAucThresholds_AS1,
     regulons_AS1_df,
     file = 'regulon_AUC_thresholds_500bp_AS1.Rdata')
# load('regulon_AUC_thresholds_500bp_AS1.Rdata')

embeddings<-get_embeddings(loom_AS1)  
close_loom(loom_AS1)

# rownames(regulonAUC_AS1)
# colnames(regulonAUC_AS1)
# names(regulons_AS1)

AS1$sub_celltype<-AS1$celltype
Idents(AS1)<-AS1$sub_celltype

sub_regulonAUC_AS1<-regulonAUC_AS1[,match(colnames(AS1),colnames(regulonAUC_AS1))]
row.names(sub_regulonAUC_AS1)<-unique(regulons_AS1_df$regulon_name)
row.names(sub_regulonAUC_AS1@assays@data@listData$AUC)<-unique(regulons_AS1_df$regulon_name)
dim(sub_regulonAUC_AS1)

identical(colnames(sub_regulonAUC_AS1),colnames(AS1))

cellClusters_AS1<-data.frame(row.names = colnames(AS1), 
                             seurat_clusters = as.character(AS1$seurat_clusters))
cellTypes_AS1<-data.frame(row.names = colnames(AS1), 
                          celltype = AS1$sub_celltype)
# head(cellTypes_AS1)
# head(cellClusters_AS1)
# sub_regulonAUC_AS1[1:4,1:4]

save(sub_regulonAUC_AS1,cellTypes_AS1,cellClusters_AS1,AS1,
     file = 'for_rss_and_visual_500bp_AS1.Rdata')

# load('for_rss_and_visual_500bp_AS1.Rdata')

regulonsToPlot<-c('EOMES (119g)')
regulonsToPlot
AS1@meta.data<-data.frame(AS1@meta.data ,
                          sub_regulonAUC_AS1@assays@data$AUC[regulonsToPlot,])
colnames(AS1@meta.data)[length(colnames(AS1@meta.data))]<-regulonsToPlot
Idents(AS1)<-AS1$sub_celltype
table(Idents(AS1))

DotPlot(AS1, features = unique(regulonsToPlot))+ RotatedAxis()
RidgePlot(AS1, features = regulonsToPlot,ncol = 1)
VlnPlot(AS1, features = regulonsToPlot,pt.size = 0) 
FeaturePlot(AS1, features = regulonsToPlot)

#ras_mean
cellsPerGroup<-split(rownames(cellTypes_AS1), 
                     cellTypes_AS1[,'celltype']) 
sub_regulonAUC_AS1<-sub_regulonAUC_AS1[onlyNonDuplicatedExtended(rownames(sub_regulonAUC_AS1)),] 
dim(sub_regulonAUC_AS1)

regulonActivity_byGroup<-sapply(cellsPerGroup,
                                function(cells) 
                                  rowMeans(getAUC(sub_regulonAUC_AS1)[,cells]))

regulonActivity_byGroup_Scaled<-t(scale(t(regulonActivity_byGroup),
                                        center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

png('ras_AS1_500bp.png',width = 1400,height = 5200,res = 300)
pheatmap(regulonActivity_byGroup_Scaled,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue","white", "firebrick3"))(50),
         cluster_cols = TRUE,cluster_rows = TRUE,
         border=FALSE,
         fontsize = 5,cellwidth = 20,cellheight = 5,
         annotation_names_col = TRUE,annotation_legend = TRUE,
         fontfamily = 'A')
dev.off()

save(cellsPerGroup,regulonActivity_byGroup,regulonActivity_byGroup_Scaled,
     file = 'ras_mean_500bp_AS1.Rdata')

# load('ras_mean_500bp_AS1.Rdata')

#rss
rss_AS1<-calcRSS(AUC=getAUC(sub_regulonAUC_AS1), 
                 cellAnnotation=cellTypes_AS1[colnames(sub_regulonAUC_AS1),'celltype'])
rss_AS1=na.omit(rss_AS1) 

rss_AS1_top5<-data.frame()
for(i in 1:length(colnames(rss_AS1))){
  rss_AS1_test1<-data.frame(rss_AS1[,i])
  rss_AS1_test1$regulons<-row.names(rss_AS1)
  rss_AS1_test1$celltype<-rep(colnames(rss_AS1)[i],times=length(row.names(rss_AS1_test1)))
  rss_AS1_test2<-rss_AS1_test1[order(rss_AS1_test1[,1],decreasing = TRUE),]
  names(rss_AS1_test2)[1]<-'rss_values'
  rss_AS1_top5<-rbind.data.frame(rss_AS1_top5,rss_AS1_test2[1:5,])
}
fwrite(x = rss_AS1_top5,file = 'rss_AS1_500bp_top5.csv',row.names = FALSE)

save(rss_AS1,rss_AS1_top5,file = 'rss_500bp_AS1.Rdata')
# load('rss_500bp_AS1.Rdata')

#binary
regulonAucThresholds_AS1_df<-data.frame(regulonAucThresholds_AS1)
head(regulonAucThresholds_AS1_df)

regulonAucThresholds_AS1_df$thresholds<-as.numeric(row.names(regulonAucThresholds_AS1_df))
row.names(regulonAucThresholds_AS1_df)<-as.numeric(1:length(row.names(regulonAucThresholds_AS1_df)))
names(regulonAucThresholds_AS1_df)[1]<-'regulons'
regulonAucThresholds_AS1_df$regulons<-unique(regulons_AS1_df$regulon_name)
sub_regulonAUC_AS1_df<-as.data.frame(sub_regulonAUC_AS1@assays@data@listData$AUC)

head(regulonAucThresholds_AS1_df)
sub_regulonAUC_AS1_df[1:5,1:5]

if(identical(row.names(sub_regulonAUC_AS1_df),
             regulonAucThresholds_AS1_df$regulons)){
  print("Consistent of content and order")
}else{
  print("Inconsistent of content and order")
}

regulonsAUC_AS1_binary_1<-cbind.data.frame(sub_regulonAUC_AS1_df,
                                           regulonAucThresholds_AS1_df[,2])
names(regulonsAUC_AS1_binary_1)[length(regulonsAUC_AS1_binary_1)]<-'thresholds'
regulonsAUC_AS1_binary_1_row<-regulonsAUC_AS1_binary_1

# regulonsAUC_AS1_binary_1[1:5,1:5]
# regulonsAUC_AS1_binary_1_row[1:5,1:5]
# regulonsAUC_AS1_binary_1_row$thresholds[1:5]

for(i in 1:(ncol(regulonsAUC_AS1_binary_1)-1)){
  regulonsAUC_AS1_binary_1[,i]<-ifelse(regulonsAUC_AS1_binary_1[,i] < regulonsAUC_AS1_binary_1$thresholds,0,1)
}
regulonsAUC_AS1_binary_1<-regulonsAUC_AS1_binary_1[,-length(colnames(regulonsAUC_AS1_binary_1))]
regulonsAUC_AS1_binary_1[1:5,1:5]

fwrite(x = regulonsAUC_AS1_binary_1,row.names = TRUE,
       file = 'regulonsAUC_AS1_binary_500bp.csv')

regulonsAUC_AS1_binary_1_t<-data.frame(t(regulonsAUC_AS1_binary_1))
regulonsAUC_AS1_binary_1_t[1:5,1:5]
colnames(regulonsAUC_AS1_binary_1_t)<-row.names(regulonsAUC_AS1_binary_1)

pheatmap(regulonsAUC_AS1_binary_1_t[1:10,1:10],
         show_rownames = TRUE,show_colnames = TRUE,
         color = c('white','black'),
         cluster_cols = TRUE,cluster_rows = TRUE,
         border=FALSE,
         fontsize = 10,cellwidth = 10,cellheight = 10,
         annotation_names_col = TRUE,annotation_legend = TRUE,
         fontfamily = 'A')

#UMI group
cellTypes_AS1_again<-cellTypes_AS1%>%arrange(celltype)
head(cellTypes_AS1_again)
regulonsAUC_AS1_binary_2<-regulonsAUC_AS1_binary_1_t
regulonsAUC_AS1_binary_2[1:5,1:5]

regulonsAUC_AS1_binary_3<-regulonsAUC_AS1_binary_2[row.names(cellTypes_AS1_again),]
regulonsAUC_AS1_binary_3[1:5,1:5]
if(identical(row.names(regulonsAUC_AS1_binary_3),
             row.names(cellTypes_AS1_again))){
  print("Consistent of content and order")
}else{
  print("Inconsistent of content and order")
}

regulonsAUC_AS1_binary_4<-as.data.frame(t(regulonsAUC_AS1_binary_3))
regulonsAUC_AS1_binary_4[1:5,1:5]
regulonsAUC_AS1_binary_group<-data.frame(table(cellTypes_AS1_again$celltype))
head(regulonsAUC_AS1_binary_group)
regulonsAUC_AS1_binary_group_again<-data.frame()
for(i in 1:length(row.names(regulonsAUC_AS1_binary_group))){
  data_middle_1<-data.frame(celltype=factor(rep(regulonsAUC_AS1_binary_group$Var1[i],
                                                each=regulonsAUC_AS1_binary_group$Freq[i])))
  regulonsAUC_AS1_binary_group_again<-rbind.data.frame(regulonsAUC_AS1_binary_group_again,
                                                       data_middle_1)
}

if(identical(regulonsAUC_AS1_binary_group_again$celltype,
             cellTypes_AS1_again$celltype)){
  print("Consistent of content and order")
}else{
  print("Inconsistent of content and order")
}

row.names(regulonsAUC_AS1_binary_group_again)<-colnames(regulonsAUC_AS1_binary_4)
regulonsAUC_AS1_binary_group_again$celltype<-as.factor(regulonsAUC_AS1_binary_group_again$celltype)
levels(regulonsAUC_AS1_binary_group_again$celltype)

png('auc_binary_AS1_500bp.png',width = 2000,height = 5000,res = 300)
pheatmap(regulonsAUC_AS1_binary_4,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = c('white','black'),
         cluster_cols = TRUE,border=FALSE,
         annotation_col = regulonsAUC_AS1_binary_group_again,
         cluster_rows = TRUE,fontsize = 5,cellwidth = 0.1,cellheight = 5,
         annotation_names_col = TRUE,annotation_legend = TRUE)
dev.off()

save(regulonsAUC_AS1_binary_1,regulonsAUC_AS1_binary_3,
     regulonsAUC_AS1_binary_4,regulonsAUC_AS1_binary_group,
     regulonsAUC_AS1_binary_group_again,
     file = 'binary_visual_AS1_500bp.Rdata')
# load('binary_visual_AS1_500bp.Rdata')
