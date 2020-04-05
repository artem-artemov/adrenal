library(Seurat)
library(dplyr)

do_seurat_preprocess=function(env){
    with(env, {

        SR[["percent.mt"]] <- PercentageFeatureSet(SR, pattern = "^mt-")
        VlnPlot(SR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

        plot1 <- FeatureScatter(SR, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(SR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        plot(CombinePlots(plots = list(plot1, plot2)))

        SR <- subset(SR, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25)

        SR <- NormalizeData(SR) #, normalization.method = "LogNormalize", scale.factor = 10000)
    })
}

do_seurat_regress = function(env, regress=c("nFeature_RNA", "percent.mt")){
    env$regress = regress
    with(env, {
        SR[["percent.mt"]] <- PercentageFeatureSet(SR, pattern = "^MT-")
        SR[["percent.hb"]] <- PercentageFeatureSet(SR, pattern = "^HB")

        VlnPlot(SR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        plot1 <- FeatureScatter(SR, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(SR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        plot(CombinePlots(plots = list(plot1, plot2)))

        SR <- subset(SR, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25)

        SR <- NormalizeData(SR)
        
        SR <- FindVariableFeatures(SR, selection.method = "vst", nfeatures = 2000)
        top10 <- head(VariableFeatures(SR), 20)
        plot1 <- VariableFeaturePlot(SR)
        plot(LabelPoints(plot = plot1, points = top10, repel = TRUE))

        all.genes <- rownames(SR)
        
        SR <- ScaleData(SR, vars.to.regress = regress)

        SR <- RunPCA(SR, features = VariableFeatures(object = SR))
        print(SR[["pca"]], dims = 1:5, nfeatures = 5)
        plot(VizDimLoadings(SR, dims = 1:2, reduction = "pca"))

        DimHeatmap(SR, dims = 1:10, cells = 500, balanced = TRUE)

        SR <- JackStraw(SR, num.replicate = 100)
        SR <- ScoreJackStraw(SR, dims = 1:20)
        plot(ElbowPlot(SR))

        SR <- FindNeighbors(SR, dims = 1:10)
        SR <- FindClusters(SR, resolution = 0.5)

        SR <- RunUMAP(SR, dims = 1:15)
        plot(DimPlot(SR, reduction = "umap", label=T))

        SR.markers <- FindAllMarkers(SR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        
        top10 <- SR.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

    })
}


do_seurat_cellcycle = function(env){
    with(env, {
    
        SR <- CellCycleScoring(SR, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

        # view cell cycle scores and phase assignments
        #head(marrow[[]])
        plot(RidgePlot(SR, features = toupper(c("Pcna", "Top2a", "Mcm6", "Mki67")), ncol = 2))
        SR <- RunPCA(SR, features = c(s.genes, g2m.genes))
        plot(DimPlot(SR))

        SR <- ScaleData(SR, vars.to.regress = c("S.Score", "G2M.Score")) #, features = rownames(SR)
        SR <- RunPCA(SR, features = VariableFeatures(SR), nfeatures.print = 10)
        plot(DimPlot(SR))

        SR <- FindNeighbors(SR, dims = 1:10)
        SR <- FindClusters(SR, resolution = 0.5)

        SR <- RunUMAP(SR, dims = 1:15)
        plot(DimPlot(SR, reduction = "umap", label=T))
    })
    
}

do_seurat_cellcycle_remove=function(env){
    with(env, {
        SR[["percent.mt"]] <- PercentageFeatureSet(SR, pattern = "^MT-")
        VlnPlot(SR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

        plot1 <- FeatureScatter(SR, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(SR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        plot(CombinePlots(plots = list(plot1, plot2)))

        SR <- subset(SR, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25)

        SR <- NormalizeData(SR) #, normalization.method = "LogNormalize", scale.factor = 10000)

        SR <- FindVariableFeatures(SR, selection.method = "vst", nfeatures = 2000)

        # Identify the 10 most highly variable genes
        top10 <- head(VariableFeatures(SR), 20)

        # plot variable features with and without labels
        #VariableFeaturePlot(SR)
        plot1 <- VariableFeaturePlot(SR)
        plot(LabelPoints(plot = plot1, points = top10, repel = TRUE))
        #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        #CombinePlots(plots = list(plot1, plot2))

        sel.genes <- setdiff(rownames(SR), cell_cycle_genes)
        #SR <- ScaleData(SR, features = all.genes)
        
        SR <- ScaleData(SR, vars.to.regress = c("nFeature_RNA", "percent.mt"), features = sel.genes)
        #"nCount_RNA", 

        SR <- RunPCA(SR, features = VariableFeatures(object = SR))
        print(SR[["pca"]], dims = 1:5, nfeatures = 5)
        plot(VizDimLoadings(SR, dims = 1:2, reduction = "pca"))

        DimHeatmap(SR, dims = 1:10, cells = 500, balanced = TRUE)###

        SR <- JackStraw(SR, num.replicate = 100)
        SR <- ScoreJackStraw(SR, dims = 1:20)

        #!plot(JackStrawPlot(SR, dims = 1:15))

        plot(ElbowPlot(SR))

        SR <- FindNeighbors(SR, dims = 1:10)
        SR <- FindClusters(SR, resolution = 0.5)

        SR <- RunUMAP(SR, dims = 1:15)
        plot(DimPlot(SR, reduction = "umap", label=T))

        #head(SR@meta.data)
        #DimPlot(SR, reduction = "umap", label=T, group.by = "orig.ident")

        #saveRDS(SR, file="SR.both_30_39.UMAP.rds")
        #SR = readRDS("SR.both_30_39.UMAP.rds")

        SR.markers <- FindAllMarkers(SR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        library(dplyr)
        top10 <- SR.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
        plot(DoHeatmap(SR, features = top10$gene) + NoLegend())

        #DoHeatmap(SR, features = top10$gene) + NoLegend()
        #ggsave('heatmap_seurat.adrenal_30_39.pdf', width=10, height=30)

        #plot(FeaturePlot(SR, features = c("MKI67", "TH", "STAR", "DLK1", "PRRX1", "COL2A1")) )
    })
}

library(conos)
library(data.table)

get.scrublet.scores <- function(mat, min.molecules.per.gene=10) {
    # write out a CSV file
    tf <- tempfile()
    dtf <- paste(tf,'doubletScores',sep='.')
    dt <- data.table(as.matrix(t(mat[Matrix::rowSums(mat)>=min.molecules.per.gene,])))
    data.table::fwrite(dt,file=tf)
    cmd <- paste("python -c 'import sys; import pandas; import scrublet; df = pandas.read_csv(\"",tf,"\"); scrub = scrublet.Scrublet(df); doublet_scores, predicted_doublets = scrub.scrub_doublets(); pandas.DataFrame(doublet_scores).to_csv(\"",dtf,"\");'",sep='')
    tmp <- system(cmd, intern=T)
    #system(cmd);
    x <- as.data.frame(data.table::fread(dtf,sep=',',header=F,skip=1))
    x <- as.numeric(as.data.frame(data.table::fread(dtf,sep=',',header=F,skip=1))[,2])
    names(x)<-colnames(mat)
    file.remove(tf)
    file.remove(dtf)
    x
}


save_seurat_as_pagoda=function(SR, fout, app.title = 'adrenal_sr_10K'){
    library(pagoda2)
    library(igraph)

    p2 <- basicP2proc(SR@assays$RNA@counts, n.cores = 1)

    go.env <- p2.generate.human.go(p2)


    p2$clusters$PCA$seurat_cluster = as.factor(SR@meta.data$seurat_cluster)
    names(p2$clusters$PCA$seurat_cluster) = rownames(SR@meta.data)

    p2$embeddings$PCA$tSNE = as.matrix(SR@reductions$umap@cell.embeddings)
        #p2$embeddings$PCA = as.matrix(p2$embeddings$PCA@cell.embeddings)

    p2$clusters$PCA$timepoint = as.factor(SR@meta.data$orig.ident)
    names(p2$clusters$PCA$timepoint) = rownames(SR@meta.data)

    p2$embeddings$PCA$tSNE = as.matrix(SR@reductions$umap@cell.embeddings)
        #p2$embeddings$PCA = as.matrix(p2$embeddings$PCA@cell.embeddings)


    n.cores=1

    cat('Calculating hdea...\n')
    hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='multilevel',z.threshold=3, n.cores = n.cores)

    extraWebMetadata = NULL
    
    metadata.forweb <- list();
    metadata.forweb$timepoint <- p2.metadata.from.factor(p2$clusters$PCA$timepoint,displayname='timepoint')
    metadata.forweb$leiden <- p2.metadata.from.factor(p2$clusters$PCA$seurat_cluster,displayname='seurat_cluster')
    metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel,displayname='multilevel')
    metadata.forweb <- c(metadata.forweb, extraWebMetadata)
    genesets <- hierDiffToGenesets(hdea)
    appmetadata = list(apptitle=app.title)
    cat('Making KNN graph...\n')
    #p2$makeGeneKnnGraph(n.cores=n.cores)
    p2w = make.p2.app(p2, additionalMetadata = metadata.forweb, geneSets = genesets, dendrogramCellGroups = p2$clusters$PCA$multilevel, show.clusters=F, appmetadata = appmetadata)

    p2w$serializeToStaticFast(binary.filename = fout)    
}
