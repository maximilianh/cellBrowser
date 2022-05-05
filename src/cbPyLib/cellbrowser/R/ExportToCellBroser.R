findMarkers_new = function(object, field){
  
  library(scran)
  data = object@assays$RNA@data
  vec = object@meta.data[, field]
  
  a = findMarkers(as.matrix(data), vec, direction = "up", lfc = 1)
  
  ord = c('p_val', 'avg_log2FC', 'pct.1',	'pct.2', 
          'p_val_adj',	'cluster',	'gene')
  
  
  for(i in 1:length(a)){
    print(i)
    clusterInfo = as.data.frame(a[[i]])
    
    clusterInfo$avg_log2FC = apply(clusterInfo[, 5:ncol(clusterInfo)], 1, mean)
    
    clusterInfo$pct.1 = apply(clusterInfo[, 5:ncol(clusterInfo)], 1, max)
    clusterInfo$pct.2 = apply(clusterInfo[, 5:ncol(clusterInfo)], 1, min)
    
    clusterInfo = clusterInfo[, c('p.value', 'FDR', 'avg_log2FC', 
                                  'pct.1', 'pct.2')]
    colnames(clusterInfo) = c('p_val', 'p_val_adj', 'avg_log2FC', 'pct.1', 'pct.2')
    clusterInfo$cluster = names(a)[[i]]
    clusterInfo$gene = rownames(clusterInfo)
    clusterInfo = clusterInfo[, ord]
    
    clusterInfo = subset(clusterInfo, clusterInfo$p_val_adj < 0.01)
    
    if(nrow(clusterInfo) == 0){
      next
    }
    
    if(i == 1){
      markers = clusterInfo
    } else {
      markers = rbind.data.frame(markers, clusterInfo)
    }
  }
  return(markers)
}

exportData = function(object, dataName, outputDir, markerField, clusterField, use.scaleData){
  
  library(Seurat)
  library(scran)
  library(Matrix)
  library("R.utils")
  library(tibble)
  library(matrixStats)
  
  print(Sys.time())
  print('Reading Object')
  data = readRDS(object)
  print(Sys.time())
  
  #create output directory or check if it already exists
  if(dir.exists(outputDir)){
    print("Output directory found")
  } else {
    dir.create(outputDir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    print("Output directory created")
  }
  setwd(outputDir)
  
  #writing common barcodes and features files
  print(Sys.time())
  print('Writing barcodes')
  write.table(colnames(data@assays$RNA@counts), file = 'barcodes.tsv', 
            row.names = F, col.names = F, quote = F)
  print(Sys.time())
  
  print('Writing features')
  write.table(rownames(data@assays$RNA@counts), file = 'features.tsv', 
              row.names = F, col.names = F, quote = F)
  
  print(Sys.time())
  print('Writing count.mtx')
  #extract count data and write as mtx
  counts = Matrix(data@assays$RNA@counts, sparse = T)
  colnames(counts) = NULL
  rownames(counts) = NULL
  writeMM(counts, file = 'counts_exprMatrix.mtx')
  
  #gzip all the files
  print(Sys.time())
  print('zipping all the files')
  gzip('barcodes.tsv', overwrite = T, remove = T)
  gzip('features.tsv', overwrite = T, remove = T)
  gzip('counts_exprMatrix.mtx', overwrite = T, remove = T)
  
  #extract cpm data and write it out as mtx
  print(Sys.time())
  print('Writing cpm files')
  
  cpm = Matrix(data@assays$RNA@data, sparse = T)
  colnames(cpm) = NULL
  rownames(cpm) = NULL
  writeMM(cpm, file = 'data_exprMatrix.mtx')
  
  #gzip the cpm data file
  print(Sys.time())
  print('zip cpm file')
  gzip('data_exprMat.mtx', overwrite = T, remove = T)
  
  if(use.scaleData){

      print(Sys.time())
      print('Writing scale data file')
      cpm = data@assays$RNA@data
      zscore = (cpm - rowMeans(cpm))/(rowSds(cpm))[row(cpm)]

      zscore = Matrix(zscore, sparse = T)
      writeMM(zscore, file = 'scale_exprMatrix.mtx')
      gzip('scale_exprMatrix.mtx', overwrite = T, remove = T)
  }
  
  #print reduction coordinates
  len = length(names(data@reductions))
  red = names(data@reductions)
  
  for(reduction in 1:len){
    if(names(data@reductions)[reduction] %in% c('pca', 'harmony')){
      next
    }
    Umap = data@reductions[[red[reduction]]]
    umap = as.data.frame(Umap@cell.embeddings)
    
    umap = tibble::rownames_to_column(umap, "cellId")
    colnames(umap) = c('cellId', 'x', 'y')
    umap_path = paste0(red[reduction], ".coords.tsv")
    write.table(umap, file = umap_path, sep = "\t", quote = F,
                row.names = F)
  }
  
  if(length(markerField) > 1){
    print("Writing Markers file")
    for(i in 1:length(markerField)){
      markers = findMarkers_new(data, markerField[i])
      
      if(i == 1){
        combinedMarkers = markers
      } else {
        combinedMarkers = rbind.data.frame(combinedMarkers, markers)
      }
    }
  } else {
    combinedMarkers = findMarkers_new(data, markerField)
  }
  combinedMarkers = cbind.data.frame(rownames(combinedMarkers), combinedMarkers)
  colnames(combinedMarkers)[1] = 'geneID'
  write.table(combinedMarkers, file = 'markers.tsv', sep = '\t', 
              quote = F, row.names = F)
  
  meta = data@meta.data
  uniqueLen = apply(meta, 2, function(x) length(unique(x)))
  num = which(uniqueLen >= 32000)
  if(length(num) > 0){
    meta = meta[, -num]
  }
  meta = cbind.data.frame(rownames(meta), meta)
  colnames(meta)[1] = 'cellID'
  write.table(meta, file = 'meta.tsv', sep = '\t',
              quote = FALSE, row.names = F)
  
  ## Writing cellbrowser.conf
  
  config = '
name="%s"
shortLabel="%s"
exprMatrix="%s"
%s
meta="meta.tsv"
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
markers = %s
coords=%s
  '
  
  
  dataset.name = dataName
  matrixOutPath = sprintf("counts_exprMat.mtx.gz")
  
  if(use.scaleData){
  matrices.string = paste0("matrices=[{'label':'counts','fileName':'counts_exprMatrix.mtx.gz'},\n",
                            "{'label':'data','fileName':'data_exprMatrix.mtx.gz'},\n",
                            "{'label':'scale','fileName':'scale_exprMatrix.mtx.gz'}]" )
  } else {
    matrices.string = paste0("matrices=[{'label':'counts','fileName':'counts_exprMatrix.mtx.gz'},\n",
                             "{'label':'data','fileName':'data_exprMatrix.mtx.gz'}]")
  }
  cluster.field = clusterField
  
  
  enum.fields = colnames(meta)
  enum.string = paste0( "[", paste(paste0('"', enum.fields, '"'), collapse = ","), "]" )
  
  
  embeddings.conf = sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
                              red[-which(red %in% c("pca","harmony"))])
  coords.string = paste0( "[", paste(embeddings.conf, collapse = ",\n"), "]" )
  
  markers.string <- sprintf('markers = [{"file": "markers.tsv", "shortLabel": "Seurat Cluster Markers"}]')
  
  
  
  
  config <- sprintf(
    config,
    dataset.name,
    dataset.name,
    matrixOutPath,
    matrices.string,
    cluster.field,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  
  message("Writing cellbrowser config ")
  cat(config, file = 'cellbrowser.conf')
}

#Example Input
exportData(object = "../ExportToCellbrowser/rosmap_fib.rds",dataName = "ROSMAP_fib",
           outputDir = "../ExportToCellBrowser/check",
           markerField = c("Level_1","subtype"),clusterField = "Level_1",use.scaleData = FALSE)



