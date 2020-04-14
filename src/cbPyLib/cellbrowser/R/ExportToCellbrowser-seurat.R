# The markers "---" are used to mark the part that is used by the command line tool cbImportSeurat2
require(reticulate)

# ---
#' Write sparse matrix to .tsv.gz file by writing chunks, concating them with the Unix cat command,
#' then gziping the result. This does not work on Windows, we'd have to use the copy /b command there.
#'
#' @param inMat input matrix
#' @param outFname output file name, has to end with .gz
#' @param sliceSize=1000, size of each chunk in number of lines.
#' @examples
#' \dontrun{
#' writeSparseMatrix( pbmc_small@data, "exprMatrix.tsv.gz")
#' }
#'
writeSparseMatrix = function (inMat, outFname, sliceSize=1000) { 
    require(data.table)
    fnames = c(); 
    setDTthreads(8);  # otherwise this would use dozens of CPUs on a fat server
    mat = inMat; 
    geneCount = nrow(mat); 
    message("Writing expression matrix to ", outFname); 
    startIdx = 1; 
    while (startIdx<geneCount) { 
        endIdx=min(startIdx+sliceSize-1, geneCount); 
        matSlice = mat[startIdx:endIdx,]; 
        denseSlice = as.matrix(matSlice); 
        dt <- data.table(denseSlice); 
        dt = cbind(gene=rownames(matSlice), dt); 
        writeHeader = FALSE; 
        if (startIdx==1) { 
            writeHeader = TRUE 
        }; 
        sliceFname=paste0("temp", startIdx,".txt"); 
        fwrite(dt, sep="\t", file=sliceFname, quote = FALSE, col.names = writeHeader); 
        fnames=append(fnames, sliceFname); startIdx=startIdx+sliceSize}; 
        message("Concatenating chunks"); 
        system(paste("cat", paste(fnames, collapse=" "), "| gzip >", outFname, sep=" ")); 
        unlink(fnames); 
}

#' Return the correct matrix object from a Seurat object
#'
#' @param object Seurat object
#' @param matrix.slot the name of the slot
#'
findMatrix = function( object, matrix.slot ) {
      if (matrix.slot=="counts") {
          counts <- GetAssayData(object = object, slot="counts")
      } else if (matrix.slot=="scale.data") {
              counts <- GetAssayData(object = object, slot="scale.data")
      }
      else if (matrix.slot=="data") {
              counts <- GetAssayData(object = object)
      } else {
          error("matrix.slot can only be one of: counts, scale.data, data")
      }
}

#' Export Seurat object for UCSC cell browser
#'
#' @param object Seurat object
#' @param dir output directory path
#' @param dataset.name name of the dataset
#' @param meta.fields vector of metadata fields to export. By default all are exported.
#' @param meta.fields.names list that defines metadata field names
#'                          after the export. Should map metadata
#'                          column name to export field name
#' @param use.mtx for datasets > 100k cells, this is necessary, as otherwise R will report "problem too big"
#' @param matrix.slot one of "counts", "scale.data", "data". Defaults to "counts".
#' @param embeddings vector of embedding names to export. By default all are exported.
#' @param markers.file path to file with marker genes
#' @param cluster.field name of the metadata field containing cell cluster
#' @param cb.dir in which dir to create UCSC cell browser content root
#' @param port \emph{experimental} on which port to run UCSC cell browser after export
#'
#' @export
#'
#' @importFrom reticulate py_module_available
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#' ExportToCellbrowser(object = pbmc_small, dataset.name = "PBMC", dir = "out")
#' }
#'
ExportToCellbrowser <- function(
  object,
  dir,
  dataset.name,
  meta.fields = NULL,
  meta.fields.names = NULL,
  embeddings = NULL,
  matrix.slot = "counts",
  markers.file = NULL,
  cluster.field = NULL,
  port = NULL,
  cb.dir = NULL,
  markers.n = 100,
  skip.expr.matrix = FALSE,
  skip.markers = FALSE,
  all.meta = FALSE,
  use.mtx = FALSE
) {
  if (!require("Seurat",character.only = TRUE)) {
          stop("This script requires that Seurat (V2 or V3) is installed")
  }
  message("Seurat Version installed: ", packageVersion("Seurat"))
  message("Object was created with Seurat version ", object@version)

  if (substr(object@version, 1, 1)!='2' && substr(object@version, 1, 1)!='3') {
          stop("can only process Seurat2 or Seurat3 objects, version of rds is ", object@version)
  }

  if (substr(object@version, 1, 1)!=substr( packageVersion("Seurat"), 1, 1)) {
          stop("The installed version of Seurat is different from Seurat object loaded. You have to down- or upgrade your installed Seurat version, see the Seurat documentation")
  }


  # compatibility layer for Seurat 2 vs 3 
  # see https://satijalab.org/seurat/essential_commands.html
  if (substr(object@version, 1, 1)=='2') {
      # Seurat 2 data access
      idents <- object@ident # Idents() in Seurat3
      meta <- object@meta.data
      cellOrder <- object@cell.names
      if (matrix.slot=="counts") {
          counts <- object@raw.data
      } else if (matrix.slot=="scale.data") {
          counts <- object@scale.data
      }
      else if (matrix.slot=="data") {
          counts <- object@data
      } else {
          error("matrix.slot can only be one of: counts, scale.data, data")
      }
      genes <- rownames(x = object@data)
      dr <- object@dr
  } else {
      # Seurat 3 functions
      idents <- Idents(object)
      meta <- object@meta.data
      cellOrder <- colnames(object)
      counts <- findMatrix(object, matrix.slot)
      if (dim(counts)[1]==0) { 
          message(paste0("The Seurat data slot '", matrix.slot, "' contains no data. Trying default assay."))
          defAssay = DefaultAssay(object)
          assay = GetAssay(object, defAssay)
          message(paste0("Default assay is ", defAssay))
          counts <- findMatrix(assay, matrix.slot)
          genes <- rownames(counts)

          if (dim(counts)[1]==0) { 
              stop("Could not find an expression matrix",
                     "Please select the correct slot where the matrix is stored, possible ",
                     "values are 'counts', 'scale.data' or 'data'. To select a slot, ",
                   "use the option 'matrix.slot' from R or the cbImportSeurat option -s from the command line.")
          }
      }
      else {
          genes <- rownames(x = object)
      }
      dr <- object@reductions
  }

  if (is.null(cluster.field)) {
          cluster.field = "Cluster"
  }

  if (is.null(meta.fields)) {
    meta.fields <- colnames(meta)
    if (length(levels(idents)) > 1) {
      meta.fields <- c(meta.fields, ".ident")
    }
  }

  if (!is.null(port) && is.null(cb.dir)) {
    stop("cb.dir parameter is needed when port is set")
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  if (!dir.exists(dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  }

  enum.fields <- c()

  # Export expression matrix

  if (!skip.expr.matrix) { 
      if (use.mtx) {
            require(Matrix)
            require(R.utils)
            matrixPath <- file.path(dir, "matrix.mtx")
            genesPath <- file.path(dir, "features.tsv")
            barcodesPath <- file.path(dir, "barcodes.tsv")
            message("Writing expression matrix to ", matrixPath)
            writeMM(counts, matrixPath)
            # easier to load if the genes file has at least columns even though seurat objects
            # don't have yet explicit geneIds/geneSyms data, we just duplicate whatever the matrix has now
            write.table(as.data.frame(cbind(rownames(counts), rownames(counts))), file=genesPath, sep="\t", row.names=F, col.names=F, quote=F)
            write(colnames(counts), file = barcodesPath)
            message("Gzipping expression matrix")
            gzip(matrixPath)
            gzip(genesPath)
            gzip(barcodesPath)
      } else {
          gzPath <- file.path(dir, "exprMatrix.tsv.gz")
          if ((((ncol(counts)/1000)*(nrow(counts)/1000))>2000) && is(counts, 'sparseMatrix')) {
              writeSparseMatrix(counts, gzPath);
          } else {
              mat = as.matrix(counts)
              df <- as.data.frame(mat, check.names=FALSE)
              df <- data.frame(gene=genes, df, check.names = FALSE)
              z <- gzfile(gzPath, "w")
              message("Writing expression matrix to ", gzPath)
              write.table(x = df, sep="\t", file=z, quote = FALSE, row.names = FALSE)
              close(con = z)
          }
      }
  }

  # Export cell embeddings
  if (is.null(embeddings)) {
      embeddings = names(dr)
      message("Using all embeddings contained in the Seurat object: ", embeddings)
  }

  embeddings.conf <- c()
  for (embedding in embeddings) {
    emb <- dr[[embedding]]
    if (is.null(emb)) {
        message("Embedding ",embedding," does not exist in Seurat object. Skipping. ")
        next
    }
    #df <- slot(emb, "cell.embeddings")
    df <-  emb@cell.embeddings
    if (ncol(df) > 2) {
      warning('Embedding ', embedding, ' has more than 2 coordinates, taking only the first 2')
      df <- df[, 1:2]
    }
    colnames(df) <- c("x", "y")
    df <- data.frame(cellId = rownames(df), df, check.names=FALSE)
    fname <- file.path(
      dir,
      sprintf("%s.coords.tsv", embedding)
    )
    message("Writing embeddings to ", fname)
    write.table(df[cellOrder, ], sep="\t", file=fname, quote = FALSE, row.names = FALSE)
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      embedding
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }

  # Export metadata
  df <- data.frame(row.names = cellOrder, check.names=FALSE)
  for (field in meta.fields) {
    if (field == ".ident") {
      df$Cluster <- idents
      enum.fields <- c(enum.fields, "Cluster")
    } else {
      name <- meta.fields.names[[field]]
      if (is.null(name)) {
        name <- field
      }
      #df[[name]] <- object[[field]][, 1]
      #df[[name]] <- FetchData(object, field)[, 1]
      df[[name]] <- meta[[field]]
      if (!is.numeric(df[[name]])) {
        enum.fields <- c(enum.fields, name)
      }
    }
  }
  df <- data.frame(Cell=rownames(df), df, check.names=FALSE)

  fname <- file.path(dir, "meta.tsv")
  message("Writing meta data to ", fname)
  write.table(df[cellOrder, ], sep="\t", file=fname, quote = FALSE, row.names=FALSE)

  # Export markers
  markers.string <- ''
  if (is.null(markers.file)) {
    ext <- "tsv"
  } else {
    ext <- tools::file_ext(markers.file)
  }
  file <- paste0("markers.", ext)
  fname <- file.path(dir, file)

  if (!is.null(markers.file) && !skip.markers) {
    message("Copying ", markers.file, " to ", fname)
    file.copy(markers.file, fname)
  }
  if (is.null(markers.file) && skip.markers) {
    file <- NULL
  }

  if (is.null(markers.file) && !skip.markers) {
    if (length(levels(idents)) > 1) {

          markers.helper <- function(x) {
            partition <- markers[x,]
            ord <- order(partition$p_val_adj < 0.05, -partition$avg_logFC)
            res <- x[ord]
            naCount <- max(0, length(x) - markers.n)
            res <- c(res[1:markers.n], rep(NA, naCount))
            return(res)
          }

      if (.hasSlot(object, "misc") && !is.null(object@misc["markers"][[1]])) {
          message("Found precomputed markers in obj@misc['markers']")
          markers <- object@misc["markers"]$markers
      } else {
          message("Running FindAllMarkers(), using wilcox test, min logfc diff 0.25")
          markers <- FindAllMarkers(object, do.print=TRUE, print.bar=TRUE, test.use="wilcox", logfc.threshold = 0.25)
      }

      message("Writing top ", markers.n, ", cluster markers to ", fname)
      markers.order <- ave(rownames(markers), markers$cluster, FUN=markers.helper)
      top.markers <- markers[markers.order[!is.na(markers.order)],]
      write.table(top.markers, fname, quote=FALSE, sep="\t", col.names=NA)

    } else {

      message("No clusters found in Seurat object and no external marker file provided, so no marker genes can be computed")
      file <- NULL
    }
  }

  if (!is.null(file)) {
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      file
    )
  }

  matrixOutPath <- "exprMatrix.tsv.gz"
  if (use.mtx)
      matrixOutPath <- "matrix.mtx.gz"

  config <- '
# This is a bare-bones cellbrowser config file auto-generated from R.
# Look at https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf
# for a full file that shows all possible options
name="%s"
shortLabel="%1$s"
exprMatrix="%s"
#tags = ["10x", "smartseq2"]
meta="meta.tsv"
# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"
geneIdType="auto"
# file with gene,description (one per line) with highlighted genes, called "Dataset Genes" in the user interface
# quickGenesFile="quickGenes.csv"
clusterField="%s"
labelField="%s"
enumFields=%s
%s
coords=%s'

  enum.string <- paste0(
    "[",
    paste(paste0('"', enum.fields, '"'), collapse = ", "),
    "]"
  )
  coords.string <- paste0(
    "[",
    paste(embeddings.conf, collapse = ",\n"),
    "]"
  )
  config <- sprintf(
    config,
    dataset.name,
    matrixOutPath,
    cluster.field,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )

  confPath = file.path(dir, "cellbrowser.conf")
  message("Writing cellbrowser config to ", confPath)
  cat(config, file=confPath)

  message("Prepared cellbrowser directory ", dir)

  if (!is.null(cb.dir)) {
    if (!py_module_available("cellbrowser")) {
          stop(
            "The Python package `cellbrowser` is required to prepare and run ",
            "Cellbrowser. Please install it ",
            "on the Unix command line with `sudo pip install cellbrowser` (if root) ",
            "or `pip install cellbrowser --user` (as a non-root user). ",
            "To adapt the Python that is used, you can either set the env. variable RETICULATE_PYTHON ",
            "or do `require(reticulate) and use one of these functions: use_python(), use_virtualenv(), use_condaenv(). ",
            "See https://rstudio.github.io/reticulate/articles/versions.html . ",
            "At the moment, R's reticulate is using this Python: ",import('sys')$executable,". "
           )
    }

    if (!is.null(port)) {
      port <- as.integer(port)
    }
    message("Converting cellbrowser directory to html/json files")
    cb <- import("cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
    message("HTML files are ready in ", cb.dir)

    if (!is.null(port)) {
      message("Starting http server")
      cb$cellbrowser$stop()
      cb$cellbrowser$serve(cb.dir, port)
      Sys.sleep(0.4)
      utils::browseURL(sprintf("http://localhost:%d", port))
    }
  }
}
# ---

#' Stop Cellbrowser web server
#'
#' @export
#'
#' @importFrom reticulate py_module_available
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#' StopCellbrowser()
#' }
#'
StopCellbrowser <- function() {
  if (py_module_available("cellbrowser")) {
    cb <- import("cellbrowser")
    cb$cellbrowser$stop()
  } else {
    stop("The `cellbrowser` package is not available in the Python used by R's reticulate")
  }
}
