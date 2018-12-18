# the following is copied from Seurat3's utilities.R and generics.R
# it was modified to work with Seurat2. At some point we probably should merge the two code bases.
# the markers "---" are used to mark the part that is used by the command line tool cbImportSeurat2
require(reticulate)

# ---
#' Export Seurat object for UCSC cell browser
#'
#' @param object Seurat object
#' @param dir output directory path
#' @param dataset.name name of the dataset
#' @param meta.fields vector of metadata fields to export
#' @param meta.fields.names list that defines metadata field names
#'                          after the export. Should map metadata
#'                          column name to export name
#' @param embeddings vector of embedding names to export
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
  embeddings = c("tsne"),
  markers.file = NULL,
  cluster.field = NULL,
  port = NULL,
  cb.dir = NULL,
  markers.n = 100,
  skip.expr.matrix = FALSE,
  skip.markers = FALSE
) {
  idents <- FetchData(object,c("ident"))$ident;

  if (is.null(meta.fields)) {
    meta.fields <- c("nGene", "nUMI")
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

  order <- object@cell.names
  enum.fields <- c()

  # Export expression matrix
  if (!skip.expr.matrix) {
      df <- as.data.frame(as.matrix(GetAssayData(object = object)), check.names=FALSE)
      df <- data.frame(gene=rownames(object@data), df, check.names = FALSE)
      gzPath <- file.path(dir, "exprMatrix.tsv.gz")
      z <- gzfile(gzPath, "w")
      message("Writing expression matrix to ", gzPath)
      write.table(x = df, sep="\t", file=z, quote = FALSE, row.names = FALSE)
      close(con = z)
  }

  # Export cell embeddings
  embeddings.conf <- c()
  dr <- object@dr
  for (embedding in embeddings) {
    #df <- Embeddings(object = object, reduction = embedding)
    emb <- dr[[embedding]]
    df <- slot(emb, "cell.embeddings")
    if (ncol(df) > 2) {
      warning(
        'Embedding ', embedding,
        ' has more than 2 coordinates, taking only the first 2',
      )
      df <- df[, 1:2]
    }
    colnames(df) <- c("x", "y")
    df <- data.frame(cellId = rownames(df), df, check.names=FALSE)
    fname <- file.path(
      dir,
      sprintf("%s.coords.tsv", embedding)
    )
    message("Writing embeddings to ", fname)
    write.table(df[order, ], sep="\t", file=fname, quote = FALSE, row.names = FALSE)
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      embedding
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }

  # Export metadata
  df <- data.frame(row.names = object@cell.names, check.names=FALSE)
  for (field in meta.fields) {
    if (field == ".ident") {
      df$Cluster <- idents # no Idents() in Seurat2
      enum.fields <- c(enum.fields, "Cluster")
      cluster.field <- "Cluster"
    } else {
      name <- meta.fields.names[[field]]
      if (is.null(name)) {
        name <- field
      }
      #df[[name]] <- object[[field]][, 1]
      df[[name]] <- FetchData(object, field)[, 1]
      if (!is.numeric(df[[name]])) {
        enum.fields <- c(enum.fields, name)
      }
    }
  }
  df <- data.frame(Cell=rownames(df), df, check.names=FALSE)
  fname <- file.path(dir, "meta.tsv")
  message("Writing meta data to ", fname)
  write.table(df[order, ], sep="\t", file=fname, quote = FALSE, row.names=FALSE)

  # Export markers
  markers.string <- ''
  if (is.null(markers.file)) {
    ext <- "tsv"
  } else {
    ext <- tools::file_ext(markers.file)
  }
  file <- paste0("markers.", ext)
  fname <- file.path(dir, file)
  if (is.null(markers.file) & !skip.markers) {
    if (length(levels(idents)) > 1) {
      message("Running FindAllMarkers() and writing cluster markers to ", fname)
      markers <- FindAllMarkers(object, do.print=TRUE, print.bar=TRUE)
      top.markers <- markers %>% group_by(cluster) %>% top_n(markers.n, avg_logFC)
      write.table(top.markers, fname, quote=FALSE, sep="\t", row.names=FALSE)
    } else {
      file <- NULL
    }
  } else if (!skip.markers) {
    file.copy(markers.file, fname)
  }
  if (!is.null(file)) {
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      file
    )
  }

  config <- 'name="%s"
shortLabel="%1$s"
exprMatrix="exprMatrix.tsv.gz"
#tags = ["10x", "smartseq2"]
meta="meta.tsv"
# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"
geneIdType="auto"
clusterField="%s"
labelField="%2$s"
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
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  cat(config, file=file.path(dir, "cellbrowser.conf"))
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
