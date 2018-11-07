#' clusterProfiler.R
#' 
#' @description clusterProfiler
#' 
#' run clusterProfiler
#' 
#' @param GTF character. GTF file
#' @param inPath character vector. input
#' @param outPath character vector. output
#' @usage postAna::clusterProfiler(inPath, outPath, ...)
#' @return NULL
#' @details TBC
#' @note TBC
#' @section Main Warning:
#' TBC
#' \subsection{Other Warning} {
#'   TBC
#' }
#' 
#' @references TBC
#' @example
#' postAna::clusterProfiler("overlap.txt", "enrichGO.ALL.txt")



rm(list = ls())


# IMPORT logging library

# https://github.com/zatonovo/futile.logger

if(! "futile.logger" %in% rownames(installed.packages())) {
  if(! "devtools" %in% rownames(installed.packages())) {
    install.packages(
      "devtools",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  
  library("devtools")
  
  devtools::install_github(
    repo = "zatonovo/futile.logger"
  )
}
library("futile.logger") # help(package = "futile.logger") # ls(pos = "package:futile.logger")


# (v.log4j <- c(
#   OFF = NULL,
#   FATAL = futile.logger::FATAL,
#   ERROR = futile.logger::ERROR,
#   WARN = futile.logger::WARN,
#   INFO = futile.logger::INFO,
#   DEBUG = futile.logger::DEBUG,
#   TRACE = futile.logger::TRACE,
#   ALL = NULL
# ))


flog.info("PARSE arguments")

# https://docs.python.org/3/howto/argparse.html

if(! "argparse" %in% rownames(installed.packages())) {
  library("devtools")
  
  devtools::install_github(
    repo = "trevorld/argparse"
  )
}
library("argparse") # help(package = "argparse") # ls(pos = "package:argparse")

v.desc <- c(
  "clusterProfiler.R",
  "TBC",
  "0.0.1",
  "LGPLv3",
  "Wen-Ching Chan <wchan10@bsd.uchicago.edu>"
)

names(v.desc) <- c(
  "Program",
  "Description",
  "Version",
  "License",
  "Contact"
)

# https://docs.python.org/3/library/argparse.html#the-add-argument-method

parser <- argparse::ArgumentParser(
  description = paste(
    paste(
      names(v.desc),
      v.desc,
      sep = ": "
    ),
    collapse = "\\n"
  )
)

parser$add_argument(
  '--inPath',
  dest = "c.v.in.file.path",
  type = "character",
  required = TRUE,
  help = "input (default: required arg)",
  metavar = "IN"
)

parser$add_argument(
  '--compPair',
  dest = "c.v.comp.pair",
  type = "character",
  required = TRUE,
  help = "input (default: required arg)",
  metavar = "COMPARISON_PAIR"
)

parser$add_argument(
  '--outPath',
  dest = "c.out.file.path",
  type = "character",
  required = TRUE,
  help = "output (default: required arg)",
  metavar = "OUT"
)

parser$add_argument(
  '--assembly',
  dest = "c.assembly",
  type = "character",
  required = FALSE,
  help = "ASSEMBLY (default: NULL)",
  metavar = "ANCHOR"
)

parser$add_argument(
  '--logPath',
  dest = "c.log.file.path",
  type = "character",
  required = FALSE,
  help = "Log file (default: NULL)",
  metavar = "LOG"
)

args <- try(
  parser$parse_args(
    commandArgs(TRUE)
  ),
  silent = TRUE
)


if(methods::is(args, "try-error")) {
  geterrmessage()
  parser$print_help()
  cat("\n")
  print(sessionInfo())
  q(save = "no")
} else if(methods::is(args, "list")) {
  dput(args)
}


str2bool <- function(s) {
  if(any(tolower(s) %in% c('yes', 'true', 't', 'y', '1')))
    return(TRUE)
  else if(any(tolower(s) %in% c('no', 'false', 'f', 'n', '0')))
    return(FALSE)
  else
    stop(paste("ArgumentTypeError('Boolean value expected.'), but", s))
}

(args$c.out.dir.path <- dirname(args$c.out.file.path))
if(! file.exists(args$c.out.dir.path)) dir.create(args$c.out.dir.path, showWarnings = TRUE, recursive = TRUE)


if(is.null(args$c.log.file.path)) {
  args$c.log.file.path <- file.path(
    args$c.out.dir.path,
    paste(
      format(
        Sys.time(),
        '%Y-%m-%d_%H-%M-%S'
      ),
      "log",
      sep = "."
    )
  )
}

flog.logger(
  name <- "ROOT",
  threshold = INFO,
  appender = appender.tee(
    args$c.log.file.path
  )
)

for(c.name in names(args)) {
  flog.info(paste("VAR:", c.name, args[[c.name]], sep = "\t"))
}

flog.debug("PARSE arguments - DONE")


flog.info("PREPROCESS")


flog.info(paste("DETERMINE Species:\t", c.species <- "Homo sapiens"))
flog.info(paste("DETERMINE AnnotHub name of OrgDb:\t", c.AnnotHub.org.db.name <- 'AH61777'))
flog.info(paste("DETERMINE OrgDb:\t", c.org.db <- 'org.Hs.eg.db'))

if(! c.org.db %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install(c.org.db, version = "devel")
}

flog.info(paste("DETERMINE AnnotHub name of TxDb:\t", c.AnnotHub.txdb.name <- 'AH52260'))
flog.info(paste("DETERMINE TxDb:\t", c.txdb <- 'TxDb.Hsapiens.UCSC.hg38.knownGene'))

if(! c.txdb %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install(c.txdb, version = "devel")
}

if(! "TxDb.Hsapiens.UCSC.hg19.knownGene" %in% rownames(installed.packages())) {
  if(! "BiocUpgrade" %in% rownames(installed.packages())) {
    if(! "BiocManager" %in% rownames(installed.packages()))
      install.packages(
        "BiocManager",
        dependencies = TRUE,
        repos = "https://cloud.r-project.org"
      )
    BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "devel")
  }
  
  if(! "devtools" %in% rownames(installed.packages())) {
    install.packages(
      "devtools",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  library(devtools)
  
  install_github(
    "GuangchuangYu/ChIPseeker",
    build_vignettes = FALSE,
    repos = BiocInstaller::biocinstallRepos(),
    dependencies = TRUE
  )
}
library("ChIPseeker")

clusterProfiler::idType(OrgDb = c.org.db)



if(! "KEGGREST" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("KEGGREST", version = "devel")
}
library("KEGGREST")

flog.info(paste("DETERMINE KEGG organism:\t", c.org.kegg <- 'hsa'))


v.in.file.path <- unlist(strsplit(x = args$c.v.in.file.path, split = ","))

if(length(v.in.file.path) == 1) {
  names(v.in.file.path) <- args$c.v.comp.pair
}

lst.gene.id <- list()

for(i.idx in seq_len(length(v.in.file.path))) {
  # (i.idx <- 1)
  c.in.file.path <- v.in.file.path[i.idx]
  c.label <- names(v.in.file.path)[i.idx]
  flog.debug(paste("PROCESS", c.in.file.path, "with", c.label))
  
  df.in <- read.csv(
    file = c.in.file.path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  c.idx.gene <- 1
  
  df.in.uniq <- df.in[! duplicated(df.in[, c.idx.gene]), ]
  flog.debug(paste("REMOVE duplicate KEY from", nrow(df.in), "to", nrow(df.in.uniq)))
  v.from <- df.in.uniq[, c.idx.gene, drop = TRUE]
  
  c.key <- "ENSEMBL"
  v.from <- gsub("\\.\\d+$", "", v.from)
  
  df.key.2.feat <- try(clusterProfiler::bitr(geneID = v.from, fromType = c.key, toType = c("ENTREZID", "SYMBOL", "GENENAME"), OrgDb = c.org.db, drop = TRUE))
  
  i.idx.log2FC <- which(colnames(df.in.uniq) == "log2FoldChange")
  df.symb.2.log2FC <- data.frame("KEY" = v.from, "log2FC" = as.numeric(df.in.uniq[, i.idx.log2FC]))
  
  df.t <- merge(df.key.2.feat, df.symb.2.log2FC, by.x = c.key, by.y = "KEY")
  
  v.t <- df.t$log2FC
  names(v.t) <- df.t$ENTREZID
  # <
  
  lst.gene.id[[c.label]] <- v.t
  
  flog.info(paste(length(lst.gene.id[[c.label]]), "entries in", c.label))
}

c.keytype <- "ENTREZID"

flog.debug("PREPROCESS - DONE")


v.ont <- c(
  "BP" = "Biological Process",
  "CC" = "Cellular Component",
  "MF" = "Molecular Function",
  "ALL" = "Pooled GO"
)

for(i.idx in seq_len(length(lst.gene.id))) {
  # (i.idx <- 1)
  c.name <- names(lst.gene.id)[i.idx]
  flog.debug(paste("PROCESS", c.name))
  
  c.out.dir.path <- args$c.out.dir.path
  
  c.ont <- "ALL"
  
  flog.info(paste("GO Classification in", v.ont[c.ont], "(", names(v.ont[c.ont]), ")"))
  
  flog.info(paste("GO Over-representation Test in", v.ont[c.ont]))
  
  system.time(
    enrichRes.go <- try(
      clusterProfiler::enrichGO(
        gene = names(lst.gene.id[[c.name]]), 
        OrgDb = c.org.db,
        keyType = c.keytype, 
        ont = c.ont,
        readable = T,
        pool = (c.ont == "ALL")
      )
    )
  )
  
  
  (c.out.file.t.path <- args$c.out.file.path)
  flog.info(paste("SAVE GO enrichment list in", v.ont[c.ont], "to", c.out.file.t.path))
  
  write.table(
    x = as.data.frame(enrichRes.go),
    file = c.out.file.t.path,
    quote = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  c.out.file.t.path <- gsub("\\.txt$", ".dotplot.pdf", args$c.out.file.path)
  flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
  print(
    enrichplot::dotplot(
      object = enrichRes.go,
      font.size = 12, 
      title = paste(c.name, c.ont, "GO")
    )
  )
  dev.off()
  
  c.out.file.t.path <- gsub("\\.txt$", ".emapplot.pdf", args$c.out.file.path)
  flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 8, onefile = TRUE)
  print(
    enrichplot::emapplot(
      x = enrichRes.go,
      showCategory = 10
    )
  )
  dev.off()
  
  c.out.file.t.path <- gsub("\\.txt$", ".cnetplot.pdf", args$c.out.file.path)
  flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 8, onefile = TRUE)
  print(
    enrichplot::cnetplot(
      x = enrichRes.go,
      showCategory = 5, 
      foldChange = lst.gene.id[[i.idx]], 
      layout = "kk"
    )
  )
  dev.off()
  
  flog.debug(paste("GO Over-representation Test in", v.ont[c.ont], "- DONE"))
  
  
  flog.info(paste("GO Gene Set Enrichment Analysis in", v.ont[c.ont]))
  
  system.time(
    gseaRes.go <- try(
      clusterProfiler::gseGO(
        geneList = sort(lst.gene.id[[c.name]], decreasing = TRUE),
        ont = as.character(c.ont),
        OrgDb = c.org.db,
        keyType = c.keytype,
        verbose = T
      )
    )
  )
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", ".enrichGSEAGO.ALL.txt", args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  write.table(
    x = as.data.frame(gseaRes.go),
    file = c.out.file.t.path,
    quote = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  c.geneSetID <- head(subset(gseaRes.go@result, NES > 0)$ID, min(5, nrow(subset(gseaRes.go@result, NES > 0))))[1]
  c.top <- paste("pos", sprintf(fmt = "%03d", match(c.geneSetID, subset(gseaRes.go@result, NES > 0)$ID)), sep = "")
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", paste0(".enrichGSEAGO.ALL.", c.top, ".pdf"), args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
  print(
    enrichplot::gseaplot(
      x = gseaRes.go,
      geneSetID = c.geneSetID,
      by = "all",
      title = c.geneSetID,
      color.line = "red",
      color.vline = "blue"
    )
  )
  dev.off()
  
  
  c.geneSetID <- head(subset(gseaRes.go@result, NES < 0)$ID, min(5, nrow(subset(gseaRes.go@result, NES < 0))))[1]
  c.top <- paste("neg", sprintf(fmt = "%03d", match(c.geneSetID, subset(gseaRes.go@result, NES < 0)$ID)), sep = "")
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", paste0(".enrichGSEAGO.ALL.", c.top, ".pdf"), args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
  print(
    enrichplot::gseaplot(
      x = gseaRes.go,
      geneSetID = c.geneSetID,
      by = "all",
      title = c.geneSetID,
      color.line = "green",
      color.vline = "blue"
    )
  )
  dev.off()
  
  flog.debug(paste("GO Gene Set Enrichment Analysis in", v.ont[c.ont], "- DONE"))
  
  
  flog.info(paste("KEGG over-representation test"))
  
  system.time(
    enrichRes.kegg <- try(
      clusterProfiler::enrichKEGG(
        gene = names(lst.gene.id[[c.name]]), 
        organism = c.org.kegg, 
        keyType = "kegg"
      )
    )
  )
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", ".enrichKEGG.txt", args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  write.table(
    x = as.data.frame(enrichRes.kegg),
    file = c.out.file.t.path,
    quote = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", ".enrichKEGG.dotplot.pdf", args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
  print(
    enrichplot::dotplot(
      object = enrichRes.kegg, 
      font.size = 12, 
      title = paste(c.name, "KEGG")
    )
  )
  dev.off()
  
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", ".enrichKEGG.emapplot.pdf", args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
  print(
    enrichplot::emapplot(
      x = enrichRes.kegg
    )
  )
  dev.off()
  
  if(identical(enrichRes.kegg@gene, enrichRes.go@gene)) {
    enrichRes.kegg@gene2Symbol <- enrichRes.go@gene2Symbol
    enrichRes.kegg@result$SYMBOL <- sapply(
      enrichRes.kegg@result$geneID,
      function(x) {
        paste(
          enrichRes.kegg@gene2Symbol[as.character(unlist(strsplit(x, "/", fixed = TRUE)))],
          collapse = "/"
        )
      }
    )
    enrichRes.kegg@readable <- TRUE
  }
  
  enrichRes.kegg.t <- enrichRes.kegg
  enrichRes.kegg.t@result$geneID <- enrichRes.kegg.t@result$SYMBOL

  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", ".enrichKEGG.cnetplot.pdf", args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

  pdf(file = c.out.file.t.path, width = 16, height = 16, onefile = F)
  print(
    enrichplot::cnetplot(
      x = enrichRes.kegg.t,
      showCategory = 5,
      foldChange = lst.gene.id[[i.idx]],
      layout = "kk"
    )
  )
  dev.off()
  
  flog.debug(paste("KEGG over-representation test - DONE"))
  
  
  flog.info(paste("KEGG Gene Set Enrichment Analysis in"))
  
  system.time(
    gseaRes.kegg <- try(
      clusterProfiler::gseKEGG(
        geneList = sort(lst.gene.id[[c.name]], decreasing = TRUE),
        organism = c.org.kegg,
        keyType = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot")[1],
        use_internal_data = FALSE,
        verbose = T
      )
    )
  )
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", ".enrichGSEAKEGG.txt", args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  write.table(
    x = as.data.frame(gseaRes.kegg),
    file = c.out.file.t.path,
    quote = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  c.geneSetID <- head(subset(gseaRes.kegg@result, NES > 0)$ID, min(5, nrow(subset(gseaRes.kegg@result, NES > 0))))[1]
  c.top <- paste("pos", sprintf(fmt = "%03d", match(c.geneSetID, subset(gseaRes.kegg@result, NES > 0)$ID)), sep = "")
  
  c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", paste0(".enrichGSEAKEGG.", c.top, ".pdf"), args$c.out.file.path)
  flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
  print(
    enrichplot::gseaplot(
      x = gseaRes.kegg,
      geneSetID = c.geneSetID,
      by = "all",
      title = c.geneSetID,
      color.line = "red",
      color.vline = "blue"
    )
  )
  dev.off()
  
  
  c.geneSetID <- head(subset(gseaRes.kegg@result, NES < 0)$ID, min(5, nrow(subset(gseaRes.kegg@result, NES < 0))))[1]
  if(! is.na(c.geneSetID)) {
    c.top <- paste("neg", sprintf(fmt = "%03d", match(c.geneSetID, subset(gseaRes.kegg@result, NES < 0)$ID)), sep = "")
    
    c.out.file.t.path <- gsub("\\.enrichGO.ALL.txt$", paste0(".enrichGSEAKEGG.", c.top, ".pdf"), args$c.out.file.path)
    flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
    
    pdf(file = c.out.file.t.path, width = 8, height = 6, onefile = F)
    print(
      enrichplot::gseaplot(
        x = gseaRes.kegg,
        geneSetID = c.geneSetID,
        by = "all",
        title = c.geneSetID,
        color.line = "green",
        color.vline = "blue"
      )
    )
    dev.off()
  }
  
  flog.debug(paste("GO Gene Set Enrichment Analysis in", v.ont[c.ont], "- DONE"))
}


print(sessionInfo())
q(save = "no")
