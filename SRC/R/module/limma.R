#' limma.R
#' 
#' @description limma
#' 
#' run limma
#' 
#' @param metadata character. meta data file
#' @param GTF character. GTF file
#' @param inPath character vector. input
#' @param outPath character vector. output
#' @usage callLoci::limma(metadata, GTF, inPath, ...)
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
#' callLoci::limma("metadata.txt", "Ensembl.gtf", "input_dir", "output_dir")



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
  "limma.R",
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
  '--metadata',
  dest = "c.metadata.path",
  type = "character", # c('logical', 'integer', 'double', 'character')
  required = TRUE,
  help = "meta data file (default: required arg)",
  metavar = "METADATA"
)

parser$add_argument(
  '--GTF',
  dest = "c.gtf.file.path",
  type = "character",
  required = TRUE,
  help = "GTF (default: required arg)",
  metavar = "GTF"
)

parser$add_argument(
  '--inPath',
  dest = "c.in.dir.path",
  type = "character",
  required = TRUE,
  help = "input (default: required arg)",
  metavar = "IN"
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
  '--filterLowExpr',
  dest = "b.filter.low.expr",
  default = "TRUE",
  type = "character",
  required = FALSE,
  help = "filtering for lowly expressed genes",
  metavar = "RM_LOW_EXPR"
)

parser$add_argument(
  '--criteria',
  dest = "c.criteria",
  default = "fc_q",
  type = "character",
  required = FALSE,
  help = "Criteria c('fc_q', 'fc', 'fc_p', 'p') (default: fc_q)",
  metavar = "CRITERIA"
)

parser$add_argument(
  '--foldChange',
  dest = "d.fold.change",
  default = 1.5,
  type = "double",
  required = FALSE,
  help = "Fold Change (default: 1.5)",
  metavar = "FC"
)

parser$add_argument(
  '--pval',
  dest = "d.pval",
  default = 0.1,
  type = "double",
  required = FALSE,
  help = "p-value (default: 0.1)",
  metavar = "P"
)

parser$add_argument(
  '--Case2CtrlPairs',
  dest = "c.v.case.2.ctrl",
  type = "character",
  required = FALSE,
  help = "Comparison pairs 'Case1_vs_Ctrl1~Cast2_vs_Ctrl2' (default: NULL)",
  metavar = "COMPARISON_PAIR"
)

parser$add_argument(
  '--application',
  dest = "c.appl",
  type = "character",
  required = FALSE,
  help = "Application c('RNAseq') (default: NULL)",
  metavar = "APPL"
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

args$b.filter.low.expr <- str2bool(as.character(args$b.filter.low.expr))
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

if(! file.exists(args$c.metadata.path)) {
  flog.error(paste("CANNOT find meta data file with", args$c.metadata.path)); print(sessionInfo()); q(save = "no")
}
df.meta.tbl <- read.csv(
  args$c.metadata.path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

df.meta.tbl <- df.meta.tbl[! duplicated(df.meta.tbl$Library), ]

if(args$c.appl == "RNAseq") {
  v.in.file.path <- list.files(
    path = args$c.in.dir.path,
    pattern = "\\.count$",
    recursive = TRUE,
    full.names = TRUE
  )
  if(length(v.in.file.path) == 0) {
    flog.error(paste("CANNOT find any count tables under", args$c.in.dir.path)); print(sessionInfo()); q(save = "no")
  }
  
  flog.info(paste(length(v.in.file.path), "input file(s) identified"))
}

flog.debug("PREPROCESS - DONE")


flog.info("limma: Linear Models for Microarray Data")

if(! "limma" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("limma", version = "devel")
}
library("limma") # help(package = "limma") # ls(pos = "package:limma")


flog.info("The count table + The metadata")

df.t <- read.csv(
  v.in.file.path[1],
  sep = "\t",
  header = FALSE,
  check.names = FALSE,
  nrows = 2,
  comment.char = "#",
  stringsAsFactors = FALSE
)

df.counts <- NULL
df.gene.annot <- NULL

if(length(v.in.file.path) > 1) {
  for(c.in.file.path in v.in.file.path) {
    df.count <- read.csv(
      c.in.file.path,
      sep = "\t",
      header = ifelse(df.t[1, 1] == "Geneid", TRUE, FALSE),
      check.names = FALSE,
      comment.char = "#",
      stringsAsFactors = FALSE
    )
    
    if(any(apply(df.count[1, -c(1), drop = FALSE], 2, is.character))) {
      if(is.null(df.gene.annot)) {
        df.gene.annot <- df.count[, c(1:(ncol(df.count)-1)), drop = FALSE]
        rownames(df.gene.annot) <- df.gene.annot[, c(1)]
      }
    }
    
    if(is.null(df.counts)) {
      df.counts <- df.count[, c(1, ncol(df.count))]
    } else {
      df.counts <- merge(df.counts, df.count[, c(1, ncol(df.count))], by = colnames(df.counts)[1], sort = FALSE)
    }
  }
  
  rownames(df.counts) <- df.counts[, c(1)]
  df.counts <- df.counts[, -c(1)]
  
  (v.fn.2.library <- sapply(colnames(df.counts), function(x) { basename(dirname(x)) }))
  (v.fn.2.group <- df.meta.tbl[match(v.fn.2.library, df.meta.tbl$Library), "Group"])
  names(v.fn.2.group) <- names(v.fn.2.library)
  colnames(df.counts) <- v.fn.2.library
}


rownames(df.meta.tbl) <- df.meta.tbl$Library
df.meta.tbl <- df.meta.tbl[match(rownames(df.meta.tbl), colnames(df.counts)), ]
df.meta.tbl$Group <- as.factor(df.meta.tbl$Group)
mat.design <- stats::model.matrix( ~ 0 + Group, df.meta.tbl)


if(! "edgeR" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("edgeR", version = "devel")
}
library("edgeR") # help(package = "edgeR") # ls(pos = "package:edgeR")


system.time(
  DGEList.s <- edgeR::DGEList(
    counts = round(df.counts),
    samples = df.meta.tbl,
    genes = c(NULL, df.gene.annot)[1],
    remove.zeros = TRUE
  )
)


c.out.file.t.path <- gsub("\\.count.txt$", ".RData", args$c.out.file.path)
save(
  args,
  DGEList.s,
  file = c.out.file.t.path
)


flog.debug("The count table + The metadata - DONE")


if(args$b.filter.low.expr) {
  flog.info("FILTER for lowly expressed genes")
  
  DGEList.t <- DGEList.s[apply(edgeR::getCounts(DGEList.s), 1, function(x) { all(x >= 10) }), , keep.lib.size = FALSE]
  flog.info(paste("From", nrow(edgeR::getCounts(DGEList.s)), "to", nrow(edgeR::getCounts(DGEList.t))))
  
  flog.debug("FILTER for lowly expressed genes - DONE")
} else {
  DGEList.t <- DGEList.s
}


flog.info("Normalziation & Variance estimation")

system.time(
  DGEList.sizeFact <- edgeR::calcNormFactors(
    DGEList.t,
    method = "TMM"
  )
)

mat.norm <- t(t(DGEList.sizeFact$counts) * DGEList.sizeFact$samples$norm.factors)

flog.debug(paste("OUT_FILE_PATH:", args$c.out.file.path, sep = "\t"))
write.table(
  mat.norm,
  args$c.out.file.path,
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)


c.out.file.t.path <- gsub("\\.count.txt$", ".voom.mean-variance.pdf", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
pdf(c.out.file.t.path, width = 8, height = 8)
system.time(
  EList.voom <- limma::voom(
    counts = DGEList.sizeFact,
    design = mat.design,
    plot = TRUE,
    save.plot = TRUE,
    normalize.method = "quantile",
    span = 0.5
  )
)
dev.off()

mat.obj <- EList.voom


c.out.file.t.path <- gsub("\\.count.txt$", ".count.voom.meanSdPlot.pdf", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

pdf(c.out.file.t.path, width = 8, height = 8)
vsn::meanSdPlot(as.matrix(EList.voom))
dev.off()

flog.debug("Normalziation & Variance estimation - DONE")


flog.info("Inference: Calling differential expression")

system.time(
  MArrayLM.t <- limma::lmFit(
    object = mat.obj,
    design = mat.design,
    method = "ls"
  )
)

v.t <- unlist(strsplit(args$c.v.case.2.ctrl, split = "_vs_"))

system.time(
  mat.contrast.t <- limma::makeContrasts(
    contrasts = paste(paste0("Group", v.t), collapse = "-"),
    levels = mat.design
  )
)

MArrayLM.c <- MArrayLM.t

system.time(
  MArrayLM.s <- limma::contrasts.fit(
    fit = MArrayLM.c,
    contrasts = mat.contrast.t
  )
)

system.time(
  MArrayLM.s <- limma::eBayes(
    fit = MArrayLM.s,
    proportion = 0.01,
    stdev.coef.lim = c(0.1, 4),
    trend = FALSE,
    robust = FALSE,
    winsor.tail.p = c(0.05, 0.1)
  )
)

flog.debug("Inference: Calling differential expression - DONE")


flog.info("Extracting results from a limma analysis")

system.time(
  df.limma.t <- limma::topTable(
    fit = MArrayLM.s,
    number = Inf,
    adjust.method = "BH",
    sort.by = "B",
    p.value = 1,
    lfc = 0,
    confint = FALSE
  )
)


c.out.file.t.path <- gsub("\\.count.txt$", ".plotMA.pdf", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

df.t <- as.data.frame(df.limma.t) # Geneid, logFC, AveExpr, t, P.Value, adj.P.Val, B

colnames(df.t) <- gsub(
  "AveExpr",
  "baseMean",
  colnames(df.t)
)

colnames(df.t) <- gsub(
  "logFC",
  "log2FoldChange",
  colnames(df.t)
)

colnames(df.t) <- gsub(
  "adj.P.Val",
  "FDR",
  colnames(df.t)
)

pdf(c.out.file.t.path, width = 8, height = 8)
if(args$c.criteria == 'q') {
  flog.info(paste("FILTER by q value (FDR) < ", args$d.pval))
  print(table(DEG = df.t$DEG <- as.integer(df.t$FDR < args$d.pval) * sign(df.t$log2FoldChange), sign = sign(df.t$log2FoldChange)))
  
  print(plot(x = df.t$log2FoldChange, y = df.t$FDR, col = ifelse((df.t$FDR < args$d.pval), "red3", "gray32"), pch = 19, cex = .3, xlim = c(-4, 4)))
  DESeq::plotMA(df.t, col = ifelse(df.t$FDR < args$d.pval, "red3", "gray32"))
} else if(args$c.criteria == 'fc_q') {
  flog.info(paste("FILTER by abs(log2(fc)) >= ", log2(args$d.fold.change), "& q value (FDR) < ", args$d.pval))
  print(table(DEG = df.t$DEG <- as.integer((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) & (df.t$FDR < args$d.pval)) * sign(df.t$log2FoldChange), sign = sign(df.t$log2FoldChange)))
  
  print(plot(x = df.t$log2FoldChange, y = df.t$FDR, col = ifelse(((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) & (df.t$FDR < args$d.pval)), "red3", "gray32"), pch = 19, cex = .3, xlim = c(-4, 4)))
  DESeq::plotMA(df.t, col = ifelse(((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) & (df.t$FDR < args$d.pval)), "red3", "gray32"))
} else if(args$c.criteria == 'fc_p') {
  flog.info(paste("FILTER by abs(log2(fc)) >= ", log2(args$d.fold.change), "& p value < ", args$d.pval))
  print(table(DEG = df.t$DEG <- as.integer((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) & (df.t$PValue < args$d.pval)) * sign(df.t$log2FoldChange), sign = sign(df.t$log2FoldChange)))
  
  print(plot(x = df.t$log2FoldChange, y = df.t$PValue, col = ifelse(((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) & (df.t$PValue < args$d.pval)), "red3", "gray32"), pch = 19, cex = .3, xlim = c(-4, 4)))
  DESeq::plotMA(df.t, col = ifelse(((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) & (df.t$PValue < args$d.pval)), "red3", "gray32"))
} else if(args$c.criteria == 'fc') {
  flog.info(paste("FILTER by abs(log2(fc)) >= ", log2(args$d.fold.change)))
  print(table(DEG = df.t$DEG <- as.integer(abs(df.t$log2FoldChange) >= log2(args$d.fold.change)) * sign(df.t$log2FoldChange), sign = sign(df.t$log2FoldChange)))
  
  print(plot(x = df.t$log2FoldChange, y = df.t$PValue, col = ifelse((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)), "red3", "gray32"), pch = 19, cex = .3, xlim = c(-4, 4)))
  DESeq::plotMA(df.t, col = ifelse((abs(df.t$log2FoldChange) >= log2(args$d.fold.change)), "red3", "gray32"))
}
dev.off()


c.out.file.t.path <- gsub("\\.count.txt$", ".test.txt", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
write.table(
  df.t,
  file = c.out.file.t.path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

c.out.file.t.path <- gsub("\\.count.txt$", ".test.DEG.txt", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
write.table(
  subset(df.t, DEG != 0),
  file = c.out.file.t.path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

flog.debug("Extracting results from a limma analysis - DONE")


flog.debug("limma: Empirical Analysis of Digital Gene Expression Data in R - DONE")



print(sessionInfo())
q(save = "no")
