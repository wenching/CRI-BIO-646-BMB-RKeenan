#' pheatmap.R
#' 
#' @description pheatmap
#' 
#' run pheatmap
#' 
#' @param GTF character. GTF file
#' @param inPath character vector. input
#' @param outPath character vector. output
#' @usage postAna::pheatmap(inPath, outPath, ...)
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
#' postAna::pheatmap("overlap.txt", "heatmap.pdf", "count.txt")



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
  "pheatmap.R",
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
  dest = "c.in.file.path",
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
  '--exprTblPath',
  dest = "c.expr.tbl.path",
  type = "character",
  required = TRUE,
  help = "expression table (default: required arg)",
  metavar = "EXPR_TBL"
)

parser$add_argument(
  '--assembly',
  dest = "c.assembly",
  type = "character",
  required = FALSE,
  help = "assembly (default: NULL)",
  metavar = "ASSEMBLY"
)

parser$add_argument(
  '--logPath',
  dest = "c.log.file.path",
  type = "character",
  required = FALSE,
  help = "log file (default: NULL)",
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


df.in <- read.csv(
  file = args$c.in.file.path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

df.expr <- read.csv(
  file = args$c.expr.tbl.path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

df.expr.vst <- as.data.frame.matrix(
  DESeq2::vst(
    round(as.matrix.data.frame(df.expr)),
    blind = TRUE
  )
)

df.expr.vst$Geneid <- row.names(df.expr.vst)

df.expr.vst.subset <- df.expr.vst[which(df.expr.vst$Geneid %in% df.in$Geneid), -ncol(df.expr.vst)]

df.expr.vst.subset.scaled <- as.data.frame.matrix(
  t(
    apply(
      df.expr.vst.subset,
      1,
      scale
    )
  )
)

colnames(df.expr.vst.subset.scaled) <- colnames(df.expr.vst.subset)

if(! "pheatmap" %in% rownames(installed.packages())) {
  install.packages(
    "pheatmap",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}
library("pheatmap")

h <- pheatmap::pheatmap(
  mat = df.expr.vst.subset.scaled,
  show_rownames = FALSE,
  filename = NA
)

c.out.file.t.path <- args$c.out.file.path
flog.info(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

pdf(c.out.file.t.path, width = 6, height = 8)
print(
  h
)
dev.off()

flog.info("PREPROCESS - DONE")


print(sessionInfo())
q(save = "no")
