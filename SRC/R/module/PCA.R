#' PCA.R
#' 
#' @description PCA
#' 
#' run PCA
#' 
#' @param metadata character. meta data file
#' @param GTF character. GTF file
#' @param inPath character vector. input
#' @param outPath character vector. output
#' @usage callLoci::PCA(inPath, outPath, ...)
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
#' callLoci::PCA("expr.txt", "pca.pdf")



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
  "PCA.R",
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
  '--ids',
  dest = "c.v.id",
  type = "character",
  required = FALSE,
  help = "keys (default: NULL)",
  metavar = "LIB1,LIB2,..."
)

parser$add_argument(
  '--groups',
  dest = "c.v.group",
  type = "character",
  required = FALSE,
  help = "values (default: NULL)",
  metavar = "GRP1,GRP2,..."
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
if(args$c.v.id != "NULL") {
  args$v.library.2.group <- structure(
    unlist(
      strsplit(
        args$c.v.group,
        split = ","
      )
    ),
    .Names = unlist(
      strsplit(
        args$c.v.id,
        split = ","
      )
    )
  )
}

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

if(! file.exists(args$c.in.file.path)) {
  flog.error(paste("CANNOT find expression table with", args$c.in.file.path)); print(sessionInfo()); q(save = "no")
}
df.expr <- read.csv(
  args$c.in.file.path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

if(! is.null(args$v.library.2.group)) {
  if(length(v.t <- setdiff(colnames(df.expr), names(args$v.library.2.group))) != 0) {
    flog.error("cannot find the group information of ", paste(v.t, collapse = ", "), "in expression table"); print(sessionInfo()); q(save = "no")
  }
}

flog.debug("PREPROCESS - DONE")


flog.info("MAKE PCA plot")

if(all(df.expr >= 0)) df.epxr.log <- log(df.expr + 1)

if(! "FactoMineR" %in% rownames(installed.packages())) {
  install.packages(
    "FactoMineR",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}
library("FactoMineR")

system.time(
  lst.pca <- FactoMineR::PCA(
    t(df.epxr.log), 
    scale.unit = FALSE,
    graph = FALSE
  )
)



if(! "factoextra" %in% rownames(installed.packages())) {
  install.packages(
    "factoextra",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}
library("factoextra")


c.out.file.t.path <- gsub("\\.pca.pdf$", ".eigenvalue_vs_variance.pdf", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

pdf(c.out.file.t.path, width = 8, height = 8)
factoextra::fviz_screeplot(
  lst.pca,
  addlabels = TRUE,
  ylim = c(0, 50)
)
dev.off()


c.out.file.t.path <- gsub("\\.pca.pdf$", ".contribution_of_PC1_n_PC2.pdf", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

pdf(c.out.file.t.path, width = 8, height = 8)
factoextra::fviz_contrib(
  lst.pca,
  choice = "var",
  axes = 1,
  top = 10
)
factoextra::fviz_contrib(
  lst.pca,
  choice = "var",
  axes = 2,
  top = 10
)
dev.off()


if(! is.null(args$v.library.2.group)) {
  v.group <- unique(args$v.library.2.group)
  v.group.2.colour <- RColorBrewer::brewer.pal(n = max(3, length(v.group)), name = "Set1")[1:length(v.group)]
  names(v.group.2.colour) <- v.group
  
  f.habillage <- factor(args$v.library.2.group[as.character(colnames(df.epxr.log))], levels = v.group)
  
  df.2d <- as.data.frame.matrix(lst.pca$ind$coord[, 1:2])
  colnames(df.2d) <- c("x", "y")
  df.2d$group <- f.habillage
  df.2d$label <- rownames(df.2d)
}


c.out.file.t.path <- args$c.out.file.path
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

pdf(c.out.file.t.path, width = 8, height = 8)
if(! is.null(args$v.library.2.group)) {
  factoextra::fviz_pca_ind(
    lst.pca,
    habillage = f.habillage,
    palette = v.group.2.colour
  ) + 
    stat_ellipse(mapping = aes(colour = group), data = rbind(df.2d, df.2d), geom = "polygon", alpha = 1/8) + 
    theme_classic()
} else {
  factoextra::fviz_pca_ind(
    lst.pca
  )
}
dev.off()


flog.debug("MAKE PCA plot - DONE")


print(sessionInfo())
q(save = "no")
