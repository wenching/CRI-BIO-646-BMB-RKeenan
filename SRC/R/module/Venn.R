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
  dest = "c.v.in.file.path",
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
  '--anchorPath',
  dest = "c.anchor.file.path",
  type = "character",
  required = FALSE,
  help = "Anchor file (default: NULL)",
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

if(! "Biobase" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
    
  BiocManager::install("Biobase", version = "devel")
}
library("Biobase")


v.in.file.path <- unlist(strsplit(x = args$c.v.in.file.path, split = ","))

flog.info(paste("It seems like to compare across different methods."))
flog.info(paste("COMPARISON_AMONG_DIFFERENT_METHODS\t", b.comp.among.method <- TRUE))
names(v.in.file.path) <- gsub(pattern = Biobase::lcPrefix(basename(v.in.file.path)), replacement = "", x = gsub(pattern = Biobase::lcSuffix(basename(v.in.file.path)), replacement = "", x = basename(v.in.file.path)))
flog.info(paste("#VAR:\t", "OUT_BASE\t", c.out.base <- gsub(pattern = "\\.$", replacement = "", Biobase::lcPrefix(basename(v.in.file.path)))))

lst.label.2.df <- list()
lst.label.2.df.id <- list()

for(c.label in names(v.in.file.path)) {
  (c.in.file.path <- v.in.file.path[c.label])
  
  flog.debug(paste("LABEL\t[", c.label, "]"))
  if(! file.exists(c.in.file.path)) {
    flog.warn(paste("cannot find", c.in.file.path))
    next
  }
  df.t <- read.csv(
    c.in.file.path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  if(nrow(df.t) == 0) {
    flog.warn(paste("SKIP", c.label, "due to no data"))
    next
  }
  
  v.idx.anc <- 1
  
  lst.label.2.df[[c.label]] <- df.t
  
  lst.label.2.df.id[[c.label]] <- df.t[, v.idx.anc, drop = FALSE]
  
  flog.info(paste("LABEL\t[", c.label, "] (", nrow(lst.label.2.df.id[[c.label]]), "entries ) - DONE"))
}


v.label.2.anchorCol <- structure(
  rep("Geneid", length.out = length(lst.label.2.df.id)),
  .Names = names(lst.label.2.df.id)
)

lst.label.2.id <- list()

if(b.comp.among.method) {
  df.id.2.DEG.rbind <- data.frame()
}

for(c.label in names(v.label.2.anchorCol)) {
  flog.debug(paste("LABEL:\t[", c.label, "]"))
  
  head(df.t <- lst.label.2.df[[c.label]])
  (i.idx.anc <- which(v.label.2.anchorCol[c.label] %in% names(df.t)))
  
  if(any(duplicated(df.t[, i.idx.anc]))) flog.info(paste("REMOVE duplicated ids in", colnames(df.t)[i.idx.anc], "from", nrow(df.t), "to", nrow(df.t <- df.t[! duplicated(df.t[, i.idx.anc]), ]), c.label))
  
  head(lst.label.2.df[[c.label]] <- df.t)
  
  lst.label.2.id[[c.label]] <- sort(as.vector(lst.label.2.df[[c.label]][, i.idx.anc, drop = T]))
  
  colnames(df.t)[i.idx.anc] <- "..anchor.."
  
  if(b.comp.among.method) {
    df.id.2.DEG.rbind <- rbind(df.id.2.DEG.rbind, df.t[, c("..anchor..", "DEG")])
  }
}

flog.info(paste("#STAT:\t", "RANGE STATISTICS"))
sapply(lst.label.2.id, length)


flog.info(paste("#STAT:\t", "Up/Down regulated DEGs per methods"))
addmargins(sapply(lst.label.2.df, function(x) { table(factor(as.character(x$DEG), levels = c("-1", "0", "1"))) }), 1)


if(b.comp.among.method) {
  head(lst.id.2.DEG <- sapply(split(df.id.2.DEG.rbind$DEG, df.id.2.DEG.rbind$..anchor..), function(x) { names(sort(table(x), decreasing = T))[1] }))
  head(df.id.2.DEG <- data.frame(..anchor.. = names(lst.id.2.DEG), DEG = as.numeric(lst.id.2.DEG), stringsAsFactors = F))
  
  # Sanity Checking
  for(i.idx in 1:min(c(5, nrow(df.id.2.DEG)))) {
    flog.debug(paste("The consensus trend of", df.id.2.DEG$..anchor..[i.idx], "gene is", c("1" = "up-regulated", "-1" = "down-regulated", "0" = "no change")[as.character(df.id.2.DEG$DEG[i.idx])]))
    flog.debug(paste("The individual trends of", df.id.2.DEG$..anchor..[i.idx], "gene are",  paste(c("1" = "up-regulated", "-1" = "down-regulated", "0" = "no change")[as.character(df.id.2.DEG.rbind[which(df.id.2.DEG.rbind$..anchor.. == df.id.2.DEG$..anchor..[i.idx]), "DEG"])], collapse = " | "), "\n"))
  }
}


flog.debug("PREPROCESS - DONE")


if(! "Vennerable" %in% rownames(installed.packages())) {
  if(! "RBGL" %in% rownames(installed.packages())) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("RBGL")
  }
  
  if(! "graph" %in% rownames(installed.packages())) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("graph")
  }
  
  if(! "reshape" %in% rownames(installed.packages())) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("reshape")
  }
  
  if(! "devtools" %in% rownames(installed.packages())) {
    install.packages(
      "devtools",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  library("devtools")
  
  install_github(repo = 'js229/Vennerable')
  #browseURL("http://github.com/js229/Vennerable")
}
library("Vennerable")
# ls(pos = "package:Vennerable")

s4.venn <- Vennerable::Venn(Sets = lst.label.2.id)
sapply(lst.venn <- s4.venn@IntersectionSets, length)


source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")

system.time(
  OLlist <- overLapper(
    setlist = lst.label.2.id,
    sep = "&",
    type = "vennsets"
  )
)

df.venn.list <- data.frame(
  Methods = names(OLlist$Venn_List),
  Method.Num = unlist(
    lapply(
      names(OLlist$Venn_List),
      function(x) { length(unlist(strsplit(x, split = "&"))) }
    )
  ),
  ID = unlist(
    lapply(
      OLlist$Venn_List,
      function(x) { paste(x, collapse = "~")}
    ),
    use.names =  FALSE
  ),
  ID.Num = unlist(
    lapply(
      OLlist$Venn_List,
      length
    ),
    use.names =  FALSE
  ),
  stringsAsFactors = FALSE
)


for(c.set in c("CONSENSUS", "OVERLAP", "UNION")) {
  if(c.set == "CONSENSUS") {
    v.t <- which(df.venn.list$Method.Num == max(df.venn.list$Method.Num))
  } else if(c.set == "OVERLAP") {
    v.t <- which(df.venn.list$Method.Num > 1)
  } else {
    v.t <- which(df.venn.list$Method.Num > 0)
  }
  
  flog.info(paste(c.set, "IDs\t", sum(df.venn.list[v.t, "ID.Num"])))
}


c.out.file.t.path <- gsub("\\.overlap.txt$", ".VennList.txt", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

write.table(
  df.venn.list,
  file = c.out.file.t.path,
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)



if(! "VennDiagram" %in% rownames(installed.packages())) {
  install.packages(
    "VennDiagram",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}
library("VennDiagram")

venn.plot <- VennDiagram::venn.diagram(
  lst.label.2.id,
  filename = NULL,
  col = "transparent", 
  fill = RColorBrewer::brewer.pal(n = max(3, length(lst.label.2.id)), name = "Set1")[1:length(lst.label.2.id)],
  category.names = paste(names(lst.label.2.id), unlist(sapply(lst.label.2.id, length)), sep = ": "),
  cex = 3,
  margin = 0.1
)


c.out.file.t.path <- gsub("\\.overlap.txt$", ".VennDiagram.pdf", args$c.out.file.path)
flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))

pdf(file = c.out.file.t.path, width = 8, height = 8)
grid::grid.draw(venn.plot)
dev.off()


if(!is.null(args$c.anchor.file.path) && file.exists(args$c.anchor.file.path)) {
  df.anchor <- read.csv(
    file = args$c.anchor.file.path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  i.idx.anc <- 1
  
  if(b.comp.among.method) {
    flog.info(paste("UPDATE consensus trends of DEGs to anchor files"))
    
    if(length(setdiff(df.id.2.DEG[,1], df.anchor[, i.idx.anc, drop = TRUE])) != 0) {
      flog.info(paste(paste(setdiff(df.id.2.DEG[,1], df.anchor[, i.idx.anc, drop = T]), collapse = " | "), "anchros couldn't be found in the anchor file")); print(sessionInfo()); q(save = "no")
    }
    
    flog.debug(paste("RESET the DEG column in the anchor file to NA"))
    df.anchor$DEG <- NA
    flog.debug(paste("REASSIGN consensus trends of DEGs among different methods"))
    head(df.anchor[match(df.id.2.DEG[,1], df.anchor[, i.idx.anc]), "DEG"] <- df.id.2.DEG$DEG)
    
    table(df.anchor$DEG, useNA = "ifany")
  }
  
  
  if((nchar(names(lst.venn)[length(lst.venn)]) > 1) && ! (is.na(lst.venn[[length(lst.venn)]]) || is.null(lst.venn[[length(lst.venn)]]))) {
    flog.info(paste(length(v.id.cons <- sort(unlist(lst.venn[[length(lst.venn)]]))), "consensus id(s).")) # consensus
    
    flog.info(paste("#STAT[CONSENSUS]:\t", "Up/Down regulated DEGs"))
    if(b.comp.among.method) print(table(factor(c("-1" = "DOWN", "1" = "UP")[as.character(df.id.2.DEG[which(df.id.2.DEG[,1] %in% v.id.cons), "DEG"])], levels = c("DOWN", "UP"))))
    
    flog.info(paste("SELECT", nrow(df.cons <- df.anchor[which(df.anchor[, i.idx.anc] %in% v.id.cons), ]), "out of", nrow(df.anchor), "as consensus id(s)."))
  }
  
  c.out.file.t.path <- gsub("\\.overlap.txt$", ".consensus.txt", args$c.out.file.path)
  flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  options(scipen = 999)
  write.table(
    df.cons,
    file = c.out.file.t.path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  options(scipen = 0)
  
  
  flog.info(paste(length(v.id.union <- sort(unlist(lst.venn, use.names = F))), "union id(s).")) # all/union
  
  flog.info(paste("#STAT[UNION]:\t", "Up/Down regulated DEGs"))
  if(b.comp.among.method) print(table(factor(c("-1" = "DOWN", "1" = "UP")[as.character(df.id.2.DEG[which(df.id.2.DEG[,1] %in% v.id.union), "DEG"])], levels = c("DOWN", "UP"))))
  
  flog.info(paste("SELECT", nrow(df.union <- df.anchor[which(df.anchor[, i.idx.anc] %in% v.id.union), ]), "out of", nrow(df.anchor), "as union id(s)."))
  
  c.out.file.t.path <- gsub("\\.overlap.txt$", ".union.txt", args$c.out.file.path)
  flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  options(scipen = 999)
  write.table(
    df.union,
    file = c.out.file.t.path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  options(scipen = 0)
  
  
  flog.info(paste(length(v.id.overlap <- sort(unlist(lst.venn[sapply(names(lst.venn), function(x) { sum(as.numeric(unlist(strsplit(x, split = "")))) > 1 })], use.names = F))), "overlapped id(s).")) # overlapped
  
  flog.info(paste("#STAT[OVERLAP]:\t", "Up/Down regulated DEGs"))
  if(b.comp.among.method) print(table(factor(c("-1" = "DOWN", "1" = "UP")[as.character(df.id.2.DEG[which(df.id.2.DEG[,1] %in% v.id.overlap), "DEG"])], levels = c("DOWN", "UP"))))
  
  flog.info(paste("SELECT", nrow(df.overlap <- df.anchor[which(df.anchor[, i.idx.anc] %in% v.id.overlap), ]), "out of", nrow(df.anchor), "as overlapped id(s)."))
  
  
  c.out.file.t.path <- args$c.out.file.path
  flog.debug(paste("OUT_FILE_PATH:", c.out.file.t.path, sep = "\t"))
  
  options(scipen = 999)
  write.table(
    df.overlap,
    file = c.out.file.t.path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  options(scipen = 0)
  
  flog.info("")
}

print(sessionInfo())
q(save = "no")
