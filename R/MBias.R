### =========================================================================
### MBias: An S4 class to store M-bias data.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### data.frame('read', 'readPos', 'methylationType', 'M', 'U', 'total',
###            'percM')
###             read: A factor encoding whether the read is read-1 (R1) or 
###                   read-2 (R2). Single-end data are all encoded as read-1.
###         readPos: An integer giving the read-position.
### methylationType: A factor encoding the methylation type as CG, CHG or CHH.
###                M: An integer encoding the number of times that 
###                   read-position was methylated.
###                U: An integer encoding the number of times that 
###                   read-position was unmethylated.
###            total: An integer encoding the number of times that 
###                   read-position was observed.
###           percM: The percentage of times that read-positions was 
###                   methylated (M / (M + U)).
### Basically just a data.frame. Probably a little heavy-duty, but in keeping
### with the rest of the package which is S4-based.

#' MBias class
#'
#' An S4 class to store M-bias data.
#'
#' @aliases MBias
#'
#' @export
setClass("MBias", 
         contains = "data.frame"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.MBias.columns <- function(object) {
  msg <- NULL
  if (!is.data.frame(object)) {
    msg <- Biobase::validMsg(msg, "'Mbias' should extend 'data.frame'.")
  } else {
    if (!identical(colnames(object), c("sampleName", "read", "readPos", 
                                       "methylationType",  "M", "U", "total", 
                                       "percM"))) {
      msg <- Biobase::validMsg(msg, paste0("Column names should be ",
                                           "'sampleName', 'read', ",
                                           "'readPos', 'methylationType', ",
                                           "'M', 'U', 'total' and 'percM'"))
    }
    else {
      if (!is.character(object[["sampleName"]])) {
        msg <- Biobase::validMsg(msg, paste0("'sampleName' shouuld be a ",
                                             "character."))
      }
      if (!is.factor(object[["read"]]) || 
          !identical(levels(object[["read"]]), c("R1", "R2"))) {
        msg <- Biobase::validMsg(msg, paste0("'read' column should be a ", 
                                             "factor with levels 'R1' and ", 
                                             "'R2'."))
      }
      if (!is.integer(object[["readPos"]])) {
        msg <- Biobase::validMsg(msg, paste0("'readPos' column should be ", 
                                             "an integer."))
      }
      if (!is.factor(object[["methylationType"]]) ||
          !identical(levels(object[["methylationType"]]), 
                     c("CG", "CHG", "CHH"))) {
        msg <- Biobase::validMsg(msg, paste0("'methylationType' column should ", 
                                             "be a factor with levels 'CG', ", 
                                             "'CHG' and 'CHH'."))
      }
      if (!is.integer(object[["M"]])) {
        msg <- Biobase::validMsg(msg, "'M' column should be an integer.")
      }
      if (!is.integer(object[["U"]])) {
        msg <- Biobase::validMsg(msg, "'U' column should be an integer.")
      }
      if (!is.integer(object[["total"]])) {
        msg <- Biobase::validMsg(msg, "'total' column should be an integer.")
      }
      if (!is.numeric(object[["percM"]]) || 
          (min(object[["percM"]]) < 0) ||
          (max(object[["percM"]]) > 100)) {
        msg <- Biobase::validMsg(msg, paste0("'percM' column should be ", 
                                             "numeric and between 0 and 100."))
      }
    }
  }
  msg
}

.valid.MBias <- function(object) {
  # Include all .valid.MBias.* functions in this vector
  msg <- c(.valid.MBias.columns(object))
  
  if (is.null(msg)) {
    return(TRUE)
  } else{
    return(msg)
  }
}

S4Vectors::setValidity2("MBias", .valid.MBias)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# None because I don't want the user constructing these manually, rather they
# should be constructed by readMBias().

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readMBIas: Read .M-bias.txt file produced by Bismark
###

# TODO: This is an inefficient function. However, since it only reads small 
# files, performance isn't an issue. Nonetheless, it'd be nice to re-factor 
# to simplify code maintainability.
# 
#' Read Bismark's M-bias file.
#' 
#' @param file A Bismark M-bias file (with \texttt{.M-bias.txt} extension).
#' @param sampleName The name of the sample. Will guess based on \code{file} 
#' if argument is missing.
#' 
#' @return A \code{\link{MBias}} object.
#' @export
readMBias <- function(file, sampleName) {
  
  if (missing(sampleName)) {
    sampleName <- strsplit(basename(file), ".M-bias.txt")[[1]]
  }
  
  # Echo some information
  message(paste0("Reading ", file, " ..."))
  
  # Read the file.
  x <- read.table(file, sep = '\n', as.is = TRUE)
  
  # Find 'header' lines, which are dispersed throughout the file.
  header_lines <- sapply(X = x, function(xx) {
    grepl(pattern = 'context', x = xx)
  })
  y <- vector(mode = "list", length = sum(header_lines))
  names(y) <- x[header_lines]
  
  # z indexes the lines containing M-bias values.
  z <- sapply(X = x, function(xx) {
    grep(pattern = "context|^=|position", x = xx, invert = TRUE)
  })
  dz <- diff(z)
  next_context_idx <- c(which(dz > 1), nrow(z))
  colnames_idx <- sapply(X = x, function(xx) {
    grep(pattern = "position", x = xx)
  })[1]
  colnames <- strsplit(x = x[colnames_idx, ], split = "\t")[[1]]
  
  # Loop over the lines in the file.
  for (i in seq_along(next_context_idx)) {
    # Note the "+4" to skip the intermediate 'header' rows
    start_idx <- ifelse(i == 1, 4, z[next_context_idx[i - 1]] + 4)
    end_idx <- z[next_context_idx[i]]  
    yy <- x[seq(from = start_idx, to = end_idx, by = 1), ]
    y[[i]] <- do.call("rbind", lapply(yy, function(yyy, colnames) {
      tmp <- strsplit(yyy, split = "\t")
      val <- data.frame(as.integer(tmp[[1]][1]),
                        as.integer(tmp[[1]][2]),
                        as.integer(tmp[[1]][3]),
                        as.numeric(tmp[[1]][4]),
                        as.integer(tmp[[1]][5]))
      colnames(val) <- colnames
      return(val)
    }, colnames = colnames))
  }
  y <- do.call("rbind", y)
  y$Context <- strtrim(rownames(y), width = 3)
  
  # Add R1/R2 annotation if paired-end data
  if (grepl(pattern = "R1", x = rownames(y)[1])) {
    y$Read <- ifelse(grepl(pattern = "R1", x = rownames(y)), "R1", "R2")
  } else{
    y$Read <- "R1"
  }
  rownames(y) <- NULL
  
  # Simplify column names
  colnames(y) <- c('readPos', 'M', 'U', 'percM', 'total',
                   'methylationType', 'read')
  
  # Add sampleName
  y$sampleName <- rep(sampleName, nrow(y))
  
  # Re-order columns
  y <- y[, c('sampleName', 'read', 'readPos', 'methylationType', 'M', 'U', 'total',
             'percM')]
  
  # Replace CpG with CG
  y$methylationType <- ifelse(y$methylationType == 'CpG', 'CG',
                               y$methylationType)
  
  # Set column types
  y$read <- factor(y$read, levels = c("R1", "R2"))
  y$methylationType <- factor(y$methylationType, levels = c("CG", "CHG", "CHH"))
  
  # Convert to an MBias object.
  new("MBias", y)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### plot: Plot an MBias object
###

setMethod("plot", 
          signature(x = "MBias", y = "missing"),
          function(x, ...) {
            
            # TODO: Plot raw or normalised M-bias
            
            ggplot(aes(x = readPos, y = percM, colour = methylationType), 
                   data = x) +
              geom_line(lwd = 1.3) + 
              facet_grid(sampleName ~ read) +
              ggtitle("M-bias") +
              ylim(c(0, 100)) + 
              ylab("Percentage methylation") + 
              xlab("Read-position") +
              scale_color_discrete(guide = 
                                     guide_legend(title = "Methylation type"))
            }
)

