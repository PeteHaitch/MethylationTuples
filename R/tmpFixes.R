### -------------------------------------------------------------------------
### Temporary fixes: Patches submitted that are yet to be implemented into pkg
###

# Contributed (as an S3 method) to R.utils 
# (https://github.com/HenrikBengtsson/R.utils/pull/1)
#' @export
isBzipped <- function(filename, method = c("extension", "content"), ...) {
  # Argument 'method':
  method <- match.arg(method)
  
  # Argument 'filename':
  filename <- Arguments$getReadablePathname(filename, 
                                            mustExist=(method == "content"))
  
  if (method == "extension") {
    res <- (regexpr("[.]bz2$", filename, ignore.case=TRUE) != -1L)
  } else if (method == "content") {
    con <- file(filename)
    on.exit(close(con))
    res <- (summary(con)$class == "bzfile")
  }
  
  res
}
