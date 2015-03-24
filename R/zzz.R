.onUnload <- function(libpath) {
  library.dynam.unload("MethylationTuples", libpath)
}