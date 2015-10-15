# unix305 /usr/local/work/hickey/tmp/
# 
library(data.table)
library(Matrix)
library(S4Vectors)
library(repete)
library(microbenchmark)

x <- fread("test.4.tsv")
# TODO: Does as.matrix() incur a copy?
y <- as.matrix(x[, -c(1:6), with = FALSE])
rm(x)

y_dense <- Matrix(y, sparse = FALSE)
y_sparse <- Matrix(y, sparse = TRUE)
y_DF <- DataFrame(lapply(seq_len(ncol(y)), 
                                function(i, y) {
                                  Rle(y[, i])
                                }, y = y))
colnames(y_DF) <- colnames(y)

lsos()

# Name      Type       Size         PrettySize    Rows          Columns
# y_dense   dgeMatrix  5658884864     5.3 Gb      44210022      16
# y         matrix     2829442712     2.6 Gb      44210022      16
# y_sparse  dgCMatrix   794213424   757.4 Mb      44210022      16
# y_DF      DataFrame   718667376   685.4 Mb      44210022      16

# `Matrix::Matrix` objects are designed to be worked on as whole objects, not 
# element-wise (https://stat.ethz.ch/pipermail/r-help/2012-June/316535.html). 
# This, at least partially, explains why accessing individual elements is so 
# slow. The same thread suggests using Eigen matrices instead.

microbenchmark(y[666666, 4],
               y_dense[666666, 4],
               y_sparse[666666, 4],
               y_DF[666666, 4], times = 10)

# colMeans(y_DF) doesn't work out-of-the-box
microbenchmark(colMeans(y),
               colMeans(y_dense),
               colMeans(y_sparse), times = 10)