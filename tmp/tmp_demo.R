#methpat <- mp1
methpat_order <- order(methpat)
methpat_rd_sorted <- rowData(methpat)[methpat_order]
meth_level <- methLevel(methpat, min_cov = 5L)
# ipd <- 1:2000
ipd <- sort(unique(diff(start(methpat_rd_sorted))))
ipd <- ipd[ipd > 0]
in_feature <- rep(NA, nrow(methpat))
in_feature_levels <- unique(in_feature)
pair_feature_status <- sort(unique(rowSums(expand.grid(in_feature_levels, 
                                                       in_feature_levels))),
                            na.last = FALSE)
id_dt <- setDT(expand.grid(IPD = ipd, 
                           strand = levels(strand(methpat)),
                           pair_feature_status = pair_feature_status))
id_dt[, c("KEY", "ID") := list(paste(IPD, strand, pair_feature_status, 
                                     sep = ''),
                               seq_len(nrow(id_dt)))]
setkey(id_dt, ID)

seqnames <- as.character(seqnames(methpat_rd_sorted))
strand <- as.character(strand(methpat_rd_sorted))
pos <- start(methpat_rd_sorted)
method <- 'pearson'

system.time({pairs_idx <- .makeAllPairsCpp(
  methpat_order, seqnames, strand, pos, in_feature, ipd, id_dt)})

system.time({pairs_idx <- .makeAdjacentPairsCpp(
    methpat_order, seqnames, strand, pos, in_feature, id_dt)})

cors_list <- lapply(colnames(methpat), function(sample_name, 
                                                pairs_idx, betas) {
  beta_pairs <- data.table(ID = pairs_idx[["ID"]], 
                           sample = sample_name,
                           beta1 = betas[pairs_idx[["i"]], sample_name], 
                           beta2 = betas[pairs_idx[["j"]], sample_name])
  beta_pairs[, .my_cor(beta1, beta2, method = method, 
                       conf.level = conf.level), by = list(ID, sample)]
}, pairs_idx = pairs_idx, betas = betas)
cors <- setkey(rbindlist(cors_list), ID, sample)

val <- id_dt[cors]
val[, c("ID", "KEY") := list(NULL, NULL)]