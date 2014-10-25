#methpat <- mp1
methpat_order <- order(methpat)
methpat_rd_sorted <- rowData(methpat)[methpat_order]
betas <- betaVal(methpat)
#ipd <- 1:2000
ipd <- sort(unique(diff(start(methpat_rd_sorted))))
ipd <- ipd[ipd > 0]
id_dt <- setDT(expand.grid(IPD = ipd, 
                           strand = levels(strand(methpat)), 
                           feature_status = 0:3))
id_dt[, c("KEY", "ID") := list(paste(IPD, strand, feature_status, sep = ''),
                               seq_len(nrow(id_dt)))]
setkey(id_dt, ID)
method <- 'pearson'

seqnames <- as.character(seqnames(methpat_rd_sorted))
strand <- as.character(strand(methpat_rd_sorted))
pos <- start(methpat_rd_sorted)
feature_status <- rep(FALSE, length(methpat_rd_sorted))

system.time({pairs <- setkey(setDT(
  .makeAllPairsCpp(
    methpat_order, seqnames, strand, pos, feature_status, ipd, betas, id_dt)), 
  ID, sample)})

system.time({cors <- pairs[, list(cor = cor(beta1, beta2, use = "na.or.complete", 
                                            method = method)), by = list(ID, sample)]})

system.time({pairs <- setkey(setDT(
  .makeAdjacentPairsCpp(
    methpat_order, seqnames, strand, pos, feature_status, betas, id_dt)), 
  ID, sample)})

system.time({cors <- pairs[, list(cor = cor(beta1, beta2, use = "na.or.complete", 
                                            method = method)), by = list(ID, sample)]})


pairs <- setkey(setDT(.makeAdjacentPairsCpp(
  methpat_order,
  as.character(seqnames(methpat_rd_sorted)), 
  as.character(strand(methpat_rd_sorted)),
  start(methpat_rd_sorted),   
  feature_status, 
  betas, id_dt)), ID, sample)

pairs <- setkey(setDT(makeAdjacentPairs(
  methpat_order,
  as.character(seqnames(methpat_rd_sorted)), 
  as.character(strand(methpat_rd_sorted)),
  start(methpat_rd_sorted),   
  feature_status, 
  betas, id_dt)), ID, sample)


a <- betaCor(methpat[1:1000], min_cov = 5L, pair_type = 'adjacent', feature = cgi, feature_name = 'CGI')