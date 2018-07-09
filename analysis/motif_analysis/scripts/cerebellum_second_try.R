# this function works on factor eleven (cerebellum factor)
# the idea is to extract cerebellum-specific eQTLs

getPosNegEQTL <- function(flash, mash) {
  l.five <- flash$EL[, 11]
  l.one <- flash$EL[, 1]
  idx <- 1 : length(l.five)
  cerebellum.idx <- c(10, 9)
  # cond1 <- abs(l.five) > 8
  mash.mean <- mash$posterior$result$PosteriorMean
  mash.std <- mash$posterior$result$PosteriorSD
  mash.ratio <- abs(mash.mean) / mash.std
  mash.sig <- mash.ratio > 2
  mash.super <- mash.ratio > 4
  o <- rowSums(mash.sig[, !(1 : 44 %in% cerebellum.idx)])
  h <- rowSums(mash.super[, (1 : 44 %in% cerebellum.idx)])
  cerebellum.idx <- idx[o < 3 & h > 1 | o < 5 & h > 1 & abs(l.five) > 5]
  
  control.idx <- idx[(abs(l.five) < 1 & abs(l.one) > 50)]
  return(list(pos = cerebellum.idx, neg = control.idx))
}
