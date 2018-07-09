# this function works on factor three (skin factor)
# the idea is to extract skin-specific eQTLs

getPosNegEQTL <- function(flash, mash) {
  l.five <- flash$EL[, 3]
  l.one <- flash$EL[, 1]
  idx <- 1 : length(l.five)
  skin.idx <- c(36, 35)
  # cond1 <- abs(l.five) > 8
  mash.mean <- mash$posterior$result$PosteriorMean
  mash.std <- mash$posterior$result$PosteriorSD
  mash.ratio <- abs(mash.mean) / mash.std
  mash.sig <- mash.ratio > 2
  mash.super <- mash.ratio > 4
  o <- rowSums(mash.sig[, !(1 : 44 %in% skin.idx)])
  h <- rowSums(mash.super[, (1 : 44 %in% skin.idx)])
  skin.idx <- idx[o < 3 & h > 1 | o < 5 & h > 1 & abs(l.five) > 5]
  
  control.idx <- idx[(abs(l.five) < 1 & abs(l.one) > 50)]
  return(list(pos = skin.idx, neg = control.idx))
}
