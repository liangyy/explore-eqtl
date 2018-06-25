# this function works on factor five (heart factor)
# the idea is to extract heart-specific eQTLs

getPosNegEQTL <- function(flash, mash) {
  l.five <- flash$EL[, 5]
  l.one <- flash$EL[, 1]
  heart.idx <- c(25, 26, 29)
  cond1 <- abs(l.five) > 8
  mash.mean <- mash$posterior$result$PosteriorMean
  o.mean <- rowMeans(mash.mean[, !(1 : 44 %in% heart.idx)])
  h.mean <- Biobase::rowMax(mash.mean[, (1 : 44 %in% heart.idx)])
  cond2 <- abs(o.mean) < abs(h.mean) * 0.2
  mash.heart.ind <- cond2 & cond1
  e <- seq(1 : length(mash.heart.ind))[mash.heart.ind]
  exclude.idx <- c(384, 828, 3218, 4811, 5396, 6932, 7361, 9557, 10602, 11901)
  heart.idx <- e[!e %in% exclude.idx]
  idx <- 1 : length(l.five)
  control.idx <- idx[(abs(l.five) < 1 & abs(l.one) > 50)]
  return(list(pos = heart.idx, neg = control.idx))
}
