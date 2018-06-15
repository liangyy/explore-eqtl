colorGTEx <- function(col, color.name, n) {
  col <- factor(col, levels = 1 : n, labels = as.character(color.name) )
  return(col)
}

# obtained from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# end

getDistance <- function(pos, bed) {
  e <- apply(pos, 1, function(x) {
    chr <- x[1]
    bed.sub <- bed[bed[,1] == chr, ]
    if(nrow(bed) == 0) {
      return(NA)
    }
    return(min(abs(bed.sub[, c(2, 3)] - as.numeric(x[2]))))
  })
  return(e)
}

getDistanceToTSS <- function(x) {
  # x: chr, pos, gene.chr, gene.start, gene.end, strand
  tts <- x[, 4]
  neg.ind <- x[, 6] == '-' | x[, 6] == -1
  tts[neg.ind] <- x[neg.ind, 5]
  d <- x[, 2] - tts
  d[neg.ind] <- -d[neg.ind]
  d[x[, 1] != x[, 3]] <- NA
  return(d)
}
