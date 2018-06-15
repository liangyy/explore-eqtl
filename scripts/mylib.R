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