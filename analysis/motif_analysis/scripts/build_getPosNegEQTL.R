library(mashr)
source('cerebellum_second_try.R')

## get gtex color
gtex.color <- read.table(url('https://github.com/stephenslab/gtexresults/raw/master/data/GTExColors.txt'), sep = '\t', comment.char = '')
gtex.color <- gtex.color[c(1:6,9:18,21:23,26:30,32,33,35,36,38:53), 1:2]
tissue.color <- as.character(gtex.color[,2])
ntissue <- nrow(gtex.color)
## END

flash <- readRDS('../../../output/gtex_flash.rds')
mash <- readRDS('../../../output/gtex_mash.rds')
o <- getPosNegEQTL(flash, mash)
for(i in o$pos) {
  mash_plot_meta(mash$posterior, i, color = meta.colors(lines = tissue.color), main = paste0('index = ', i))
}