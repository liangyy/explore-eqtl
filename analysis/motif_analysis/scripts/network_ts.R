# this function returns genes which has tissue-specific TF edge in
# this paper: https://doi.org/10.1016/j.celrep.2017.10.001

# ../../data/STable1_GTEx_TSInfo_Edges.txt is downloaded from https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/gtex-networks


getPosNegEQTL <- function(flash, mash) {
  myfile <- '../../data/STable1_GTEx_TSInfo_Edges.txt'
  temp <- 'network_ts.temp'
  cmd <- paste('cat', myfile, '| grep heart|grep muscle|awk -F"\t" \'{if($4 < 3) print $0}\'|cut -f2|sort|uniq >', temp)
  system(cmd)
  df <- read.table(temp, header = F)
  eqtl <- read.table('../../output/strong_eqtl_names.txt', header = F)
  e <- 1 : nrow(eqtl)
  cond1 <- stringr::str_match(eqtl[, 1], '([A-Za-z0-9]+).[A-Za-z0-9_]+')[, 2]  %in% df[, 1]
  l.five <- flash$EL[, 5]
  l.one <- flash$EL[, 1]
  cond2 <- abs(l.five) > 3
  pos.idx <- e[cond1 & cond2]
  neg.idx <- e[(abs(l.five) < 1 & abs(l.one) > 50)]
  cmd <- paste('rm -f', temp)
  system(cmd)
  return(list(pos = pos.idx, neg = neg.idx))
}
