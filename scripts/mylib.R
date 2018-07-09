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

getFasta <- function(df, obj) {
  e <- getSeq(obj, df[, 1], start = df[, 2], end = df[, 3], as.character = TRUE)
  return(e)
}

getPosition <- function(pos, repos) {
  if(is.null(nrow(repos))) {
    return(data.frame(start = pos[1] + repos[1] - 1, end = pos[1] + repos[2])) 
  } else {
    n <- nrow(repos)
    out <- data.frame(start = rep(NA, n), end = rep(NA, n), strand = rep(NA, n))
    for(i in 1 : n) {
      out[i, 1] <- pos[1] + repos[i, 1] - 1
      out[i, 2] <- pos[1] + repos[i, 2]
    }
    return(out)
  }
}

getLDmat <- function(chr, pos1, pos2, extend.size) {
  start <- min(pos1, pos2)
  end <- max(pos1, pos2)
  query <- paste0(chr, ':', start - extend.size, '-', end + extend.size)
  cmd <- 'mkdir temp/'
  system(cmd)
  cmd <- paste0('tabix -fh ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz ', query, ' > temp/genotypes.vcf')
  system(cmd)
  cmd <- "vcftools --vcf temp/genotypes.vcf --keep <(wget -qO- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel | grep 'CEU\\|TSI\\|FIN\\|GBR\\|IBS'|cut -f1) --recode > temp/genotypes_eur.vcf --out temp/genotypes_eur"
  system(paste0('echo \"', cmd, '\"'))
  system(paste0('echo \"', cmd, '\" | bash'))
  cmd <- '/Users/yanyu_liang/Documents/programs/bin/plink2 --vcf temp/genotypes_eur.recode.vcf --r2 --out temp/ld_out --ld-window-kb 10000000 --ld-window 1000000 --ld-window-r2 0'
  system(paste0('echo \"', cmd, '\" | bash'))
  df <- read.table('temp/ld_out.ld', header = T)
  cmd <- 'rm -r temp/'
  system(cmd)
  cmd <- 'rm -f ALL.2of4intersection.20100804.genotypes.vcf.gz.tbi'
  system(cmd)
  return(df)
}

getUnion <- function(s, e) {
  gr <- GRanges(rep('a', length(s)), IRanges(s, e), strand="*") 
  gr <- as.data.frame(reduce(gr))
  out <- paste0(gr$start, '-', gr$end)
  return(paste(out, collapse = '|'))
}

expendUnion <- function(df) {
  df.new <- data.frame()
  for(i in 1 : nrow(df)) {
    e <- strsplit(df$union[i], '|', fixed = T)[[1]]
    e <- unlist(strsplit(e, '-'))
    n <- length(e) / 2
    s <- as.numeric(e[seq(1, 2 * n, by = 2)])
    e <- as.numeric(e[seq(2, 2 * n, by = 2)])
    temp <- data.frame(c(df[i, ]))
    temp <- temp[rep(seq_len(nrow(temp)), each=n),]
    temp$interval.start <- s
    temp$interval.end <- e
    df.new <- rbind(df.new, temp)
  }
  return(df.new)
}

getSubsetByIndex <- function(df, idx) {
  df.sub <- df[df$idx %in% idx, ]
  return(df.sub)
}

subchunkify <- function(g, fig_asp = 1) {
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')
  
  sub_chunk <- paste0("`","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.asp=", fig_asp, ", echo=FALSE}",
                      "\n(",
                      g_deparsed
                      , ")()",
                      "\n`","``")
  
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}

annotateEQTL <- function(eqtl, idx, ensembl) {
  ensg <- data.frame(t(sapply(eqtl[idx], function(str) {
    temp <- strsplit(as.character(str), '_')[[1]]
    x <- temp[1]
    chr <- temp[2]
    s <- temp[3]
    return(c(strsplit(x, '.', fixed = T)[[1]][1], chr, s))
  })))
  colnames(ensg) <- c('ensembl_gene_id', 'eqtl_chr', 'eqtl_start')
  ensg$idx <- idx
  ensg$ensembl_gene_id <- as.character(ensg$ensembl_gene_id)
  re <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position', 'strand', 'go_id', 'name_1006', 'uniprot_genename', 'uniprot_swissprot_accession', 'transcript_start', 'transcript_end'), 
    filters = 'ensembl_gene_id',
    values = ensg[, 1],
    mart = ensembl)
  re <- inner_join(re, ensg, by = 'ensembl_gene_id')
  return(re)
}
genGene2GO <- function(g2g, re) {
  gene <- unique(re$ensembl_gene_id)
  for(x in gene) {
    g2g[[x]] <- re[re$ensembl_gene_id == x, 'go_id']
  }
  return(g2g)
}