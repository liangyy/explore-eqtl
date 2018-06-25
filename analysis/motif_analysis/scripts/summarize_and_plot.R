library(optparse)

option_list = list(
  make_option("--inputs", type = "character", help = "list of compute_ld outputs"),
  make_option("--dbnames", type = "character", help = "list of names of the motif databases used"),
  make_option("--dbpaths", type = "character", help = "list of paths of the motif databases used"),
  make_option("--mash", type = "character", help = "mash RDS"),
  make_option("--eqtl", type = "character", help = "strong eQTL file containing meta-data"),
  make_option("--out_motif", type = "character", help = "output directory for motif logo plots"),
  make_option("--out_mash", type = "character", help = "output directory for mash forest plots"),
  make_option("--out_table", type = "character", help = "output filename for summarized table"),
  make_option("--fimo_prefix", type = "character", help = "prefix of fimo TSV")
);

args_parser = OptionParser(option_list=option_list);
opt = parse_args(args_parser);

library(mashr)
library(rmeta)

opt$inputs <- strsplit(opt$inputs, ',')[[1]]
opt$dbnames <- strsplit(opt$dbnames, ',')[[1]]
opt$dbpaths <- strsplit(opt$dbpaths, ',')[[1]]

## loading gtex color
gtex.color <- read.table(url('https://github.com/stephenslab/gtexresults/raw/master/data/GTExColors.txt'), sep = '\t', comment.char = '')
gtex.color <- gtex.color[c(1:6,9:18,21:23,26:30,32,33,35,36,38:53), 1:2]
tissue.color <- as.character(gtex.color[,2])
ntissue <- nrow(gtex.color)
## END


eqtl <- read.table(opt$eqtl, header = T)
eqtl <- eqtl[!duplicated(eqtl$ensembl_gene_id), ]
mash <- readRDS(opt$mash)

motifdb <- list()
for(i in 1 : length(opt$dbnames)) {
  motifdb[[opt$dbnames[i]]] <- opt$dbpaths[i]
}

df.all <- data.frame()

for(f in opt$inputs) {
  info <- strsplit(f, '/')[[1]]
  tissue <- info[4]
  filename <- info[5]
  info2 <- strsplit(filename, '.', fixed = T)[[1]]
  from_info <- info2[2]
  db_str <- info2[3]
  if(from_info == 'from_enrichment') {
    db <- db_str
  } else if(from_info == 'from_predefined_list') {
    db <- strsplit(db_str, '__', fixed = T)[[1]][2]
  }
  fimo.df <- try(read.table(paste0(opt$fimo_prefix, '/', from_info, '/', db_str, '/fimo.tsv'), header = T, sep = '\t'))
  if(class(fimo.df) == 'try-error') {
    next
  }
  motif_snp.df <- try(read.table(f, header = T))
  if(class(motif_snp.df) == 'try-error') {
    next
  }
  motif_id <- motif_snp.df$motif_id
  motif_snp.df <- motif_snp.df[, -which(colnames(motif_snp.df) %in% 'motif_id')]
  motif_id <- paste(motif_snp.df$region_id, motif_id, motif_snp.df$relative_start, motif_snp.df$relative_end, motif_snp.df$relative_strand)
  fimo_motif_id <- paste(fimo.df$sequence_name, fimo.df$motif_id, fimo.df$start, fimo.df$stop, fimo.df$strand)
  motif_snp.df <- cbind(motif_snp.df, fimo.df[match(motif_id, fimo_motif_id), c('motif_id', 'motif_alt_id', 'p.value', 'matched_sequence', 'q.value', 'start', 'stop', 'strand')])
  motif_snp.df$motif_snp.pos <- motif_snp.df$motif_snp_start - motif_snp.df$motif_start + 1
  motif_snp.df$motif.length <- motif_snp.df$stop - motif_snp.df$start + 1
  motif_snp.df$motif_snp.pos[motif_snp.df$strand == '-'] <- motif_snp.df$motif_end[motif_snp.df$strand == '-'] - motif_snp.df$motif_snp_end[motif_snp.df$strand == '-'] + 1
  motif_snp.df$hgnc_symbol <- eqtl[match(motif_snp.df$idx, eqtl$idx), 'hgnc_symbol']
  motif_snp.df$ensembl_gene_id <- eqtl[match(motif_snp.df$idx, eqtl$idx), 'ensembl_gene_id']
  df <- motif_snp.df[, c('hgnc_symbol', 'ensembl_gene_id', 'strong_eqtl_id', 'motif_snp_id', 'r2', 'motif_id', 'motif_alt_id', 'q.value', 'matched_sequence', 'motif_snp.pos', 'motif.length', 'idx', 'strand')]
  colnames(df) <- c('gene_name', 'ensembl_gene_id', 'leading_eqtl', 'motif_snp', 'r2', 'motif_id', 'motif_name', 'motif_qval', 'motif_matched_sequence', 'motif_snp.pos', 'motif.length', 'idx', 'strand')
  df$tissue <- tissue
  df$motif_from <- from_info
  df$dbname <- db
  # print(head(df))
  df.all <- rbind(df.all, df)
}

write.table(df.all, file = opt$out_table, sep = '\t', quote = F, row.names = F)

# print(df.all)
motif.df <- df.all[!duplicated(df.all[, c('motif_id', 'dbname')]), ]

for(i in 1 : nrow(motif.df)) {
  ofile <- paste0(opt$out_motif, '/', motif.df[i, 'dbname'], '__', motif.df[i, 'motif_id'], '.png')
  cmd <- paste0('ceqlogo -i ', motifdb[[motif.df[i, 'dbname']]], ' -m ', motif.df[i, 'motif_id'], ' -o ', ofile, ' -f PNG')
  print(cmd)
  system(cmd)
}

eqtl.df <- df.all[!duplicated(df.all$idx), ]

for(i in 1 : nrow(eqtl.df)) {
  ofile <- paste0(opt$out_mash, '/mash__idx', eqtl.df$idx[i], '.png')
  png(ofile)
  mash_plot_meta(mash$posterior, eqtl.df$idx[i], color = meta.colors(lines = tissue.color))
  dev.off()
}
