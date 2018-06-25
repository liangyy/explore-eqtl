library(optparse)

option_list = list(
  make_option(c("-e", "--eqtl_txt"), type="character", default=NULL,
              help="text file containing meta data of eQTL (matched rows)", metavar="character"),
	make_option(c("-l", "--list_rds"), type="character",
              help="eQTL list RDS", metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character",
              help="prefix of output", metavar="character"),
  make_option(c("-w", "--window_size"), type="numeric",
              help="extended size upstream of TSS", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source('../../scripts/mylib.R')
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(seqinr)

eqtl_annotated <- read.table(opt$eqtl_txt, header = T)
out <- readRDS(opt$list_rds)
pos.idx <- out$pos
neg.idx <- out$neg

df.sub <- getSubsetByIndex(eqtl_annotated, c(pos.idx, neg.idx))
df.sub$specific <- TRUE
df.sub$specific[df.sub$idx %in% neg.idx] <- FALSE
df.sub$tss <- df.sub$transcript_start
df.sub$tss[df.sub$strand == -1] <- df.sub$transcript_end[df.sub$strand == -1]
tss_up1kb <- df.sub$tss - df.sub$strand * opt$window_size
df.sub$tss_big <- df.sub$tss
df.sub$tss_big[df.sub$strand == -1] <- tss_up1kb[df.sub$strand == -1]
df.sub$tss_small <- tss_up1kb
df.sub$tss_small[df.sub$strand == -1] <- df.sub$tss[df.sub$strand == -1]
e <- df.sub %>% group_by(ensembl_gene_id) %>% summarise(union = getUnion(tss_small, tss_big), idx = unique(idx))
e <- e[order(e$idx), ]
df.sub <- df.sub[order(df.sub$idx), ]
df.sub <- df.sub[!duplicated(df.sub$ensembl_gene_id), ]
df.sub[, 'union'] <- e[, 'union']
df.sub <- expendUnion(df.sub)
df.sub$seq <- getFasta(data.frame(paste0('chr', df.sub$chromosome_name), df.sub$interval.start, df.sub$interval.end), BSgenome.Hsapiens.UCSC.hg19)
df.sub.pos <- df.sub[df.sub$specific, ]
write.fasta(as.list(df.sub.pos$seq), paste0(df.sub.pos$idx, ':', df.sub.pos$chromosome_name, ':', df.sub.pos$interval.start, '-', df.sub.pos$interval.end), file.out = paste0(opt$output_prefix, '/positive.fa'))
df.sub.neg <- df.sub[!df.sub$specific, ]
write.fasta(as.list(df.sub.neg$seq), paste0(df.sub.neg$idx, ':', df.sub.neg$chromosome_name, ':', df.sub.neg$interval.start, '-', df.sub.neg$interval.end), file.out = paste0(opt$output_prefix, '/negative.fa'))
