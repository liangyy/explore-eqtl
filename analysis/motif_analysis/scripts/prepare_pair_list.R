library(optparse)

option_list = list(
  make_option(c("-m", "--motif_snp"), type="character", default=NULL,
              help="list of motif variant in genotype data", metavar="character"),
	make_option(c("-s", "--eqtl_name"), type="character", default="out.txt",
              help="eQTL meta-data containing strong/leading eQTL", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.txt",
              help="name of output LD", metavar="character")
);

args_parser = OptionParser(option_list=option_list);
opt = parse_args(args_parser);

library(stringr)

motif.snp <- try(read.table(opt$motif_snp))
if(class(motif.snp) == 'try-error') {
  write.table(data.frame(), file = opt$output, row.names = F, col.names = F, sep = '\t', quote = F)
  quit()
}
strong.eqtl <- read.table(opt$eqtl_name)
motif.snp$idx <- as.numeric(str_match(as.character(motif.snp$V11), '^([0-9]+):')[, 2])
motif.snp$strong_eqtl_meta <- strong.eqtl[motif.snp$idx, 1]
motif.snp$strong_eqtl <- str_match(motif.snp$strong_eqtl_meta, '^[A-Za-z0-9.]+_([A-Za-z0-9_]+)$')[, 2]
# out <- data.frame(snp1 = motif.snp$V6, snp2 = motif.snp$strong_eqtl)
colnames(motif.snp) <- c('motif_snp_chr', 'motif_snp_start', 'motif_snp_end', 'motif_snp_a1', 'motif_snp_a2', 'motif_snp_id', 'motif_chr', 'motif_start', 'motif_end', 'strand_null', 'region_id', 'motif_id', 'relative_start', 'relative_end', 'relative_strand', 'idx', 'strong_eqtl_name', 'strong_eqtl_id')
write.table(motif.snp, file = opt$output, row.names = F, col.names = T, sep = '\t', quote = F)
