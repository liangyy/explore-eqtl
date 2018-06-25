library(optparse)

option_list = list(
  make_option(c("-f", "--func"), type="character", default=NULL,
              help="name of file used for source", metavar="character"),
	make_option(c("-m", "--mash"), type="character", default="out.txt",
              help="mash RDS", metavar="character"),
  make_option(c("-l", "--flash"), type="character", default="out.txt",
              help="flash RDS", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.txt",
              help="name of output eQTL index list", metavar="character")
);

args_parser = OptionParser(option_list=option_list);
opt = parse_args(args_parser);
source(opt$func)
flash <- readRDS(opt$flash)
mash <- readRDS(opt$mash)
out <- getPosNegEQTL(flash, mash)
saveRDS(out, file = opt$output)
