library(optparse)

option_list = list(
  make_option(c("-x", "--x"), type="character", default=NULL,
              help="name of directory for genotype files (by gene)", metavar="character"),
	make_option(c("-y", "--y"), type="character", default="out.txt",
              help="name of directory for phenotype files (by gene)", metavar="character"),
  make_option(c("-z", "--z"), type="character", default="out.txt",
              help="name of the file for covariate", metavar="character"),
  make_option(c("-g", "--gene_list"), type="character", default="out.txt",
              help="a list of gene names separated by ','", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default="out.txt",
              help="name of output directory", metavar="character"),
  make_option(c("-f", "--func"), type="character", default="out.txt",
              help="fine mapping function to use", metavar="character")
);

# x: 1. header: individual ID; 2. one individual per column; 3. first column: SNP ID
# y: 1. header: individual ID; 2. one individual per column; 3. first column: gene ID
# z: 1. header: individual ID; 2. one individual per column; 3. first column: covariate ID

args_parser = OptionParser(option_list=option_list);
opt = parse_args(args_parser);
source(opt$func) # obtain fine-mapping function in opt$func: method = function(x, y, z) { # some code; return(obj)}

z <- opt$z
genes <- strsplit(opt$gene, ',')[[1]]
for(g in genes) {
  print(g)
  x <- paste0(opt$x, '/', g, '.txt.gz')
  y <- paste0(opt$y, '/', g, '.txt.gz')
  o <- paste0(opt$out_dir, '/', g, '.rds')
  if(g == '') {
    next
  }
  geno.df <- read.table(x, header = T)
  expr.df <- read.table(y, header = T)
  cova.df <- read.table(z, header = T)
  geno <- geno.df[, 2 : ncol(geno.df)]
  expr <- expr.df[, 2 : ncol(expr.df)]
  cova <- cova.df[, 2 : ncol(cova.df)]
  indv <- intersect(intersect(colnames(geno), colnames(expr)), colnames(cova))
  indv <- indv[order(indv)]
  geno.cleaned <- geno[, match(indv, colnames(geno))]
  expr.cleaned <- expr[, match(indv, colnames(expr))]
  cova.cleaned <- cova[, match(indv, colnames(cova))]
  out <- fineMapping(t(data.matrix(geno.cleaned)), as.numeric(expr.cleaned[1, ]), t(data.matrix(cova.cleaned)))
  geno.cleaned$SNP_ID <- geno.df[, 1]
  saveRDS(list(fitted = out, X = geno.cleaned), file = o)
}
