#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
library("ape")



description <- "
Creates image with phylogenetic tree.
The program tries to adjust the width and height of the image automatically or
you can set them manually.
Input:  Newick tree
Output: Image"

option_list <- list(
  make_option(c("-i", "--input"), action="store",
              help="Input Newick tree"),
  make_option(c("-o", "--output"), action="store",
               help="Output file name [default tree-<input-filename>.png]")
)

opt_parser <- OptionParser(option_list=option_list, description=description)
opt <- parse_args(opt_parser)

if(is.null(opt$i)){
  print_help(opt_parser)
  stop("\ninput Newick tree is required\n")
}

myTree <- ape::read.tree(opt$i)

if(is.null(opt$o)){
  png(paste0("tree-",basename(opt$i),".png"))
} else {
    png(opt$o)
}

plot(myTree)


invisible(dev.off())
