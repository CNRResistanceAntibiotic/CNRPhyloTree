#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
library("ape")
library("Cairo")
library("showtext", quietly = T)
library("phangorn")
library("phytools", quietly = T)



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
               help="Output file name [default tree-<input-filename>.png]"),
  make_option(c("--width"), action="store", type="double", default=F,
              help="Sets the width in inches of the image"),
  make_option(c("--height"), action="store", type="double", default=F,
              help="Sets the height in inches of the image")
)

opt_parser <- OptionParser(option_list=option_list, description=description)
opt <- parse_args(opt_parser)

if(is.null(opt$i)){
  print_help(opt_parser)
  stop("\ninput Newick tree is required\n")
}

font.add("Arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")

if(!is.null(opt$i)){
  tree <- ape::read.tree(opt$i)
  if(!(is.double(opt$height) & opt$height > 0)){
    # Calculate image height
    nodes <- ape::Nnode(tree)
    nh = 10  # height of individual node (in pixels)
    c1 = 4   # distance between nodes
    res = 96 # resolution
    if (nodes <= 10) {
      c2 = 100 # additional space for plots with few nodes (empirically based)
    } else {
      c2 = 20  # additional space
    }
    opt$height <- ((nodes * (nh + c1)) + c2) / res
  }
  if(!(is.double(opt$width) & opt$width > 0)){
    opt$width <- 3 # Constant image width with x.lim = c(0,0.3)
  }
}

if(is.null(opt$o)){
  png(paste0("tree-",basename(opt$i),".png"),
      type = "cairo",
      antialias = "none",
      res = 96,
      width = as.numeric(opt$width),
      height = as.numeric(opt$height),
      units = 'in')
} else {
    png(opt$o,
      type = "cairo",
      antialias = "none",
      res = 96,
      width = as.numeric(opt$width),
      height = as.numeric(opt$height),
      units = 'in')
}

showtext.begin()

par(family = "Arial")

tree <- phytools::rotateNodes(tree, nodes = "all")
tree <- phangorn::midpoint(tree)
tree <- ape::reorder.phylo(tree)

ape::plot.phylo(tree,
                no.margin = T,
                label.offset = 0.01,
                align.tip.label = T,
                cex = 0.8,
                x.lim = c(0,0.3))

ape::add.scale.bar(length = 0.02, cex = 0.7)

showtext.end()

invisible(dev.off())