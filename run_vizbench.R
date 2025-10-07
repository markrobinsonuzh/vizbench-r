#!/usr/bin/env Rscript

## Originally sketched by Izaskun Mallona
## modified/populated by Mark Robinson
## various chunks of code from Zhiqian Zhai and Qingyang Wang will be linked
## from original repo: https://github.com/zhiqianZ/Benchmark-Normalization-Integration-Visualization

## Usage:
## to do a system call to showcase how to call other subscripts
##    Rscript thisfilename.R --what type --flavour specific-module

library(argparse)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

# define arguments
parser$add_argument("--what", 
                    choices = c("rawdata", "simulation"),
                    # choices = c("rawdata", "simulation", "normalization", "integration", "metric"), 
                    required = TRUE, 
                    help = "Module type: rawdata, simulation, normalization, integration, viz, metric")

#s <- switch(args$what, rawdata = c("mouse_pancreas"),
#                       simulation = c("scdesign3"))
# TODO: add subparser

parser$add_argument("--flavour", choices = c("mouse_pancreas", "scdesign3"),
                    required = TRUE, help = "Module to run: name depends on the 'what'")

parser$add_argument("--params", type = "character", default = "",
                    help = "Optional parameters as free-form text")

parser$add_argument("--verbose", type = "logical", default = TRUE,
                    help = "Optional parameters as free-form text")

parser$add_argument("--output_dir", "-o", dest="output_dir", type="character",
                    help="output directory where files will be saved", default=getwd(),
                    required = TRUE)
parser$add_argument("--name", "-n", dest="name", type="character", required = TRUE,
                    help="name of the dataset")

parser$add_argument('--rawdata.ad',
                    type="character",
                    help='gz-compressed H5 file containing data as AnnData')
# parser$add_argument('--data.true_labels',
#                     type="character",
#                     help='gz-compressed textfile with the true labels; used to select a range of ks.')


# send details to be logged
args <- parser$parse_args()
message("Selected category: ", args$what)
message("Routine selected: ", args$flavour)
message("Additional parameters: ", args$params)
message("name: ", args$name)
message("Verbose: ", args$verbose)

# infer the current directory
cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname( gsub("--file=","",cargs[[m]]) )
message("location: ", run_dir)
message("libPaths: ", paste0(.libPaths(),collapse=";"))
info <- Sys.info()
message("info: ", paste0(names(info),"=",info,collapse=";"))

# source helper functions
helpers <- file.path(run_dir, "utils", paste0(args$what, "_utils.R"))
if( file.exists(helpers) ) {
    message("Sourcing .. ", helpers)
    source(helpers)
} else {
    message(paste0("Helper code in ",helpers," not found. Exiting."))
    quit("no", status = 1)
}

# load packages
suppressPackageStartupMessages(load_pkgs())
#if(args$verbose) load_pkgs() else suppressPackageStartupMessages(load_pkgs())

# check if implemented: throw error if not; run if so
fun <- tryCatch(obj <- get(args$flavour), error = function(e) e)
if ( !("error" %in% class(fun)) ) {
    x <- fun(args) # execute function 
    print(x)
} else {
    message('Unimplemented functionality. Exiting.\n') # throw error?
    quit("no", status = 1)
}

if (args$what == 'rawdata') {
  # 'x' should be a SingleCellExperiment object
  fn <- file.path(args$output_dir, paste0(args$name, ".ad"))
  message(paste("Converting SCE -> AnnData."))
  x.ad <- as_AnnData(x)
  message(paste0("Writing: ", fn, "."))
  write_h5ad(x.ad, fn, mode = "w", compression = "gzip")
  message("Done.")
} else if (args$what == 'simulation') {
  fn <- file.path(args$output_dir, paste0(args$name, ".ad"))
  message(paste("Converting Seurat -> AnnData."))
  x.ad <- as_AnnData(x)
  message(paste0("Writing: ", fn, "."))
  write_h5ad(x.ad, fn, mode = "w", compression = "gzip")
  message("Done.")
  # 'x' should be a Seurat object now
} else if (args$what == 'method') {
    # 'x' is something here
} else if (args$what == 'metric') {
    # 'x' is something here
}
