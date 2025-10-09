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
                    choices = c("rawdata", "simulate", "normalize", 
                                "integrate", "visualize", "metric"),
                    required = TRUE, 
                    help = "Module type: rawdata, simulate, normalize, integrate, vizualize, metric")

#s <- switch(args$what, rawdata = c("mouse_pancreas"),
#                       simulation = c("scdesign3"))
# TODO: add subparser

parser$add_argument("--flavour", 
                    choices = c("mouse_pancreas",                         # raw data
                                "scdesign3",                              # simulate
                                "log1pCP10k", "log1pCPM", "sctransform", # normalize
                                "x",               # integrate
                                "y",               # visualize
                                "xx", "yy", "zz"), # metric
                    required = TRUE, 
                    help = "Module to run: name depends on the 'what'")

parser$add_argument("--params", type = "character", default = "",
                    help = "Optional parameters as free-form text")

parser$add_argument("--verbose", type = "logical", default = TRUE,
                    help = "TRUE/FALSE as to whether to write progress to stdout")

parser$add_argument("--output_dir", "-o", dest="output_dir", type="character",
                    help="output directory where files will be saved", default=getwd(),
                    required = TRUE)

parser$add_argument("--name", "-n", dest="name", type="character", required = TRUE,
                    help="name of the dataset")

parser$add_argument('--rawdata.ad',
                    type="character",
                    help='gz-compressed H5 file containing (raw) data as AnnData')

parser$add_argument('--simulate.ad',
                    type="character",
                    help='gz-compressed H5 file containing (simulated) data as AnnData')

parser$add_argument('--normalize.ad',
                    type="character",
                    help='gz-compressed H5 file containing (normalized) data as AnnData')

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

# infer the current directory (useful for debugging)
cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname( gsub("--file=","",cargs[[m]]) )
message("location: ", run_dir)
message("libPaths: ", paste0(.libPaths(),collapse=";"))
info <- Sys.info()
message("info: ", paste0(names(info),"=",info,collapse=";"))

# source helper functions (n.b.: args$what controls which code to source)
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

# check if implemented: throw error if not; run if so
# n.b.: args$flavour defines what 'main' function to call
fun <- tryCatch(obj <- get(args$flavour), error = function(e) e)
if ( !("error" %in% class(fun)) ) {
    x <- fun(args) # execute function 
    print(x)
} else {
    message('Unimplemented functionality. Exiting.\n') # throw error?
    quit("no", status = 1)
}

write_ad <- function(x, file) {
  if(args$verbose) message(paste("Converting", class(x), "-> AnnData."))
  x.ad <- as_AnnData(x)
  if(args$verbose) message(paste0("Writing: ", file, "."))
  write_h5ad(x.ad, fn, mode = "w", compression = "gzip")
  if(args$verbose) message("Done.")
}

# write to AnnData via anndataR
if (args$what %in% c("rawdata", "simulate", "normalize", "integrate")) {
  # here, always writing data files as AD (HDF5)
  fn <- file.path(args$output_dir, paste0(args$name,"_",args$what, ".ad"))
  write_ad(x, fn)
} 
# write memento about normalization method
if (args$what == "normalize") {
  fn <- file.path(args$output_dir, paste0(args$name,"_",args$what, ".json"))
  write(toJSON(list(normalize=what$flavour)), fn)
} 
if (args$what == "visualize") {
  # 'x' is something here
} else if (args$what == "metric") {
    # 'x' is something here
}
