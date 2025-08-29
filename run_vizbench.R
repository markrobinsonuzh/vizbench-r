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
                    help = "Module type: rawdata, simulation, normalizaiton, integration, viz, metric")

#s <- switch(args$what, rawdata = c("mouse_pancreas"),
#                       simulation = c("scdesign3"))
#cat(s, "\n")
# TODO: add subparser via 
# subparsers = parser.add_subparsers (Python)
# see https://stackoverflow.com/questions/9505898/conditional-command-line-arguments-in-python-using-argparse

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

args <- parser$parse_args()
message("Selected category: ", args$what)
message("Routine selected: ", args$flavour)
message("Additional parameters: ", args$params)
message("Verbose: ", args$verbose)

# source helper functions
helpers <- file.path("utils", paste0(args$what, "_utils.R"))
if( file.exists(helpers) ) {
    message("Sourcing .. ", helpers)
    source(helpers)
} else {
    message(paste0("Helper code in ",helpers," not found. Exiting."))
    quit("no", status = 1)
}

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
    # TODO: write it out in non-proprietary format?
    # try to match this:
    #     outputs:
    #   - id: rawdata.sce
    #    path: "{input}/{stage}/{module}/{params}/{dataset}.sce.rds"
    out_fn <- file.path(args$output_dir, paste0(args$name,".sce.rds"))
    message(paste0("Writing to ", out_fn, " .."))
    saveRDS(obj, out_fn)
    message("Done.")
} else if (args$what == 'simulation') {
    # 'x' is something here
} else if (args$what == 'simulation') {
    # 'x' is something here
} else if (args$what == 'metric') {
    # 'x' is something here
}
