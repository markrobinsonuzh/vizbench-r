#!/usr/bin/env Rscript
##
## Usage:
## to do a system call to showcase this could be used to call other scripts (well, system('date'))
##    Rscript thisfiliename.R --what method --routine date 
## to do a simulation using some freetext params
##    Rscript thisfiliename.R --what dataset --routine simulate --params "10, mean=5, sd=2"

library(argparse)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

# define arguments; params could be use to simplify passing extra free text
parser$add_argument("--what", choices = c("method", "metric", "dataset"), required = TRUE,
                    help = "Type of target: method, metric, or dataset")

parser$add_argument("--routine", choices = c("simulate", "edger", "limma", "date"), required = TRUE,
                    help = "Benchmarking routine: simulate, edger, or limma")

parser$add_argument("--params", type = "character", default = "",
                    help = "Optional parameters as free-form text")


args <- parser$parse_args()

## cat("Selected category: ", args$what, "\n")
## cat("Routine selected: ", args$routine, "\n")
## cat("Additional parameters: ", args$params, "\n")

if (args$what == 'dataset') {
    if (args$routine == 'simulate')  {
        eval_expr <- parse(text = paste0("rnorm(", args$params, ")"))
        print(eval(eval_expr))
    }
    else {
        cat('unimplemented\n')
    }
} else if (args$what == 'method') {
    if  (args$routine == 'date') {
        system('date')
    }
    else {
        cat('unimplemented\n')
    }
} else if (args$what == 'metric') {
    cat('unimplemented')
}
