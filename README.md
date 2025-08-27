# vizbench-r

This repo is where the code to call ALL R functions/packages is for the vizbench benchmark (designed by Jessica Jingyi Li's group and in particular, Zhiqian Zai and Qingyang Wang). For the omnibenchmark side, it's organized as a master repo that handles all (R) aspects, from deriving the datasets (usually, from public repos), running preprocessing methods (normalization, integration), and running visualization methods and then computing metrics.

TODO: add some docs here to this page to highlight how all the arguments are handled/designed.

The script will have arguments as follows:

```
--what     - specifies which operation will be performed (normalization, visualization, etc.) 
--flavour  - which data/method/metric is run (conditional on 'what')
--params   - additional parameters that are optionally needed
--verbose  - whether aspects of the output are not suppressed (default: TRUE)
```

The idea is also that the script will source some R scripts in the `utils` directory of the repo, one for each `what` specified (i.e. `file.path("utils", paste0(args$what, "_utils.R"))`, containing the relevant code for the `flavour` of interest.
