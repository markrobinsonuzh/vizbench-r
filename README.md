# vizbench-r

This repo is where the code to call ALL R functions/packages is for the vizbench benchmark (designed by Jessica Jingyi Li's group and in particular, Zhiqian Zai and Qingyang Wang). For the omnibenchmark side, it's organized as a master repo that handles all (R) aspects, from deriving the datasets (usually, from public repos), running preprocessing methods (normalization, integration), and running visualization methods and then computing metrics.

TODO: add some docs here to this page to highlight how all the arguments are handled/designed.

The script will have arguments as follows:

```
--what     - specifies which operation will be performed (rawdata, simulation, normalization, 
--flavour  - which data/method/metric is run
--params   - additional parameters that are optionally needed
--verbose  - whether some aspects of the output are suppressed (default: FALSE)
```

The idea is also that the script will source some R scripts in the `utils` directory, one for each `what` specified, containing the relevant code for the `flavour` of interest.
