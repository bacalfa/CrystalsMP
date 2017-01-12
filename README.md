#1. Scripts and Data Sets
The script and data set files and their description are given below:

- Directory [electronicprop](electronicprop/)
  - File [MetalOxidesRegressionMPI.R](electronicprop/MetalOxidesRegressionMPI.R)
    - R script to obtain kernel regression models (uses `npRmpi` package)
  - File [metal_oxides.txt](electronicprop/metal_oxides.txt)
    - Data set of electronic properties of metal oxides
  - File [ExhaustiveEnumerationSearch.R](electronicprop/ExhaustiveEnumerationSearch.R)
    - R script implementing the exhaustive enumeration algorithm 
  - File [npregkernel.R](electronicprop/npregkernel.R)
    -	R script for the manual implementation of the kernel regression function
- Directory [elasticprop](elasticprop/)
  - File [AllElasticityRegressionMPI.R](elasticprop/AllElasticityRegressionMPI.R)
    - R script to obtain kernel regression models (uses `npRmpi` package)
  - File [allelasticity.txt](elasticprop/allelasticity.txt)
    - Data set of elastic properties of several crystals

#2. Optimal Bandwidths#
The best bandwidths obtained with 5 multi-start optimization runs (function `npreg`), i.e., the bandwidths for which the objective function had the lowest value, for each of the two data sets are given in the following files. The first column contains the predictors (sg stands for space group number).

##2.1 Electronic Properties of Metal Oxides##
File [electronicbw.txt](electronicprop/electronicbw.txt).

##2.2 Elastic Properties##
File [elasticbw.txt](elasticprop/elasticbw.txt).
