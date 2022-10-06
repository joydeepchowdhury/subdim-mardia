# subdim-mardia

 The `Rcpp` file `SkewKurt.cpp` contains the functions **asymskewness**, **asymkurtosis**, **asymskewness_q** and **asymkurtosis_q**, which compute the p-values of the tests described in the paper [Sub-Dimensional Mardia Measures of Multivariate Skewness and Kurtosis](https://www.sciencedirect.com/science/article/pii/S0047259X22000859) [[journal link]](https://www.sciencedirect.com/science/article/pii/S0047259X22000859) [[arXiv link]](https://arxiv.org/abs/2111.14441). **asymskewness** and **asymkurtosis** compute the p-values of the tests of skewness and kurtosis, respectively, based on the complete vectors of all subdimensional Mardia measures of skewness and kurtosis. On the other hand, **asymskewness_q** and **asymkurtosis_q** compute the p-values of the tests of skewness and kurtosis based on the respective vectors of only the q-dimensional Mardia measures of skewness and kurtosis. The number of simulated values from the asymptotic null distribution computed to obtain the p-values in each of the tests is fixed to 1000.
 
 Both the functions **asymskewness** and **asymkurtosis** have a single argument:
 
   - `Data`: An n-by-p matrix, where each row represents a p-dimensional multivariate observation.

On the other hand, both the functions **asymskewness_q** and **asymkurtosis_q** have two single arguments:

   - `Data`: An n-by-p matrix, where each row represents a p-dimensional multivariate observation,
   - `q`: A positive integer giving the value of the dimension q corresponding to the Mardia measures of skewness and kurtosis used to form the vectors described above.

The output of each of the four functions described above is the corresponding p-value.

The code of the file `SkewKurt.cpp` is based on `RcppArmadillo`. To use the functions described above, one has to make a `R` package using `RcppArmadillo`. To make this package, one may use the facility provided in _RStudio_ to make new packages. There, choose the option `R package with RcppArmadillo`. Afterwards, put the file `SkewKurt.cpp` in the src folder of this package, and then build the package using RStudio. After the package is built, call the functions like any other `R` packages.
