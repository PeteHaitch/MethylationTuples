MethylationTuples
================================================================================
[![Build Status](https://travis-ci.org/PeteHaitch/MethylationTuples.png?branch=master)](https://travis-ci.org/PeteHaitch/MethylationTuples)

`MethylationTuples` provides tools for analysing, managing and visualising 
methylation patterns at genomic tuples. These include analyses of 
co-methylation and epipolymorphisms.

__This package is in early development__, but it 
can be installed as follows:

```R
# Install devtools if not already installed
if (suppressWarnings(!require(devtools))) {
  install.packages('devtools')  
}
# Install development version of GenomicTuples
devtools::install_github("PeteHaitch/GenomicTuples")
# Install development version of MethylationTuples
devtools::install_github("PeteHaitch/MethylationTuples")
```
