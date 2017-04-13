# Introduction
This repo will document and serve as code base for testing imputation accuracy of alternative imputation software like mimimac 3, Beagle etc.
This repo will also showcase the simple workflow one need to implement to impute with mentioned software and data formats.

## Links
[minimac3](http://genome.sph.umich.edu/wiki/Minimac3)

[Beagle4.1](https://faculty.washington.edu/browning/beagle/beagle.html)

[Shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

## 2do
* Link to relevant raw data for testing. Tim suggests Affy 50k, Illumina HD 777k and the 153 full sequence data set.
* Implement: basic QC and cleaning pipeline. 
* Implement: plink --> ConformGT strand flip --> exclude common animals --> Impute
* Implement: N-fold cross validation code for Beagle and alternative. 
* Develop R-scrips to calculate correlation R^2 between true and imputed genotypes. 

## Further thougts
@argju Get some sleep, and please fill in.
