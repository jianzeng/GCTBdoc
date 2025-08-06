
## Basic options {: .expand}

### Input and output

**\--bfile** test

Input PLINK binary PED files, e.g. test.fam, test.bim and test.bed (see PLINK user manual for details).

**\--pheno** test.phen

Input phenotype data from a plain text file, e.g. test.phen.

**\--out** test

Specify output root filename.

### Data management

**\--keep** test.indi.list

Specify a list of individuals to be included in the analysis.

**\--chr** 1

Include SNPs on a specific chromosome in the analysis, e.g. chromosome 1.

**\--extract** test.snplist

Specify a list of SNPs to be included in the analysis.

**\--exclude** test.snplist

Specify a list of SNPs to be excluded from the analysis.

**\--mpheno** 2

If the phenotype file contains more than one trait, by default, GCTB takes the first trait for analysis (the third column of the file) unless this option is specified. For example, **\--mpheno** 2 tells GCTB to take the second trait for analysis (the fourth column of the file).

**\--covar** test.qcovar

Input quantitative covariates from a plain text file, e.g. test.qcovar. Each quantitative covariate is recognized as a continuous variable.

**\--random-covar** test.randcovar

Input quantitative covariates from a plain text file which will be fitted as random effects in the model. For a categorical variable with k levels, create a matrix of 0/1 with k columns to indicate the presence of each level in a column.

### MCMC settings

**\--seed** 123

Specify the seed for random number generation, e.g. 123. Note that giving the same seed value would result in exactly the same results between two runs.

**\--chain-length** 3000

Specify the total number of iterations in MCMC, e.g. 3000 (default).

**\--burn-in** 1000

Specify the number of iterations to be discarded, e.g. 1000 (default).

**\--out-freq** 100

Display the intermediate results for every 100 iterations (default). 

**\--thin** 10

Output the sampled values for SNP effects and genetic architecture parameters for every 10 iterations (default). Only non-zero sampled values of SNP effects are written into a binary file.

**\--no-mcmc-bin**

Suppress the output of MCMC samples of SNP effects.
