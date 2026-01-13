## Tutorial

### SBayesR Tutorial

This updated version of the GCTB software (version 2.0) includes summary-data based versions of the individual-data Bayesian linear mixed models 
previously implemented. These methods require summary statistics from genome-wide association studies, which typically 
include the estimated univariate effect for each single nucleotide polymorphism (SNP), the standard error of the effect, the sample size 
used to estimate the effect for each SNP, the allele frequency of an allele, and an estimate of LD among SNPs, all of which are easily 
accessible from public databases.

This tutorial will outline how to use the summary data based methods and accompanies the manuscript 
[Improved polygenic prediction by Bayesian multiple regression on summary statistics](https://www.nature.com/articles/s41467-019-12653-0). To recreate the tutorial you will need the [PLINK 2](https://www.cog-genomics.org/plink/2.0/) software and 
the updated version of the [GCTB](https://cnsgenomics.com/software/gctb/) software [compiled](https://cnsgenomics.com/software/gctb/download/README.html) from source or the binary executable in your path. The data for the tutorial are available at [Download](https://cnsgenomics.com/software/gctb/#Download). The tutorial is designed so that each part can be reconstructed or picked up at any point with the
following directory structure. 

```
tutorial
  gctb # Static binary executable for Linux 64-bit systems 
   data # 1000 Genomes data in PLINK format
   pheno_generation # R scripts for generating simulated phenotypes 
   pheno # Simulated phenotypes
   gwas # PLINK GWAS summary statistics
   ldm # LD correlation matrices and scripts to calculate them
   ma # Summary statistics in GCTB compatible format
   sbayesr # Scripts for running and SBayesR analysis using GCTB and subsequent results
```

#### Data

The tutorial will run through all aspects of the software using 
genotype data from chromosome 22 of Phase 3 of the [1000 Genomes Project (1000G)](https://www.nature.com/articles/nature15393). The genotype
data have been filtered to exclude variants with minor allele frequency (MAF) < 0.01, which left 15,938 SNPs available for analysis. 
These data include a subset of 378 individuals of European ancestry from the CEU, TSI, GBR and FIN populations. 

N.B. The results generated from this example tutorial using data from the 1000G are designed to be lightweight and should be only used
to explore how to run summary data based analyses using GCTB. The results do not make the best use of GCTB capabilities given the small sample 
size \\$(N=378)\\$ in the example. 


#### Phenotype simulation and GWAS

Using these genotypes, phenotypes were generated under the multiple regression model \\$y\_i = \sum\_{j=1}^pw\_{ij}\beta\_j + \epsilon\_i\\$, where \\$w\_{ij} = (x\_{ij} − 2q\_j)/ \sqrt{2q\_j(1 − q\_j)}\\$
with \\$x\_{ij}\\$ being the reference allele count for the \\$i\\$th individual at the \\$j\\$th SNP, \\$q\_j\\$ the allele
frequency of the \\$j\\$th SNP and \\$\epsilon\_i\\$ was sampled from a normal distribution with mean 0 and variance 
\\$\text{Var}({\bf W}\boldsymbol\beta)(1/h\_{SNP}^2 − 1)\\$ such that \\$h\_{SNP}^2 = 0.1\\$ for each of 20 simulation replicates, 
which is much larger than the contribution to the genome-wide SNP-based heritability (\\$h^2\_{SNP}\\$) estimate for chromosome 22 
for most quantitative traits. All phenotypes were generated using the R programming language with scripts used available in the
`pheno_generation` subdirectory.

For each scenario replicate, we randomly sampled a new set of 1,500 causal variants. The genetic architecture simulated contained two causal variants of large effect explaining 3% and 2% of the phenotypic variance respectively and a polygenic tail of 1,498 causal variants sampled from a N(0, 0.05/1, 498) distribution such that the expected total genetic variance explained by all variants was 0.1.

For each of the 20 simulation replicates, simple linear regression for each variant was run using the PLINK 2 software to generate summary statistics. Code for running the GWAS and the output is available in the `gwas` subdirectory.

#### GCTB summary statistics input format

The GCTB summary-based methods have inherited the [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO) `.ma` format.

```
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
```

Columns are SNP identifier, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value 
and sample size. The headers are not keywords and will be omitted by the program. Important: "A1" needs to be the effect allele with 
"A2" being the other allele and "freq" should be the frequency of "A1".

The transformation of your summary statistics file to the `.ma` format will depend on the software used to analyse your data. The
 `ma` subdirectory contains a simple R script for constructing `.ma` files from PLINK 2 `--linear` output. There is no need to filter
 you summary statistics to match the LD reference as GCTB will perform this data management for you.

#### GCTB linkage disequilibrium (LD) correlation matrix 

The construction of an LD correlation matrix is a key component of the summary-data based methods in GCTB. Theoretically, 
the summary-data based methods are equivalent to the individual data methods if the LD matrix is constructed from all 
variants and the individuals used in the GWAS. Calculating and storing LD matrices that contain all pairwise entries is computationally
prohibitive and thus a subset of the matrix entries are stored. Typically, all inter-chromosomal LD is ignored and 
a within chromosome block diagonal matrix structure is used. Futhermore, the individual level data from
the GWAS is not available, which is often the case, and a reference population that is similar in ancestry to the GWAS
cohort is used to construct the LD correlation matrix. The use of a subset of the correlation matrix and the construction
of the LD matrix from a reference leads to an approximation to the individual data mode. The level of deterioration in the
model parameter estimates introduced by these approximations is assessed empirically through simulation, real data analysis
and comparison with individual-data method estimates.

#### *Full chromosome-wise LD matrices*

Although possible to create genome-wide LD matrices in GCTB, we recommend to only 
calculate matrices for each chromosome. The following is the basic execution of
GCTB LD matrix calculation. The `ldm` subdirectory has example BASH scripts to
perform this on a HPC system. N.B. the HPC commands may not be compatible with your HPC scheduler 
and will require alteration.

```
# Build full LD matrix for chromosome 22
gctb --bfile 1000G_eur_chr22 --make-full-ldm --out 1000G_eur_chr22
```

GCTB produces a `.bin` file, which contains the matrix elements stored in binary, and a `.info`
file, which looks like

```
Chrom ID         GenPos PhysPos  A1 A2 A2Freq   Index WindStart WindEnd WindSize WindWidth N   SamplVar  LDsum
22    rs131538   0      16871137 A  G  0.066138 0     0         15937   15938    34347869  378 41.92552 10.18005
22    rs9605903  0      17054720 C  T  0.264550 1     0         15937   15938    34347869  378 41.95537 -33.15312
22    rs5746647  0      17057138 G  T  0.052910 2     0         15937   15938    34347869  378 41.93396 15.90974
22    rs16980739 0      17058616 T  C  0.148148 3     0         15937   15938    34347869  378 41.94403 -12.93685
```

The columns headings are chromosome, SNP identifier, genetic position, physical base pair position, SNP allele 1,
SNP allele 2, frequency of the A2 allele, variant index, for this SNP (row) the first column that is non-zero
in the LD matrix, the last column that is non-zero, how many non-zero elements for this row, difference in
base pair position between WindStart and WindEnd, the sample size used
to calculate the correlation, sum of samping variance of the LD correlation estimate for the SNP and all 
other SNPs on the chromosome and the sum of the non-zero LD correlation estimates between this SNP and 
all other SNPs.

#### *Full chromosome-wise LD matrices using multiple CPUs*

To reduce the computational burden of computing LD correlation matrices, GCTB can partition the genotype file
into user-specified chunks. 

```
# Build full LD matrix for chromosome 22 for the first 5000 variants 
gctb --bfile 1000G_eur_chr22 --make-full-ldm --snp 1-5000 --out 1000G_eur_chr22
```

When the `--snp` option is invoked the GCTB software will append `.snp1-5000`, for example, to the output files
so no separate extension to `--out` is required. 

For the 1000G chromsome 22 example the following files are generated

```
1000G_eur_chr22.snp1-5000.ldm.full.bin
1000G_eur_chr22.snp1-5000.ldm.full.info
1000G_eur_chr22.snp5001-10000.ldm.full.bin
1000G_eur_chr22.snp5001-10000.ldm.full.info
1000G_eur_chr22.snp10001-15000.ldm.full.bin
1000G_eur_chr22.snp10001-15000.ldm.full.info
1000G_eur_chr22.snp15001-20000.ldm.full.bin
1000G_eur_chr22.snp15001-20000.ldm.full.info
```

To merge the files, a separate file
with the list of files to be joined is required. For our example is looks like

```
cat 1000G_eur_chr22.mldmlist
1000G_eur_chr22.snp1-5000.ldm.full
1000G_eur_chr22.snp5001-10000.ldm.full
1000G_eur_chr22.snp10001-15000.ldm.full
1000G_eur_chr22.snp15001-20000.ldm.full
```

```
# Merge the LD matrix chunks into a single file
gctb --mldm 1000G_eur_chr22.mldmlist --make-full-ldm --out 1000G_eur_chr22
```

#### *Shrunk LD matrix*

The shrinkage estimator for the LD correlation matrix was originally proposed by [Wen and Stephens (2010)](https://projecteuclid.org/euclid.aoas/1287409368).
The estimator shrinks the off-diagonal entries of 
the sample LD matrix toward zero. [Zhu and Stephens (2017)](https://projecteuclid.org/euclid.aoas/1507168840)
used the shrinkage estimator in their Regression with Summary Statistics (RSS) methodology and showed
empirically that it can provide improved inference. The shrinkage estimator overcomes some of the
issues arising from approximating the full LD matrix in the summary-data based model using a subset 
of the full LD matrix and constructing the matrix from a reference. The GCTB implementation is a C++
port from that provided with the RSS software and has been adapted for use with the GCTB software.

The calculation of the LD correlation matrix shrinkage estimate requires a genetic map. 
The genetic map files can be provided by the user but the tutorial uses interpolated map positions
for the CEU population generated from the 1000G OMNI arrays that were downloaded from
the from https://github.com/joepickrell/1000-genomes-genetic-maps. The file format looks like

```
rs149201999 16050408 0.0
rs146752890 16050612 0.0
rs139377059 16050678 0.0
rs188945759 16050984 0.0
rs6518357   16051107 0.0
rs62224609  16051249 0.0
rs62224610  16051347 0.0
rs143503259 16051453 0.0
rs192339082 16051477 0.0
rs79725552  16051480 0.0
```

and contains SNP identifier, physical base pair position and genetic position (cumulative). 
GCTB can use any genetic map provided it adheres to this format.

In this tutorial, the calculation of the shrunk LD matrix requires the effective population 
sample size, which we set to be 11,400 (as in [Zhu and Stephens (2017)](https://projecteuclid.org/euclid.aoas/1507168840)), 
the sample size of the genetic map reference, which corresponds to the 183 individuals from the CEU 
cohort of the 1000G used to create the genetic map, and the hard threshold on the shrinkage value, 
which we set to \\$10^{−5}\\$ and is the default value. These three parameters can be changed with 

```
--ne            # Effective population size. Default is 11,400
--genmap-n      # Genetic map sample size. Default is 183
--shrunk-cutoff # Shrinkage hard threshold. Default is 10^-5
```

The constuction of the shrunk LD matrix is similar to that for the full LD matrix

```
# Build shrunk LD matrix for chromosome 22 
gctb --bfile 1000G_eur_chr22 \
     --make-shrunk-ldm \
     --gen-map chr22.OMNI.interpolated_genetic_map \
     --out 1000G_eur_chr22
```

The GCTB can also perform shrunk matrix calculation using multiple CPUs with a similar
process to the above calculation of the full LD matrix.

```
# Build shrunk LD matrix for chromosome 22 for the first 5000 variants 
gctb --bfile 1000G_eur_chr22 \
     --make-shrunk-ldm \
     --gen-map chr22.OMNI.interpolated_genetic_map \
     --snp 1-5000 \
     --out 1000G_eur_chr22
```

Please see the `/ldm` subdirectory for example BASH scripts of performing these tasks
both locally and on a HPC. 

#### *Sparse LD matrix*

Motivated by the computational practicality of not storing the genome-wide LD correlation matrix we wish to set a subset of elements of the 
LD correlation to zero. The calculation of an LD correlation matrix from genotype data is an estimate of the true LD correlation matrix for 
a population. It therefore contains a mix of true LD and those values that differ from zero due to sampling variation, which we wish to 
ignore. By default we ignore interchromsomal LD between markers. Within chromosome, we set elements to zero if their chi-squared statistic 
under the sampling distribution of the correlation coefficient does not exceed a user-chosen threshold.

Specifically, from \citet{lynch1998genetics} the sampling variation of the correlation coefficient \\$r = \frac{\text{cov}(x, y)}{[\text{var}(x)\text{var}(y)]^{1/2}}\\$ is 
\\$\sigma^2(r) \approx \frac{(1-\rho^2)^2}{n}\\$, which under the null hypothesis that \\$\rho = 0\\$ is just \\$\frac{1}{n}\\$. Given this we can
generate a chi-squared statistic for each element of the LD matrix by \\$(\frac{r}{\sqrt{1/n}})^2\\$. These statistics should be \\$\chi\_1^2\\$ distributed
and thus for each LD matrix element we set it to zero if the statistic does not exceed a user-specified arbitrary value. The chi-squared
statistic can be altered with `--chisq`. The default value is 10, which we have found empirically to be a good balance balance between
matrix sparsity and retaining true positive LD values.

```
# Build sparse LD matrix for chromosome 22
gctb --bfile 1000G_eur_chr22 --make-sparse-ldm --out 1000G_eur_chr22
# Sparse matrices can also be built using multiple CPUs
gctb --bfile 1000G_eur_chr22 --make-sparse-ldm --snp 1-5000 --out 1000G_eur_chr22
```

GCTB can alter current LD matrices, for example, make a full matrix sparse.

```
# Make a full LD matrix sparse
gctb --ldm 1000G_eur_chr22.ldm.full --make-sparse-ldm --out 1000G_eur_chr22
```

We recommend that the shrunk LD matrices be stored in sparse format for efficient
storage and computation. The shrunk LD correlation matrix estimate will have zeroes
which can be eliminated as follows

```
# Make a shrunk LD matrix sparse by eliminating elements equal to 0. Output sparse
gctb --ldm 1000G_eur_chr22.ldm.shrunk --make-sparse-ldm --chisq 0 --out 1000G_eur_chr22
```

#### *Extra options for LD matrix calculation*

```
# For checking or comparison with other methods you can write the LD correlation out in 
# text format.
--write-ldm-txt # Write LD matrix as text file as well as binary. Will generate large files.
# Examples
gctb --bfile 1000G_eur_chr22 --make-shrunk-ldm --write-ldm-txt --out 1000G_eur_chr22
gctb --bfile 1000G_eur_chr22 --make-full-ldm   --write-ldm-txt --out 1000G_eur_chr22
# Read a list of LD matrix files. This would be used to read the list of chromosome wise LD
# matrices when performing genome-wide analyses. Also used for merging matrix chunks for 
# multiple CPU matrix building
--mldm
```

#### Running an SBayesR analysis

With the GWAS summary statistics in `.ma` format and the LD correlation matrix calculated, an SBayesR analysis can be performed by

```
# Calling SBayesR
gctb --sbayes R 
     --ldm ../ldm/sparse/chr22/1000G_eur_chr22.ldm.sparse
     --pi 0.95,0.02,0.02,0.01
     --gamma 0.0,0.01,0.1,1
     --gwas-summary ../ma/sim_1.ma
     --chain-length 10000
     --burn-in 2000
     --out-freq 10
     --out sim_1
```

`--sbayes R`

Specify the summary data based Bayesian alphabet model for the analysis, e.g., `R`. Different letters of the alphabet launch 
different models, which differ in the prior specification for the SNP effects. `R` and `C` are currently available.

`--pi 0.95,0.02,0.02,0.01`

Specific to the `R` model. A comma seperated string where the number of values defines the number of mixture components and each value defines the starting value for each component (the first value is reserved for the zero component). 
These must sum to 1.
The length of the \\$\pi\\$ and \\$\gamma\\$ vectors stipulates how many components are included in the finite mixture of normals distribution prior 
for the genetic effects. The default values are as specified in the example.

`--gamma 0,0.01,0.1,1`

Specific to the `R` model. Speficies the gamma values seperated by comma, each representing the scaling factor for the 
the variance of the distribution of genetic effects \\$\sigma\_{\beta}^2\\$ for each distribution. Note that the vector length should match 
that in `--pi`. The default values are as specified in the example. These differ from the original BayesR model as the
weights are for the variance of the distribution of genetic effects \\$\sigma\_{\beta}^2\\$ rather than the genetic variance \\$\sigma\_g^2\\$. 

`--gwas-summary`

Path to the GWAS summary statistics in `.ma` format.

The other options are inherited from the individual-data GCTB models.

#### *Standard output*

Upon execution of GCTB, preliminary data management information will be 
outputted to standard out. A key summary of the summary is also 
displayed, for example, 

```
Data summary:
                                             mean       sd
            GWAS SNP Phenotypic variance    1.002    0.052
                    GWAS SNP sample size      378        0
                         GWAS SNP effect    0.000    0.117
                             GWAS SNP SE    0.106    0.049
            MME left-hand-side diagonals  379.478   19.774
                     MME right-hand-side    0.130   19.983
                    LD sampling variance    1.298    0.076
                                LD score   17.359   15.140
```                             

These statistics can be used to diagnose errors between individual runs of the program. 
For each row the mean and standard deviation (sd) of the parameter distribution are displayed. The
acronym MME represents mixed model equations. The left-hand-side diagonals
are the elements of \\${\bf D}=\text{diag}({\bf X'}{\bf X})\\$, where \\$D\_j=2p\_jq\_jn\_j\\$ 
under Hardy-Weinberg equilibrium and \\$n\_j\\$ is the sample
size of variant \\$j\\$. The right-hand-side is equal to \\${\bf Db}\\$ where \\${\bf b}\\$ is the \\$p\times 1\\$ vector of marginal 
effect estimates from GWAS. LD sampling variance and LD score are 
as described in the `.info` file.

Post launch of a summary data model the sampled values of key model parameters are printed
to standard out in real time. For example, during an SBayesR analysis the following snapshot
shows what is being printed

```
Iter Pi1    Pi2    Pi3    Pi4    NnzSnp SigmaSq ResVar GenVar hsq    Rounding TimeLeft  
10   0.9513 0.0162 0.0170 0.0155 795    0.0001  0.8905 0.0260 0.0284 0.0000   0:16:39  
20   0.9507 0.0112 0.0172 0.0208 761    0.0001  0.9421 0.0213 0.0221 0.0000   0:8:19   
30   0.9479 0.0041 0.0201 0.0279 761    0.0000  0.9747 0.0249 0.0249 0.0000   0:5:32   
40   0.9448 0.0046 0.0187 0.0319 865    0.0000  0.9117 0.0203 0.0218 0.0000   0:4:9   
50   0.9445 0.0054 0.0228 0.0272 893    0.0000  1.0019 0.0223 0.0218 0.0000   0:6:38   
60   0.9409 0.0019 0.0223 0.0349 930    0.0000  1.0481 0.0265 0.0247 0.0000   0:5:31 
```   

The parameters printed are iteration number, which can be altered with `--out-freq`, proportion
of variants in distribution one (zero component) to four (can be altered by changing the \\$\pi\\$
and \\$\gamma\\$ vector lengths), number of non-zero effects in the model at iteration \\$i\\$,
the variance of the distribution of genetic effects \\$\sigma\_{\beta}^2\\$, the residual variance
\\$\sigma\_{\epsilon}^2\\$, the genetic variance \\$\sigma\_{g}^2\\$ SNP-based heritability, rounding accumulation (a value
substantially greater than zero indicates insufficent machine precision and the user should consider terminating the
program), and an estimate of
the time left to complete the analysis. Upon completion of the analysis a summary of the 
posterior mean and standard deviation of the key models parameters is reported to standard 
out

```
Posterior statistics from MCMC samples:

              Mean            SD             
       Pi1    0.942564        0.006791       
       Pi2    0.056959        0.006728       
       Pi3    0.000305        0.000274       
       Pi4    0.000172        0.000094       
    NnzSnp    1726.798706     199.570724     
   SigmaSq    0.002898        0.000375       
    ResVar    0.900596        0.004205       
    GenVar    0.103104        0.002216       
       hsq    0.102722        0.002039  
```   

The output files are the same as those from the individual data models and include:

```
test.log: # Text file of running status, intermediate output and final results
test.snpRes: # Text file of posterior statistics of SNP effects
test.parRes: # Text file of posterior statistics of key model parameters
test.mcmcsamples.SnpEffects: # Binary file of MCMC samples for the SNP effects
test.mcmcsamples.Par: # Text file of MCMC samples for the key model parameters
```   

#### *Extra commands for running SBayes analyses*

```
--unscale-genotype # Run analyses assuming genotypes are not scaled to unit variance. Default is scaled genotypes
--exclude-mhc # Exclude variants in the Major Histocompatibility Complex from the analysis. Highly recommended 
              # when performing genome-wide analyses
--exclude-region # Path to file where each line had three entries: chromosome, starting base pair position of 
                 # region and end base pair of region
--impute-n # Let GCTB impute the sample size for each variant using an algoirithm if you cannot trust 
           # the sample size reported. Filters variants with a sample size further than 3 SD from the median of 
           # the sample size distribution. The method assumes each SNP has been calculated from the same set 
           # of individuals.  
```   

#### Summary


The method implemented in GCTB assumes certain ideal data constraints such as summary data computed from a single set of individuals at fully observed genotypes
as well as minimal imputation error and data processing errors such as allele coding and frequency mismatch. Summary data in the public 
domain often substantially deviate from these ideals and can contain residual population stratification, which is not accounted for in this model. Practical
solutions to these ideal data deviations include the use of data that are imputed and the restriction of analyses to variants that are known to be imputed with 
high accuracy. We found that the simple filtering of SNPs with sample sizes that deviate substantially 
from the mean across all variants from an analysis, when using summary statistics from the public domain substantially improved 
model convergence. However, even the reported sample size cannot be trusted in some summary data sets.

We explored LD pruning of variants to remove variants in very high LD (\\$R^2>0.99\\$) but found that this did not substantially improve model 
convergence or parameter estimates although this was not formally assessed. However, removal of high LD regions, such as the MHC region 
improved model convergence for real traits. High LD regions are expected to have the potential to
be extreme sources of model misspecification with the model expecting summary data in to be very similar for 
variants in high LD. Small deviations due to data error not expected in the model likelihood at these loci thus have high potential to lead to model
divergence (see Zhu and Stephens\citep{zhu2017bayesian} for further discussion).  We strongly recommend reading the
[current manuscript](https://www.biorxiv.org/content/10.1101/522961v3?rss=1) and that of [Zhu and Stephens (2017)](https://projecteuclid.org/euclid.aoas/1507168840) as a preliminary to 
running the summary data based methods in GCTB.








### SBayesRC Tutorial

SBayesRC is a method that incorporates functional genomic annotations with high-density SNPs (> 7 millions) for polygenic prediction. Our manuscript is available at [here](https://www.nature.com/articles/s41588-024-01704-y).

This method is based on a low-rank approximation which utilises the eigenvalues and eigenvectors of block-wise LD correlation matrices. It requires SNPs that are included in the LD reference to be also present in the GWAS samples. Thus, the first step is to impute the GWAS summary statistics for any missing SNPs that are only present in the LD reference. If there are more than 30% of missing SNPs, we recommend to recompute the eigen-decomposition based on a matched LD reference sample (see below).

#### Input files
`ldm` is a folder comprising the eigen-decomposition data for each block. We have provided the data from a random sample of 20K unrelated UKB individuals of European ancestry using either ~1 million HapMap3 SNPs or 7 million imputed SNPs, which are available at [Download](https://cnsgenomics.com/software/gctb/#Download).

`test.ma` is the file of GWAS summary statistics with the following format
```
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
```

`annot.txt` is the annotation file, with columns being SNP ID, Intercept (a column of one), and annotation values. You can include an arbitary number of annotations.
```
SNP     Intercept   Anno1   Anno2   Anno3 
rs1001  1           1       0       0.2
rs1002  1           0       1       -1.5
rs1003  1           1       1       0.7
```

#### Step 1: QC & imputation of summary statistics for SNPs in LD reference but not in GWAS data.

```
gctb --ldm-eigen ldm --gwas-summary test.ma --impute-summary --out test --thread 4
```
This command will generate a text file `test.imputed.ma` which contains the summary statistics for all SNPs including those from the original GWAS file and those with imputed summary data. Use this file as input for the next step. `--thread 4` enables multi-thread (e.g., 4 threads) computing and if it is ignored, a single thread will be used. Before imputation, a QC step is carried out to match alleles between GWAS and reference samples and remove SNPs with per-SNP sample size beyond 3 standard deviation acround the median in GWAS.

**Tips for parallel computing:** The imputation step is fully available for parallel computing across blocks and the runtime is expected to decrease linearly with the number of threads, so you can use as many cores as possible in your own computing platform. Alternatively, for the best efficiency, you can impute each block separately by using `--block $i` where i starts from 1 to the total number of blocks:
```
gctb --ldm-eigen ldm --gwas-summary test.ma --impute-summary --block $i --out test
```
This command will generate a text file `test.block$i.imputed.ma` for each block. After all blocks are finished, use the following command to merge all `.ma` files across block into a single file with genome-wide SNPs:
```
gctb --gwas-summary test --merge-block-gwas-summary --out test.imputed
```

#### Step 2: Run SBayesRC with annotation data.

```
gctb --ldm-eigen ldm --gwas-summary test.imputed.ma --sbayes RC --annot annot.txt --out test --thread 4
```
This command will generate a text file for SNP effect estimates `test.snpRes`, text files for model parameters `test.parRes` and `test.parSetRes`, and a folder that stores the MCMC samples for all model parameters `test.mcmcsamples`.

#### Predicting polygenic scores.
`PLink` can be used to produce polygenic score for the target cohort, using columns of `test.snpRes` as input. For example,
```
plink --bfile target --score test.snpRes 2 5 8 header sum center --out target
```

#### Recomputing eigen-decomposition of block-wise LD matrics

In the case where eigen-decomposition needs to be performed, e.g., the GWAS data is genetically different from our LD reference, or GWAS SNPs are substantially different from our LD reference SNPs, you can use the following steps to generate the eigen-decompostion data based on your own LD reference sample:


#### *Step 1: compute block-wise LD matrices*
Skip this step if you are using the eigen-decomposition data provided by us. This step requires individual-level genotype data from LD reference `ldRef` and a text input file that defines the boundaries of each block. Here, we use text file [ref4cM_v37.pos](download/ref4cM_v37.pos) which partitions the human genome into 591 approximately independent LD blocks, each with a width of at least 4 cM, based on the results of [Berisa and Pickrell](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/).     

```
gctb --bfile ldRef --make-block-ldm --block-info ref4cM_v37.pos --out ldm --thread 4
```

This command will generate a folder called `ldm` as specified by `--out`. This folder will contain a serise of binary files for blocks `block*.ldm.bin`, a text file for genome-wide SNP information `snp.info`, and a text file for block information `ldm.info`. You can use `--write-ldm-txt` to print the content in the binary file into text file for each block, which will generate another file called `block*.ldm.txt`.

**Tips for parallel computing:** This step can be carried out by computing each chromosome or block separately. To compute chromosomes in parallel using array job submittion (`i` is the job index), you can do

```
gctb --bfile ldRef.chr$i --make-block-ldm --block-info ref4cM_v37.pos --out ldm --thread 4
```
where `ldRef.chr$i` is the bed file that only contains data for the ith chromosome (while the ldBlockInfo file `ref4cM_v37.pos` is still genome-wide). To compute all blocks in parallel, you can do

```
gctb --bfile ldRef --make-block-ldm --block $i --block-info ref4cM_v37.pos --out ldm
```
where `ldRef` needs to be the genome-wide bed file. Here using a single thread is sufficient.

After all chromosomes/blocks are done, run the command below to merge info files across blocks. Do not need to run this if you are NOT doing chromosomes/blocks separately.

```
gctb --ldm ldm --merge-block-ldm-info --out ldm
```

#### *Step 2: perform eigen-decomposition in each block*
Skip this step if you are using the eigen-decomposition data provided by us. The command below will perform eigen-decomposition for each block one by one and store the result in binary format `block*.eigen.bin` in folder `ldm`. Again, you can check the file content with `--write-ldm-txt` which will produce a text file.

```
gctb --ldm ldm --make-ldm-eigen --out ldm --thread 4
```

**Tips for parallel computing:** This step can also be done for each block in parallel by using array job submittion:

```
gctb --ldm ldm --make-ldm-eigen --block $i --out ldm
```

No merging is needed after completing all blocks.




### Genome-wide Fine-mapping analysis

The Genome-wide Bayesian Mixture Model (GBMM) implemented in GCTB (e.g., SBayesRC) can perform genome-wide fine-mapping analysis. These methods require summary-level data from genome-wide association studies (GWAS) and linkage disequilibrium (LD) data from a reference sample. Our manuscript is available at [here](https://www.medrxiv.org/content/10.1101/2024.07.18.24310667v3). 

We outline below (2 steps) on how to perform the genome-wide fine-mapping (GWFM) analysis and calculate the credible set using GCTB.

#### Step 1. Generate eigen-decomposition data for SNPs matched in the GWAS summary statistics 

You can use multi-threading to compute eigen-deomcomposition for all blocks:

```
gctb --ldm ukbEUR_13M_FullLDM --gwas-summary test.ma --make-ldm-eigen --thread 32 --out matched_ldm
```

**\--ldm** specifies a folder containing the block-wise full LD matrices. We have computed these matrices for 13 million SNPs using a sample of European ancestry. You can download this file (`ukbEUR_13M_FullLDM.zip`) at [here](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/). The block map file under genome build GRCh37 (`ref_b37_1588blocks.pos`) can be found in the same folder or [here](download/ref_b37_1588blocks.pos).

**\--gwas-summary** reads summary-level data from GWAS. The file format is as follows:
```
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
```

**\--out matched_ldm** specifies the output folder that stores the eigen-decomposition data for matched SNPs.

**\--make-ldm-eigen --thread 32** computes eigen-decomposition for all LD blocks using multi-threading (e.g., 32 threads).

Alternatively, when running on a HPC cluster, it's often more efficient to compute all blocks in parallel by using `--block [number]`:

```
gctb --ldm ukbEUR_13M_FullLDM --gwas-summary test.ma --make-ldm-eigen --block $i --out matched_ldm
```

**\--block $i** extract a specific block where `i` is the block index starting from 1.


#### Step 2. Run genome-wide fine-mapping analysis using the SBayesRC model

```
gctb --gwfm RC --ldm-eigen matched_ldm --gwas-summary test.ma --annot annot.txt --gene-map gene_map.txt --thread 32 --out test
```

**\--gwfm** specifies the method to perform genome-wide fine-mapping analysis, e.g., RC.

**\--ldm-eigen** reads the eigen-decomposition data generated from Step 1.

**\--annot** reads the annotation file. The file format is as follows, with columns being SNP ID, Intercept (a column of one), and annotation values. You can include an arbitary number of annotations. A set of 96 annotations for 13 million SNPs is available at [here](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/Annotation/). 
```
SNP     Intercept   Anno1   Anno2   Anno3 
rs1001  1           1       0       0.2
rs1002  1           0       1       -1.5
rs1003  1           1       1       0.7
```

**\--gene-map** specifies gene map file to annotate the nearest gene for the identified credible set, which is available at [here](download/gene_map_hg38_hg19.txt). This flag is optional. The file format is as follows. The genome build can be selected using **\--genome-build** (hg19 or hg38). hg19 is used as default. 
```
Ensgid          GeneName        GeneType        Chrom_hg38      Start_hg38      End_hg38        Chrom_hg19      Start_hg19      End_hg19
ENSG00000290825 DDX11L2         lncRNA          1               11869           14409           1               11869           14409
ENSG00000186092 OR4F5           protein_coding  1               65419           71585           1               65419           71585
```

**\--thread** specifies the number of threads to use. The recommended number is 32 (or any value which is a multiple of 4) because by default 4 MCMC chains are used which can be run in parallel (8 cores per chain if 32 is used).

**\--out** saves the genome-wide fine-mapping result files with the given prefix filename. The result files include (if test is fused as the filename): 

`test.snpRes` shows SNP effect and PIP estimates.
```
Index  Name         Chrom   Position   A1   A2   A1Frq      A1Effect    SE         VarExplained   PEP          PIP            GelmanRubin_R
1      rs12132974   1       801661     T    C    0.075000   -0.000026   0.000270   1.02e-08       0.010000     0.00662249     1.32981
2      rs12134490   1       801680     C    A    0.075000   0.000010    0.000138   2.67e-09       0.001250     0.00811219     1.44818
```
where columns are SNP position index, SNP name, SNP chromosome, SNP position, the eﬀect (coded) allele, the other allele, frequency of the eﬀect allele, estimated posterior eﬀect size, SE for the estimated posterior effect size, posterior SNP-based heritability enrichment probability (PEP), posterior inclusion probability (PIP), and Gelman & Rubin’s R statistics for convergence assessment.

`test.parRes` shows genetic architecture estimates. 

`test.parSetRes` shows annotation-specific parameter estimates.

`test.enrich` shows per-SNP heritability enrichment estimate in each annotation.

`test.mcmcsamples` shows MCMC samples of SNP effects in binary format, which is used to compute credible sets’ posterior heritability enrichment probability (PEP).

`test.skepticalSNPs` shows a list of SNPs whose joint effects are set to zero because their effects are likely blowing up due to poor data quality or large difference in LD between the GWAS and reference samples.

`test.lcs` shows the identified local credible sets.
```
CS   Size  PIP        PGV        PGVenrich     PEP         NumUnconvgSNPs  SNP                  ENSGID_hg19       GeneName_hg19 
1    1     1.000000   0.000379   437.109253    1.000000    0               rs2066827            ENSG00000111276   CDKN1B
2    2     1.025764   0.000162   93.409187     0.956250    0               rs2641670,rs2646108  ENSG00000103196   CRISPLD2
```
where columns are credible set (CS) index, number of SNPs in the CS, cumulative PIP across SNPs in the CS, proportion of genetic variance explained by the CS, fold enrichment of genetic variance explained, posterior enrichment probability, number of SNPs with Gelman & Rubin’s R statistic > 1.2 (indicating likely non-convergence), SNPs in the CS, ensemble ID for the nearest gene to the CS, name of the nearest gene to the CS.

`test.lcsRes` shows the summary statistics for the identified local credible sets. For example,
```
                                   PIP threshold:           0.9
                                   PEP threshold:           0.7
                   Number of 1-SNP credible sets:           313
               Number of multi-SNP credible sets:          1724
           Total number of SNPs in credible sets:          6389
                       Average credible set size:           3.1
       Estimated total number of causal variants:       42712.8
                                 Estimated power:        0.0477
  Estimated number of identified causal variants:        2039.3
      Estimated proportion of variance explained:        0.3713
```

`test.gcs` shows the identified global credible sets. The first column shows the estimated percentage of causal variants included by the global credible set.
```
Alpha          SNP          PIP
 0.01    rs4670775            1
 0.01    rs2066827            1
```

`test.gcsRes` shows the summary statistics for the identified global credible sets.
```
Threshold     CS_size     Prop_hsq
    0.01          457     0.129471
    0.05         3945     0.323758
     0.1        12160     0.458595
     0.2        41186     0.614211
     0.3        86449     0.708990
     0.4       148123     0.776145
     0.5       227297     0.830163
     0.6       326394     0.875781
     0.7       450175     0.914510
     0.8       607494     0.947969
     0.9       817326     0.976685
       1      1153584     1.000000
```

#### (Re)calculate credible sets

```
gctb --cs --pwld-file ldm/rsq0.5.pwld --pip 0.9 --pep 0.7 --gene-map gene_map.txt --flank 5000 --genome-build hg19 --mcmc-samples test --out test 
```
By default, PIP threshold of 0.9 and PEP threshold of 0.7 are used to construct local credible sets (LCSs). The following command can be used to recompute LCSs with different threshold values.

**\--cs** calculate local and global credible sets.

**\--pip** specifies the threshold for the coverage of the credible set. 

**\--pep** specifies the threshold for PEP.

**\--flank** specifies flanking window to define the nearest gene, with 5000 being the default value. 

**\--genome-build** specifies genome build of the gene map file, e.g. hg19 or hg38. By default, hg19 is used.

**\--mcmc-samples** reads the MCMC samples output from the genome-wide fine-mapping analysis.  

**\--out** saves the full local credible set results in .lcs and summary of local credible set result in .lcsRes. Saves the full global credible set results in .gcs and summary of global credible set result in .gcsRes.

**\--pwld-file** reads the LD file containing pairwise LD \\$r^2\\$ > 0.5 between SNPs.  

The required LD file with pairwise LD \\$r^2\\$ > 0.5 in each LD block is included in the `ldm` folder. This file can be obtained using GCTB and eigen-decomposition LD reference with command line below:

```
gctb --get-pwld --ldm-eigen ldm --rsq 0.5 --thread 8 --out rsq0.5
```

**\--get-pwld** compute pairwise LD file using eigen-decomposition based LD reference files.

**\--ldm-eigen** specifies a folder containing the eigen-decomposition LD reference data for each block. 

**\--rsq** specifies the \\$r^2\\$ threshold used to output the pairwise LD result. 

**\--out** saves the LD file with pairwise LD \\$r^2\\$ > 0.5.


#### Prediction of fine-mapping power 
This [shyny app](https://gctbhub.cloud.edu.au/shiny/power/) predicts the genome-wide fine-mapping power along with the sample size based on the genetic architecture parameters.
