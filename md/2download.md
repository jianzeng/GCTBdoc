
## Download {: .expand}
### Executable
[gctb\_2.5.4\_Linux.zip](download/gctb_2.5.4_Linux.zip) (*Lastest version updated in 21 April 2025*)

### Source code
[GCTB 2.5.4](download/gctb_2.5.4_scr.zip)

### Tutorial data
[GCTB tutorial data](download/gctb_2.0_tutorial.zip)

### Eigen-decomposition data of LD matrices
The eigen-decomposition data are for SBayesRC and SBayesR with the low-rank model. In the follwoing link, we provide data derived from unrelated UKB individuals of Europan (EUR), East Asian (EAS), and African (AFR) ancestires. See our [manuscript](https://www.nature.com/articles/s41588-024-01704-y) and [Tutorial](https://cnsgenomics.com/software/gctb/#Tutorial) for details.

* [1M HapMap3 SNPs](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/)
* [7M Imputed SNPs](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/Imputed/)

### Functional genomic annotations
Download the formatted data for per-SNP functional annotations derived from [S-LDSC BaselineLDv2.2](https://www.nature.com/articles/ng.3954). 
* [7M SNP annotations](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/)
* [13M SNP annotations](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/Annotation/)

Gene region annotations for genome build hg19 and hg38 are avaialble [here](download/gene_map_hg38_hg19.txt).

### LD matrices
Blockwise full LD matrices for 13 million SNPs (MAF > 0.1%) in a random sample of 10K unrelated individuals of European ancestry in the UK Biobank.
* [Blockwise full LD matrix (13 million SNPs)](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/)

The following LD matrices were computed based on 1.1 million common SNPs in a random sample of 50K unrelated individuals of European ancestry in UK Biobank dataset unless otherwise noted.

* [Shrunk sparse matrix](https://zenodo.org/record/3350914#.XyFfnC17G8o)
* [Shrunk sparse LD matrix (2.8 million common SNPs)](https://zenodo.org/record/3375373#.XyFgOS17G8o)

In the shrunk sparse matrices, described in [Lloyd-Jones et al. (2019)](https://www.nature.com/articles/s41467-019-12653-0), the observed LD correlations computed from a reference sample were shrunk toward the expected values defined by a [genetic map](https://github.com/joepickrell/1000-genomes-genetic-maps), following the algorithm in [Wen and Stephens (2010)](https://projecteuclid.org/euclid.aoas/1287409368). After shrinkage, LD correlations smaller than a threshold (default 1e-5) were set to be zero to give a sparse format, which is more efficient in storage and computation. 

* [Sparse matrix (including MHC regions)](https://cnsgenomics.com/data/GCTB/ukbEURu_imp_v3_HM3_n50k.chisq10.zip)

The sparse matrices described in [Zeng et al. (2021)](https://www.nature.com/articles/s41467-021-21446-3) were computed by setting the likely chance LD to zero based on a chi-squared test (default threshold at chi-squared test statistic of 10).

* [Banded matrix (including MHC regions)](https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip)

While the shrunk sparse matrices were used in our original SBayesR paper, [Prive et al. (2021)](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa1029/6039173) found that using a banded matrix with a window size of 3 cM per SNP can improve prediction accuracy. Therefore, we have created such a LD matrix in GCTB format for SBayesR analysis.

### Summary data and PGS weights
The summary data for the 50 (including 28 approximately independent) UKB traits analysed in [Zheng et al. 2024](https://www.nature.com/articles/s41588-024-01704-y) can be downloaded here: [summary data](https://gctbhub.cloud.edu.au/data/SBayesRC/share/v1.0/summary/), and the corresponding PGS weights can be found here: [PGS weights](https://gctbhub.cloud.edu.au/data/SBayesRC/share/v1.0/PGS/). The PGS weights are joint effect estimates derived from ~7 million genome-wide SNPs, therefore it’s important to have matched SNP set between training and validation datasets. This is because if some important SNPs present in the training are missing in the validation, the genetic effects captured by these SNPs will be lost. To maximise the utility of the joint SNP weights, we recommend considering genotype imputation or rerunning SBayesRC with the matched set of SNPs. In addition, note that the PGS weights were estimated using samples of the European ancestry, so they will perform best when applying to individuals of the European ancestry.

### Genome-wide fine-mapping (GWFM) analysis results
Results of SNP PIPs, local credible sets, and globle credible sets for 599 traits can be downloaded here: [GWFM results](https://gctbhub.cloud.edu.au/data/SBayesRC/share/Finemap/v1.1/)

### Older versions

gctb_2.5.2: [[gctb\_2.5.2\_Linux.zip](download/gctb_2.5.2_Linux.zip)] [[source code](download/gctb_2.5.2_scr.zip)]

gctb_2.5.1: [[Linux executable](download/gctb_2.5.1_Linux.zip)] [[source code](download/gctb_2.5.1_scr.zip)]

gctb_2.05beta: [[Linux executable](download/gctb_2.05beta_Linux.zip)] [[source code](download/gctb_2.05beta_scr.zip)]

gctb_2.04.3: [[Linux executable](download/gctb_2.04.3_Linux.zip)] [[source code](download/gctb_2.04.3_scr.zip)]

gctb_2.03beta: [[Linux executable](download/gctb_2.03beta_Linux.zip)]  [[source code](download/gctb_2.03beta_scr.zip)]

gctb_2.02: [[Linux executable](download/gctb_2.02_Linux.zip)]

gctb_2.0: [[Linux executable](download/gctb_2.0_Linux.zip)]  [[source code](download/gctb_2.0_scr.zip)] [[MPI version source code](download/gctb_2.0_mpi_scr.zip)]

gctb_1.0: [[Linux executable](download/gctb_1.0_Linux.zip)] [[Mac executable](download/gctb_1.0_Mac.zip)] [[source code](download/gctb_1.0_scr.zip)] [[MPI version source code](download/gctb_1.0_mpi_scr.zip)]


The MPI version implements a distributed computing strategy that scales the analysis to very large sample sizes. A significant improvement in computing time is expected for a sample size > 10,000. The MPI version needs to be compiled on user’s machine. See [README.html](download/README.html) in the tarball for instructions of compilation and usage. A testing dataset is also included in each tarball.


### Update log 

**1.**  1 Dec, 2017: first release.

**2.**  8 Feb, 2019: version 2.0 includes summary-data-based Bayesian methods.

**3.** 31 Jul, 2020: version 2.02 reports a message for convergence issue.

**4.** 15 Oct, 2021: version 2.03beta includes a robust parameterisation that attempts to address convergence issue.

**5.** 12 Dec, 2022: version 2.04.3 includes additional random effect component in BayesC and BayesR models to capture non-SNP random effects.

**6.** 11 April, 2023: version 2.05beta implements SBayesRC (SBayesR with the low-rank model) for polygenic prediction incorporating functional genomic annotations.

**7.** 29 April, 2024: version 2.5.1 adds functions for fine-mapping and to remove problematic SNPs during MCMC in SBayesRC.

**8.** 20 June, 2024: version 2.5.2 fixed bugs regarding -nan SigmaSq results and adds LD-based approach for computing credible sets.

**9.** 21 April, 2025: version 2.5.4 add --gwfm flag to perform genome-wide fine-mapping analysis.
