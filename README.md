## Random-Open-parents Association Mapping (ROAM)
#### These R scripts are designed for ROAM, a new proposed multi-parental population for large scale genetic analysis, including genotype imputation, projection, bin extraction, kinship calculation, joint linkage mapping (JLM) and genome-wide association study (GWAS) analyses. The related paper is submitted.
>## Related papers:
> #### Xiao et al. Genome-wide dissection of the maize ear genetic architecture using multiple populations. *New Phytol*. 2016; 210(3):1095-106. doi: 10.1111/nph.13814.

==========================================

##Manual

This R scripts are designed for genotype data imputation and projection, and also bin extraction, if necessary, joint linkage mapping (JLM) and genome-wide association study (GWAS).

It comprises six parts, 1. imputation, 2. projection, 3. exbin, 4. kinship, 5. JLM, 6. GWAS.
All the original input files named as "*_in ", and output result files named as "*_out ", expect for the GWAS where the required input files had already saved an RData file.

1. imputation
This program is used for filling the missing (coding “0”) and heterozygous (coding “3”) SNPs in linkage map. The filling is based on the flank SNPs which is not missing or heterozygous and also the physical position. The sample data is chromosome 1 of the first population which consists of 200 lines and 200 SNPs. You should run it in each chromosome.

2. projection
This step projects the genotype data from 200 SNPs which are contained in the output file of last step, to 300 SNPs as the simulated sample data presented. This program is similar to the imputation procedure.

3. exbin
If you want to reduce the SNP number as existing redundant SNPs will highly increase the computational load, this program will help to extract bin from SNP information. The sample data is chromosome 10 of two populations that each consists of 100 lines. The column one in the result file is the bin number, as the sample data presents that it divided into 194 bins, and column three is the original physical position, so that you can pick up the first SNP in each bin to the following analysis.

4. kinship
The kinship matrix is essential for the JLM procedure. This step is used for calculating the marker-based kinship matrix. The sample data comprises ten populations that each has 200 lines and supposing each population has two parents coding “1,2”, “3,4”, “5,6” etc., and 3000 SNPs in genome. There are two output files, one is 2000*2000 kinship matrix, and another RData file is eigenvalue decomposition result which consists of eigenvalue in the first column and the remaining are the eigenvectors. 

5. JLM
The joint linkage mapping analysis is based on a linear mixed model which contains marker effect and polygenic effect as the random effect. The input files include genotype data which same to the kinship step, phenotype data which normalized within each population, kinship matrix and eigenvalue/eigenvectors which gave in the prior step. The procedure will first solve the null model which exclude the marker term so that the output lambda file has three values, variance of polygenic effect, residual error and the ratio of the two which named lambda here, respectively. The parameter file consists of the likelihood ratio test (LRT) score which used to test the significant of each marker and the blup file consists of the effect size and corresponding standard errors.

6. GWAS
Genome-wide association study is based on a multivariate regression model. To control the polygenic background effect, the GWAS was performed per chromosome at a time. For each chromosome, we forced population effects and the effects of QTLs detected by SLM and JLM from other chromosomes to be included in a general linear model. The residual of the above GLM was then treated as the dependent variable to scan all SNPs on the current chromosome. The example data “10pop.example.RData” include four input files: gen, phe, inf and qtl, which represent the numeric genotype of joint populations, phenotype value of target trait, marker information and QTL information regarding as covariates, respectively. The program output two example files, “trait_r1.bpp.res.xlsx” and “trait_r1.Genome.bin.info”. More detailed information about the GWAS procedure can be inferred from the script comments.

For more information, please contact:
Hao Tong (tonghao0605@gmail.com) or Yingjie Xiao (shanren0179@gmail.com).
