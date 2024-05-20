# BV-LDER-GE
We propose a novel statistical method to study the genome-level Gene-Environment interaction effects using human complex trait summary statistics: BiVariate Linkage-Disequilibrium Eigenvalue Regression for Gene-Environment interactions (BV-LDER-GE). The major contribution of BV-LDER-GE comes from two aspects: joint modeling and testing of GxE interaction variance and the genetic covariance between additive genetic effect and GxE interaction effect (GxE genetic covariance); more accurate estimation of the GxE genetic covariance utilizing full LD information. 

:open_book: Citation:

To be added.

Acknowledgement: BV-LDER-GE package is changed based on the original LDER package. If you are studying narrow-sense heritability, please use and refer to https://github.com/shuangsong0110/LDER). We thank Shuang Song for sharing the original LDER code.

## Table of contents
* [Install and LD preparation](#hammer-install-and-ld-preparation)
* [Run BV-LDER-GE main function](#rocket-run-bv-lder-ge-main-function)
* [Output](#bulb-output)
* [A Simplified Pipeline](#key-a-simplified-pipeline)
* [Run BV-LDER-GE-adj function to adjust for the genetic effect of E](#key-run-bv-lder-ge-adj-function-to-adjust-for-the-genetic-effect-of-e)

  
## :hammer: Install and LD preparation
BV-LDER-GE R package requires R >= 3.5.0 and Python 3.6.
BV-LDER-GE can be installed using the command:
```r
devtools::install_github('dongzhblake/BiVariate-LDER-GE')
```

 For precise estimation, we recommend using in-sample LD or LD sample with higher number of sample sizes such as UKBB. Below we provide two pre-computed LD information using UKBB European ancestry. We include functions to calculate your own in-sample LD. For a detailed LD reference panel preparation please refer to https://github.com/dongzhblake/LDER-GE?tab=readme-ov-file#hammer-install-and-ld-preparation.

### Use the pre-computed LD information

The pre-computed LD information for 396,330 hapmap3 variants from 276,050 UK Biobank European individuals can be manually downloaded from [https://drive.google.com/file/d/1mvDA79qPAoPXUjmUC1BQw-tInklZ4gPD/view?usp=drive_link](https://drive.google.com/file/d/1CCGil-ZnXourrk5JFJeqyEKIzLCvX2vh/view?usp=drive_link)

The pre-computed LD information for 966,766 hapmap3 and array variants from 307,259 UK Biobank European individuals can be manually downloaded from [https://drive.google.com/file/d/1UF1xP1Rg1JiFMozkFJ3bbJmFggDY8-5Y/view?usp=drive_link](https://drive.google.com/file/d/1FgikyxYd_jW05aLuHgj0yChW9i2ZgtT9/view?usp=sharing)

After downloading, decompress the files:

`unzip UKB396kvariant_hm3.zip`

`unzip UKB966kvariant_hm3.zip`

## :rocket: Run BV-LDER-GE main function
The main funcion can be run with:

```r
runBV_LDER_GE(assoc_gwis,
assoc_gwas,
n.gwis,
n.gwas,
path,
n.ld,
method)

```
- assoc_gwis GWIS (GE interaction effect) summary statistics, need to include snp, chr, a0, a1, z (header is necessary)

- assoc_gwas GWAS (additive genetic effect) summary statistics, need to include snp, chr, a0, a1, z (header is necessary)
  
- n.gwis The sample size of the GWIS (GE interaction effect) summary statistics
  
- n.gwas The sample size of the GWAS (additive genetic effect) summary statistics

- path The path of LD panel directory

- n.ld The sample size of the LD reference panel

- method (optional): Default='lder'. We also provide a choice of 'both', which outputs the results for both LDER and LDSC.


## :bulb: Output

If `method='lder'`, the `runBV_LDER_GE` function returns a list with 10 elements:

`genecov`: Estimated GxE genetic covariance

`h2I`: Estimated GxE proportion

`h2g`: Estimated narrow-sense heritability

`genecov.se`: The standard error of estimated GxE genetic covariance with block-jackknife.

`h2I.se`: The standard error of estimated GxE proportion with block-jackknife.

`h2g.se`: The standard error of estimated narrow-sense heritability with block-jackknife.

`genecov.p`: The P value for testing the estimated GxE genetic covariance.

`h2I.p`: The P value for testing the estimated GxE proportion.

`h2g.p`: The P value for testing the estimated narrow-sense heritability.

`lder.BV_test`: The P value for the joint modeling test of the GxE genetic covariance and GxE proportion (Still test test GxE proportion).

If `method='both'`, the `runBV_LDER_GE` function returns a list containing the results of both LDER-GE and LDSC-based methods.


## :key: A Simplified Pipeline
Download a sample GWIS (GxE interaction sumstats) summary statistics at 

https://drive.google.com/file/d/1GoULJNPVHVsOR2dv_4JpaC5eaGOeVxMA/view?usp=drive_link

Download a sample GWAS (additive genetic effect sumstats) summary statistics at 

https://drive.google.com/file/d/1ZYefklTqBvCxkoBkrEl3Q-rimGjhloKu/view?usp=drive_link


Download the pre-computed LD information for 396,330 hapmap3 variants from 276,050 UK Biobank European individuals at

https://drive.google.com/file/d/1mvDA79qPAoPXUjmUC1BQw-tInklZ4gPD/view?usp=drive_link


`unzip UKB396kvariant_hm3.zip`


Run with R:

```r
devtools::install_github('dongzhblake/BiVariate-LDER-GE')
library(BVLDERGE)
library(data.table)
path0 <- "UKB396kvariant_hm3" # the complete system path to this LD folder
assoc <- fread('LDER_GE_exampleGWIS.txt')
assoc_gwis <- fread('BV_LDER_GE_exampleGWIS.txt')
assoc_gwas <- fread('BV_LDER_GE_exampleGWAS.txt')
n.gwis = median(assoc_gwis$n)
n.gwas = median(assoc_gwas$n)

res=runBV_LDER_GE(assoc_gwis, assoc_gwas, n.gwis, n.gwas, path = path0,  n.ld=276050 ,method='both')

> unlist(res)
lder.genecov           lder.h2I           lder.h2g    lder.genecov.se 
5.569623e-02       1.071866e-01       8.660379e-02       6.902761e-03 
lder.h2I.se        lder.h2g.se lder.cov_rrhog_h2I     lder.genecov.p 
9.046920e-03       8.844408e-03       2.741344e-05       7.105664e-16 
lder.h2I.p         lder.h2g.p       ldsc.genecov           ldsc.h2I 
2.207742e-32       1.219501e-22       5.778562e-02       1.076898e-01 
ldsc.h2g    ldsc.genecov.se        ldsc.h2I.se        ldsc.h2g.se 
8.037813e-02       8.228211e-03       1.156984e-02       1.035050e-02 
ldsc.cov_rrhog_h2I     ldsc.genecov.p         ldsc.h2I.p         ldsc.h2g.p 
4.356295e-05       2.173623e-12       1.305036e-20       8.124193e-15 
lder.BV_test       ldsc.BV_test
2.025652e-33       1.229271e-21 

```

## :key: Run BV-LDER-GE-adj function to adjust for the genetic effect of E

Download a sample GWIS (GxE interaction sumstats) summary statistics of Y at 

https://drive.google.com/file/d/1SvQZqHu-K8QNLaeuVPydsoNAOmIRjn7S/view?usp=drive_link

Download a sample GWAS (additive genetic effect sumstats) summary statistics  of Y at 

https://drive.google.com/file/d/15U0P1GvvG42h3f4VLqJZd8tglhW2gaxC/view?usp=drive_link

Download a sample GWAS (additive genetic effect sumstats) summary statistics  of E at 

https://drive.google.com/file/d/1YjzrCTezLLA0VyxNATpBMG9If_2st6Iw/view?usp=drive_link

Download the pre-computed LD information for 396,330 hapmap3 variants from 276,050 UK Biobank European individuals at

https://drive.google.com/file/d/1mvDA79qPAoPXUjmUC1BQw-tInklZ4gPD/view?usp=drive_link


`unzip UKB396kvariant_hm3.zip`


Run with R:

```r
devtools::install_github('dongzhblake/BiVariate-LDER-GE')
library(BVLDERGE)
library(data.table)
path0 <- "UKB396kvariant_hm3" # the complete system path to this LD folder
assoc_gwis_Y <- fread('BV_LDER_GE_example_adj_GWIS_Y.txt')
assoc_gwas_Y <- fread('BV_LDER_GE_example_adj_GWAS_Y.txt')
assoc_gwas_E <- fread('BV_LDER_GE_example_adj_GWAS_E.txt')
n.gwis_Y = median(assoc_gwis_Y$n)
n.gwas_Y = median(assoc_gwas_Y$n)
n.gwas_E = median(assoc_gwas_E$n)

res=runBV_LDER_GE_adj(assoc_gwis_Y, assoc_gwas_Y,assoc_gwas_E, n.gwis_Y, n.gwas_Y, n.gwas_E, path0,  n.ld=276050 ,method='both')


> unlist(res)
lder.genecov           lder.h2I    lder.genecov.se        lder.h2I.se 
-8.669521e-03       8.185539e-02       8.533855e-03       8.825768e-03 
lder.cov_rrhog_h2I     lder.genecov.p         lder.h2I.p       lder.gcov_IE 
-1.048700e-05       3.096782e-01       1.783005e-20       1.066853e-01 
lder.gcov_IE.se     lder.gcov_IE.p       ldsc.genecov           ldsc.h2I 
8.504955e-03       4.293074e-36      -7.798517e-03       9.876998e-02 
ldsc.genecov.se        ldsc.h2I.se ldsc.cov_rrhog_h2I     ldsc.genecov.p 
1.090945e-02       1.112378e-02      -2.508875e-05       4.747077e-01 
ldsc.h2I.p       ldsc.gcov_IE    ldsc.gcov_IE.se     ldsc.gcov_IE.p 
6.736042e-19       1.048404e-01       1.057090e-02       3.482327e-23 
lder.BV_test       ldsc.BV_test 
2.016574e-19       3.937359e-18

The output is largely the same with BV-LDER-GE main function, but the output of genecov and h2I are adjusted for the genetic effect of E.
There are three new output about gcov_IE term: the GxE genetic covariance between the addtive genetic effect of E and GxE interaction effect of Y. If gcov_IE is tested positive, then we suggest make adjustment. If gcov_IE is not tested positive, we suggest the unadjusted result from BV-LDER-GE for better estimation quality.

## :busts_in_silhouette: Maintainer

Please contact Zihan Dong (zihan.dong@yale.edu) if there are any problems or questions.


