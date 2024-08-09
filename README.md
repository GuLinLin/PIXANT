# PIXANT [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/GuLinLin/PIXANT/issues) [![](https://img.shields.io/badge/Release-v0.1.0-important.svg)](https://github.com/GuLinLin/PIXANT/archive/refs/heads/master.zip) [![](https://img.shields.io/badge/license-GPL3.0-blue.svg)](https://github.com/GuLinLin/PIXANT/blob/master/LICENSE)<br>

## *multi-[P](https://github.com/GuLinLin/PIXANT/)henotype [I](https://github.com/GuLinLin/PIXANT)mputation based on mi[X](https://github.com/GuLinLin/PIXANT/)ed f[A](https://github.com/GuLinLin/PIXANT/)st ra[N](https://github.com/GuLinLin/PIXANT/)dom fores[T](https://github.com/GuLinLin/PIXANT/)*<br>

![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")

## Brief introduction <br>
Deep phenotype datasets can enhance the power of genetic analysis such as
genome-wide association study (GWAS), but recurrence of missing phenotypes
compromises the potentials of such resources. Here we address the central issue of
missing phenotypes in studies with any level of relatedness between phenotypes. we
propose a multi-phenotype imputation method that is scalable to large data over a
million individuals. We call our method PIXANT. PIXANT models nonlinear and
linear effects across multi-phenotype correlation and higher-order interactions
between predictive factors, and brings out unbiased imputation of much higher
accuracy than state-of-the-art methods. Tested in a dataset of n=20,000 individuals
and p=30 phenotypes, PIXANT is ~24.45 times faster and uses only about one ten
thousandth memory compared with PHENIX( ***`Nature Genetics 2016 (48):466–472`*** ).
Moreover, real data set analysis and biologically plausible results suggest that our
method imputation can uncover new true positive results. <br>

***`PIXANT`*** is developed by [***Linlin Gu***](https://github.com/GuLinLin)

## The flow chart of PIXANT <br>
![image](https://github.com/GuLinLin/PIXANT/blob/master/PIXANT.bmp)

## Version and download <br>
* [Version 0.1.0](https://github.com/GuLinLin/PIXANT/blob/master/PIXANT_0.1.0.tar.gz) -First version released on August, 21th, 2023<br>

## The function details of PIXANT <br>
```
Table of Contents
│
├ PIXANT                            # The function imputed a large data with missing values using the PIXANT.
│ ├ Usage
│ │ ├                               # PIXANT(data, aimPhenName, maxIterations = 20, maxIterations0 = 20, num.trees = 100,
│ │ ├                                        initialLinearEffects = 0, errorTolerance = 0.001, aimPhenMissingSize = 500,
│ │ ├                                        initialImputeType = 'random', refPhenThreshold = 0.3, maxNum.refPhen = 10,
│ │ └                                        SC.Threshold = 0.6, seed = 123, decreasing = TRUE, verbose = TRUE, xtrue = NA)
│ ├ Arguments
│ │ ├ data                          # A data frame with missing values.
│ │ ├ aimPhenName                   # The column name for the imputed phenotype.
│ │ ├ maxIterations                 # Stop after how many iterations (default = 20).
│ │ ├ maxIterations0                # The maximum iteration times of mixed fast random forest (default = 20).
│ │ ├ num.trees                     # How many trees are grown in the mixed fast random forest (default = 100).
│ │ ├ initialLinearEffects          # The initial values for linear effects (default = 0).
│ │ ├ errorTolerance                # The tolerance for log-likelihood (default = 0.001).
│ │ ├ aimPhenMissingSize            # The missing values size of settings in imputed phenotype (default = 500).
│ │ ├ initialImputeType             # Initial imputation method for missing values in a data frame. Currently support: "random", "average", "median".
│ │ ├ refPhenThreshold              # Relevance threshold for the selection of reference phenotypes, the range from 0 to 1 (default = 0.3).
│ │ ├ maxNum.refPhen                # Maximum number of reference phenotypes to be selected (default = 10).
│ │ ├ SC.Threshold                  # Individual's imputation quality threshold (default = 0.6).
│ │ ├ seed                          # Random seed. Default is 123, which generates the seed from R. Set to 0 to ignore the R seed.
│ │ ├ decreasing                    # (boolean) If TRUE the columns are sorted with decreasing amount of missing values.
│ │ ├ verbose                       # (boolean) If TRUE then missForest returns error estimates, runtime and if available true error during iterations.
│ │ └ xtrue                         # The complete data(a vector, matrix or data frame).
│ ├ Value                           # Return a list, the list contains:
│ │ ├ ximp                          # The imputed data (a data frame).
│ │ ├ imputePhen.accuracy           # Imputation accuracy of imputed phenotype.
│ │ ├ imputePhen.value              # The missing data of settings in the imputed phenotype, including sample index, observed values and imputed values (a data frame).
│ │ ├ imputePhen.refPhen            # A list of reference phenotypes for impute phenotype.
│ │ ├ imputePhen.missingRate        # Original missing rate of the imputed phenotypes.
│ │ ├ imputePhen                    # The impute phenotype, including imputed values and SC for each samples (a data frame).
│ │ ├ imputePhen.filter.missingRate # Missing rate after SC quality control of the imputed phenotype.
│ │ ├ imputePhen.r2                 # The r square of fitting between observed values and imputed values in the imputed phenotype.
│ └ └ imputePhen.pValue             # The p value of fitting between observed values and imputed values in the imputed phenotype.
│
├ PIXANT.eval                       # The function is the correlation between these imputed phenotypes and their true hidden values.
│ ├ Usage
│ │ └                               # PIXANT.eval(ximp, xmis, xtrue)
│ ├ Arguments
│ │ ├ ximp                          # The imputed data(a vector, matrix or data frame).
│ │ ├ xmis                          # A vector, matrix or data frame with missing values.
│ │ └ xtrue                         # The complete true data(a vector, matrix or data frame).
│ └ Value                           # Return the correlation coefficient between the real values and the imputed values.
│
├ prodNA                            # The function to produce missing values in a given and data set completely at random.
│ ├ Usage
│ │ └                               # prodNA(x, noNA, seed)
│ ├ Arguments
│ │ ├ x                             # A vector, matrix or data frame.
│ │ ├ noNA                          # Proportion of missing values to add to x. In case x is a data frame, noNA can also be a vector of probabilities per column or a named vector (see examples).
│ │ └ seed                          # An integer seed.
│ └ Value                           # Return a vector, matrix or data frame with missing values.
│
├ simG                              # Simulated Genome Relationship Matrix.
│ ├ Usage
│ │ └                               # simG( N, k, fam_size) 
│ ├ Arguments
│ │ ├ N                             # The number of individuals and must be a positive integer.
│ │ ├ k                             # Coefficient of kinship and the value ranges from 0 to 1.
│ │ └ fam_size                      # The size of the family, fam_size must be a positive integer and must divide N.
│ └ Value                           # Return Genome Relationship Matrix.
│
├ simPhen                           # Simulated phenotypic value.
│ ├ Usage
│ │ └                               # simPhen(N=N, P=P, K=G, h2=rep(0.6, P), B, E)
│ ├ Arguments
│ │ ├ N                             # The number of individuals.
│ │ ├ P                             # The number of phenotypes.
│ │ ├ K                             # A genome relational matrix.
│ │ ├ h2                            # The heritability of each phenotype in individuals.
│ │ ├ B                             # Genetic covariance. (allow the missing).
│ │ └ E                             # Environmental or residual covariance.(allow the missing)
│ └ Value                           # Return Simulated phenotype.
│
├ imputeUnivariate                  # Fills missing values of a vector, matrix or data frame by sampling with replacement from the non-missing values. For data frames, this sampling is done within column.
│ ├ Usage
│ │ └                               # imputeUnivariate(x, initialImputeType, v = NULL, seed = NULL)
│ ├ Arguments
│ │ ├ x                             # A data frame with missing values.
│ │ ├ initialImputeType             # Initial imputation method for missing values in a data frame. Currently support: "random", "average", "median".
│ │ ├ v                             # A character vector of column names to impute (only relevant if x is a data frame). The default NULL imputes all columns.
│ │ └ seed                          # An integer seed.
│ └ Value                           # Return complete data(a data frame).
│
├ PhenAdj                           # Adjusted phenotypic values base on covariates.
│ ├ Usage
│ │ └                               # PhenAdj(Phen, Cov)
│ ├ Arguments
│ │ ├ Phen                          # Phenotype file. The missing values should be denoted by NA.
│ │ └ Cov                           # A matrix of covariates. Each row is a sample and each column corresponds to one covariate. For example, age, gender.
│ └ Value                           # Return adjusted phenotype file.
│ 
├ SC                                # Estimating the SC (Sample Score) for each phenotype of each individual.
│ ├ Usage
│ │ └                               # SC(Phen, aimPhenName, use="pairwise",method="pearson",adjust="fdr",alpha=.05)
│ ├ Arguments
│ │ ├ Phen                          # Phenotype file. The missing values should be denoted by NA.
│ │ ├ aimPhenName                   # The column name for the imputed phenotype.
│ │ ├ use                           # use="pairwise" is the default value and will do pairwise deletion of cases. use="complete" will select just complete cases.
│ │ ├ method                        # method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall". These last two are much slower, particularly for big data sets.
│ │ ├ adjust                        # What adjustment for multiple tests should be used? ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). See p.adjust for details about why to use "holm" rather than "bonferroni").
│ │ └ alpha                         # alpha level of confidence intervals.
└ └ Value                           # Return the SC (Sample Score) for imputed phenotype of each individual.
```

## GETTING STARTED
### Installation
Please install it locally as following:
```sh
$ wget https://github.com/GuLinLin/PIXANT/blob/main/PIXANT_0.1.0.tar.gz
$ R CMD INSTALL PIXANT_0.1.0.tar.gz
```
After installed successfully, the ***`PIXANT`*** package can be loaded by typing
```r
library("PIXANT")
```
## INPUT
### Phenotype file
The very basic input data for PIXANT is a data matrix, in which each row represents an individual and each column represents a trait. The input data contain a header row, which may represents the Phenotypic names. The missing values should be denoted by NA.
Notice that only the numeric values are allowed and the characters will not be recognized. However, if a phenotype 
takes only values of 0, 1 (or only two levels), ***`PIXANT`*** would consider it to be a case-control phenotype, 
and the predicted value could be directly interpreted as the probability of being a case. <br>

> `phenotype.txt`

| Phenotype1 | Phenotype2 | Phenotype3 | Phenotype4 | Phenotype5 | Phenotype6 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: |
| 0.214992 | 0.224991 | NA | 1 | NA | -0.285427 |
| -0.974543 | -0.904542 | NA | 0 | NA | -2.333531 |
| 0.195909 | 0.195909 | NA | 1 | NA | 0.046818 |
| NA | NA | NA | NA | -0.65421 | NA |
| 0.367591 | NA | 0.25987 | NA | NA | -0.385427 |
| ... | ... | ... | ... | ... | ... |
| NA | -0.05247 | NA | 0 | NA | 0.720009 |

## Running build-in data <br>
```R
> library("PIXANT")
> data(demoData)
> head(demoData)
    Phenotype1     Phenotype2     Phenotype3     Phenotype4     Phenotype5     Phenotype6     Phenotype7     Phenotype8
1  -13.07400974   -12.47319542    0.588613265   -2.574249616    6.532921148   9.798340249   -0.765743678    0.059914511
2  -2.483186418   -1.393727485   -1.493932449   -1.648425872   -6.097557807   0.010440626   -0.534436706    0.264419892
3   3.318953906    3.761264938   -0.532976741    3.474977473   -4.591645791  -0.557710801    0.993769937   -0.688824646
4   8.660696095    8.388324362    0.746899563   -3.341995972   -4.485148772   -4.41061568   -0.640711383    2.931062553
5    13.5960061    13.10111325   -1.558483568   -3.551238007   -5.724882243  -2.582869264        NA        -1.551938619
> system.time(PIXANT.imp <- PIXANT(demoData, aimPhenName = 'Phenotype7', maxIterations = 20, maxIterations0 = 20,
>                                  num.trees = 100, initialLinearEffects = 0, errorTolerance = 0.001,
>                                  aimPhenMissingSize = 500, initialImputeType = 'random', refPhenThreshold = 0.3,
>                                  maxNum.refPhen = 7, SC.Threshold = 0.6, seed = 123))
```

## OUTPUT
***`PIXANT`*** returns total 9 objects of results including: <br>
***`ximp`***-The imputed data(a data frame);<br>
***`imputePhen`***-The impute phenotype, including imputed values and SC for each samples (a data frame);<br>
***`imputePhen.refPhen`***-A list of reference phenotypes for impute phenotype;<br>
***`imputePhen.accuracy`***-Imputation accuracy of imputed phenotype;<br>
***`imputePhen.r2`***-The r square of fitting between observed values and imputed values in the imputed phenotype;<br>
***`imputePhen.pValue`***-The p value of fitting between observed values and imputed values in the imputed phenotype;<br>
***`imputePhen.value`***-The missing data of settings in the imputed phenotype, including sample index, observed values and imputed values (a data frame);<br>
***`imputePhen.missingRate`***-Original missing rate of the imputed phenotypes;<br>
***`imputePhen.filter.missingRate`***-Missing rate after SC quality control of the imputed phenotype;<br>

```r
> str(PIXANT.imp)
List of 9
 $ ximp                         :'data.frame':	10000 obs. of  8 variables:
  ..$ Phenotype7: num [1:10000] -0.766 -0.534 0.994 -0.641 -0.576 ...
  ..$ Phenotype1: num [1:10000] -13.07 -2.48 3.32 8.66 13.6 ...
  ..$ Phenotype2: num [1:10000] -12.47 -1.39 3.76 8.39 13.1 ...
  ..$ Phenotype3: num [1:10000] 0.589 -1.494 -0.533 0.747 -1.558 ...
  ..$ Phenotype4: num [1:10000] -2.57 -1.65 3.47 -3.34 -3.55 ...
  ..$ Phenotype5: num [1:10000] 6.53 -6.1 -4.59 -4.49 -5.72 ...
  ..$ Phenotype6: num [1:10000] 9.7983 0.0104 -0.5577 -4.4106 -2.5829 ...
  ..$ Phenotype8: num [1:10000] 0.0599 0.2644 -0.6888 2.9311 -1.5519 ...
 $ imputePhen.accuracy          : num 0.954
 $ imputePhen.value             :'data.frame':	500 obs. of  3 variables:
  ..$ indSample: num [1:500] 2993 3057 2700 625 5240 ...
  ..$ Observed : num [1:500] -0.8424 -0.0676 -0.2014 4.5788 0.5633 ...
  ..$ Imputed  : num [1:500] -0.876 -0.461 -0.501 5.082 0.786 ...
 $ imputePhen.r2                : num 0.911
 $ imputePhen.pValue            : Named num 0
  ..- attr(*, "names")= chr "value"
 $ imputePhen.refPhen           :'data.frame':	7 obs. of  2 variables:
  ..$ ref.Phenotype: chr [1:7] "Phenotype4" "Phenotype8" "Phenotype5" "Phenotype2" ...
  ..$ correlation  : chr [1:7] "0.908317600022565" "-0.238486577355447" "0.215227900643262" "-0.210171431049365" ...
 $ imputePhen.missingRate       : num 0.184
 $ imputePhen                   :'data.frame':	10000 obs. of  2 variables:
  ..$ Phenotype7: num [1:10000] -0.766 -0.534 0.994 -0.641 -0.576 ...
  ..$ SC        : num [1:10000] 1 1 1 1 1 1 1 1 1 1 ...
 $ imputePhen.filter.missingRate: num 0.0471
 - attr(*, "class")= chr "PIXANT"
```

## A function to correct for covariates <br>
**Note:** Please attention that NAs are not allowed in the **cov.txt**,  and all individuals should be in the same order with phenotype file.
> `cov.txt`

| PC1 | PC2 | PC3 |
| :---: | :---: |  :---: |
| 0.21 | 0.22 | 0.36 |
| -0.97 | -0.86 | 0.75 |
| 0.19 | 0.59 | -0.66 |
| 0.77 | -0.24 | 0.85 |
| 0.36 | -0.25 | 0.27 |
| ... | ... | ... |
| 0.56 | -0.05 | 0.08 |

```R
> library("PIXANT")
> phen <- read.table("phenotype.txt", header=T)
> cov <- read.table("cov.txt", header=T)
> AdjustPhen <- PhenAdj(phen, cov)
```
## How to access help <br>
If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/GuLinLin/PIXANT/issues):point_left:, or send email to contact:<br>
 
:e-mail: **Lin-Lin Gu:** linlin-gu@outlook.com <br>
## Citation <br>
