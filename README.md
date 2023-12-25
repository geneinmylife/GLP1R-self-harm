# Effect of GLP1R siganlling and self-harm using Mendelian randomization
In this study, we estimated the effect of GLP1R signalling on self-harm outcomes using Mendelian randomization (MR), colocalization and polygenic score association. 
The results highlight that the weight loss effect of GLP1R signalling was likely to increase the risk of self-harm behaviours and raised concern about mental health safety of using GLP1R agonists for weight control, potentially with higher risk in women. 

We upload the discovery MR code in `"Mendelian_randomization.R"` and the polygenic score association code in `"polygenic-score.R"`

To start using the code, you need to install `TwoSampleMR`, `Mendelianrandomization` and `ieugwasr` package.

```key
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

```key
devtools::install_github("mrcieu/ieugwasr")
```

we include the following functions in R code for reference:
* colocalization in `colocalization.R`
* LD check in `LD-check.R`

The colocalization method: PWCOCO (https://github.com/jwr-git/pwcoco/). 

# Publication
Any other related information could be referenced from our paper: Proteome-wide Mendelian randomization in global biobank meta-analysis reveals multi-ancestry drug targets for common diseases (https://www.medrxiv.org/content/10.1101/2022.01.09.21268473v1)

# Data Sources
* QTL data: [Global Biobank Meta-analysis Initiative (GBMI)](https://www.globalbiobankmeta.org/).
* Self-harm data: [Zhang et al., 2021 BioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.15.435533v1.full).
* Suicide attempt data [Zhang et al., 2021 BioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.15.435533v1.full).


