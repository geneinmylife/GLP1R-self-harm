# The association of GLP1R activity with the risk of self-harm behaviours: a population-based study applying Mendelian randomization principles

In this study, we estimated the effect of GLP1R activity (proxied by genetic predictors of gene expression, splicing, DNA methylation and protein levels) on self-harm outcomes using Mendelian randomization (MR), colocalization and polygenic score association. 
The results highlight that GLP1R-related traits show evidence of increasing the risk of self-harm behaviours, highlighting the need to further investigate the mental health consequences of using GLP1R agonists for weight control. 

We provide the following code in "code" folder for reference:
* MR code in `"Mendelian_randomization.R"` which includes F-statistics and Steiger filtering for instruments strength and direction validation.
* `"bidirectional-mr.R"` for bidirecitonal MR of self-harm traits (exposures) on protein levels of GLP1R (outcome).
* polygenic score association code in `"polygenic-score.R"`.
* `"LD-check.R"` and `"colocalization.R"` for MR cuasal evidence validation.
* `"hyprcoloc.R"` for multi-trait colocalization detection.
* `"phenospd.R"` for independent test number estimation.

To start using the code, you need to install `TwoSampleMR`, `Mendelianrandomization`, `coloc` and `ieugwasr` package.

```key
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

```key
devtools::install_github("mrcieu/ieugwasr")
```
For possible multi-trait colocalization, you need to install `gwasglue2`, `VariantAnnotation`, `gwasvcf` and `hyprcoloc` package. 

```key
install.packages("remotes")
remotes::install_github("MRCIEU/gwasglue2")
```

```key
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
```

```key
remotes::install_github("mrcieu/gwasvcf")
```

```key
remotes::install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)
```

Also, we provide datasets in the "data" folder for reference:
* Multi-omic QTL instruments, outcomelist and 8 extracted outcome datasets in `"MR"` folder.
* six prepared whole-genome datasets in 500kb around target SNP for colocalization and LDcheck, and moloc folder for multi-trait coloc in `"coloc_LDcheck"` folder.
* pQTLs for creating polygenic score in `"PGS"` folder.
* instruments for mediators in `"mediation"` folder.
* phenotypic correlation matrix for independent test estimation in `"phenospd"` folder.

Please download the whole project and set the download folder as the route for runing R code.
The whole analysis is based on R version 4.3.1 (2023-06-16) and the project is under MIT license.

We use the colocalization method: [PwCoCo](https://github.com/jwr-git/pwcoco/) to validate our colocalization findings. Please check the link for further instructions.
We used the [HyPrColoc](https://www.nature.com/articles/s41467-020-20885-8) method implemented in [```gwasglue2``` package](https://github.com/MRCIEU/gwasglue2). 

Besides, for some of the resultes shown in the supplementary table, we made the beta transformation of MR estimates when the outcome GWAS was conducted using a linear regression model. So please make the transformation with the following formula after runing the code to match the results in the table.


<image src="https://github.com/ling710/backup/blob/main/pic/Picture1.png" width="250"/>

While the variance was estimated using the following formula:  

<image src="https://github.com/ling710/backup/blob/main/pic/Picture2.png" width="450"/>

k:  represents the proportion of cases in the sampled population；
p:  effect allele frequency of all participants；                         
p0: effect allele frequency of controls；
b1: effect size of the linear regression.


# Publication
Any other related information can be referenced from our paper: The association of genetically predicted GLP1R activity with the risk of self-harm behaviours: a population-based study using Mendelian randomization

# Data Sources
* eQTL/sQTL data: [GTEX Portal](https://gtexportal.org/home/).
* pQTL data: [deCODE](https://www.decode.com).
* meQTL data: [GODMC](http://mqtldb.godmc.org.uk/downloads).
* Self-harm data: [Zhang et al., 2021 BioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.15.435533v1.full).
* Suicide attempt data: [Docherty et al.](https://ajp.psychiatryonline.org/doi/10.1176/appi.ajp.21121266).
* T2D data: [Mahajan et al.](https://www.nature.com/articles/s41588-022-01058-3).
* HbA1c data: [Chen et al.](https://www.nature.com/articles/s41588-021-00852-9).








