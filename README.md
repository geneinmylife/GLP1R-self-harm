# The association of GLP1R activity with the risk of self-harm behaviours: a population-based study applying Mendelian randomization principles

In this study, we estimated the effect of GLP1R signalling on self-harm outcomes using Mendelian randomization (MR), colocalization and polygenic score association. 
The results highlight that the weight loss effect of GLP1R signalling was likely to increase the risk of self-harm behaviours and raised concern about mental health safety of using GLP1R agonists for weight control, potentially with higher risk in women. 

We provide the following code in "code" folder for reference:
* MR code in `"Mendelian_randomization.R"` which includes F-statistics and Steiger filtering for instruments strength and direction validation.
* polygenic score association code in `"polygenic-score.R"`.
* `"LD-check.R"` and `"colocalization.R"` for MR cuasal evidence validation.
* `"hyprcoloc.R"` for multi-trait colocalization detection.

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
install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)
```

Also, we provide the datasets in "data" folder for reference:
* Multi-omic QTL instruments, outcomelist and 8 extracted outcome datasets in `"MR"` folder.
* six prepared whole-genome datasets in 500kb around target SNP for colocalization and LDcheck, and moloc folder for multi-trait coloc in `"coloc_LDcheck"` folder.
* pQTLs for creating polygenic score in `"PGS"` folder.
* instruments for mediators in `"mediation"` folder.

Please download the whole project and set the download folder as the route for runing R code.
The whole analysis is based on R version 4.3.1 (2023-06-16).

Besides, We use the colocalization method: PWCOCO (https://github.com/jwr-git/pwcoco/) for double check. Please check the link for further instructions.
The hyprcoloc method is derived from Gwasglue2 package: Gwasglue2 (https://github.com/MRCIEU/gwasglue2). Please check the link for details.

# Publication
Any other related information could be referenced from our paper: The association of genetically predicted GLP1R signaling with the risk of self-harm behaviours: a population-based study using Mendelian randomization

# Data Sources
* eQTL/sQTL data: [GTEX Portal](https://gtexportal.org/home/).
* pQTL data: [deCODE](https://www.decode.com).
* meQTL data: [GODMC](http://mqtldb.godmc.org.uk/downloads).
* Self-harm data: [Zhang et al., 2021 BioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.15.435533v1.full).
* Suicide attempt data: [Docherty et al.](https://ajp.psychiatryonline.org/doi/10.1176/appi.ajp.21121266).
* T2D data: [Mahajan et al.](https://www.nature.com/articles/s41588-022-01058-3).
* HbA1c data: [Chen et al.](https://www.nature.com/articles/s41588-021-00852-9).








