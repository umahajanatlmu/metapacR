
<img src="https://github.com/umahajanatlmu/metapacR/blob/master/inst/figures/metapacR_sticker.png" alt="alt text" width="200" height="200">

# metapacR

In-house package for metabolome analysis. It consist of following functions:

|  No |  Function | Descriptions |
|---|---|---|
| 1 | InstallScriptLibs   |  Function for installing and loading listed libraries |
| 2 | banner | Function of tidy titles and banners within scripts  |
| 3 | ImputeTransformScale |  Function for imputation, transformations and scaling of metabolome data |
| 4 |normalizeDat   | Function for normalization and anova analysis of metabolome data (fixed effects as well as mixed effect models)    |
| 5 |normalizeDat.binary   |   |
| 6 |leveneStat   |   |
| 7 |shapiroTest   |   |
| 8 |compareDiamReduction  |   |
| 9  |plotDiamRecuction  |   |
| 10  |boxPlots  |   |
| 11  |rocPlots   |   |
| 12  |piePlots   |   |
| 13  |correaltionPlots   |   |
| 14  |distributionPlots  |   |
| 15  |volcanoPlots   |   |
| 16  |plotMetaboliteAlteration   |   |
| 17  |findMarkers   |   |
| 18  |exportResults   |   |
| 19 |lipidChainLengthCorrelation   |   |
| 20  |lipidChainLengthDistribution   |   |
| 21  |compare.common.thermes  |   |

## Package installation
 ```{r}
devtools::install_github("umahajanatlmu/metapacR", ref="master", auth_token = "tokenString")
 ```
