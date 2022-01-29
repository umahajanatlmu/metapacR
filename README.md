
# metapacR 

<img src="https://github.com/umahajanatlmu/metapacR/blob/master/inst/figures/metapacR_sticker.png" alt="alt text" width="150" height="150" align="right">

In-house package for metabolome analysis.
metapacR is an easy-to-use R package implementing a complete workflow for downstream analysis of targeted and untargeted metabolome data. Metabolome results can be imported into metapacR as a numerical matrix, allowing integration into current analysis frameworks. metapacR allows data inspection, normalization, univariate and multivariate analysis, displaying informative visualizations.

## Functions

|  No |  Function | Descriptions |
|---|---|---|
| 1 | banner | function to print tidy labels and separators |
| 2 | boxPlots | function to plot box and violin plot |
| 3 | compareDiamReduction | funcion to compare different diamentionality reduction methods, namely, PCA, OPLS, TSNE and UMAP. |
| 4 |normalizeDat   | Function for normalization and pairwise anova analysis of metabolome data (fixed effects as well as mixed effect models)    |
| 5 |normalizeDat.binary   | Function for normalization and group vs rest anova analysis of metabolome data (fixed effects as well as mixed effect models)    |
| 6 |leveneStat   |Function to test variance |
| 7 |shapiroTest   |Function to test data normality |
| 8 |compareDiamReduction  |Function to test different methods of diamentionality reduction  |
| 9   |plotDiamReduction  |Fuction to plot method of choice of diamentionality reduction  |
| 10  |boxPlots  |Function to plot boxplots |
| 11  |rocPlots   |Function to plot roc curves |
| 12  |piePlots   |   |
| 13  |correlationPlots   |   |
| 14  |distributionPlots  |   |
| 15  |volcanoPlots   |   |
| 16  |plotMetaboliteAlteration   |   |
| 17  |findMarkers   |   |
| 18  |exportResults   |   |
| 19  |lipidChainLengthCorrelation   |   |
| 20  |lipidChainLengthDistribution   |   |
| 21  |compare.common.themes  |   |
| 22  |erichmentScore.KEGG |   |
| 23  |diffAbundanceScore  |   |
| 24  |enrichedNetwork  |   |
| 25  |plotSplitViolin |   |
| 26  |plotMarkerViolin |   |

## Package installation
 ```{r}
devtools::install_github("umahajanatlmu/metapacR", ref="master", auth_token = "tokenString")
 ```
