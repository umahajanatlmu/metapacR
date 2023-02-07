
# metapacR 

<img src="https://github.com/umahajanatlmu/metapacR/blob/master/inst/figures/metapacR_sticker.png" alt="alt text" width="150" height="150" align="right">

[![stability][0]][1]

In-house package for metabolome analysis.
metapacR is an easy-to-use R package implementing a complete workflow for downstream analysis of targeted and untargeted metabolome data. Metabolome results can be imported into metapacR as a numerical matrix, allowing integration into current analysis frameworks. metapacR allows data inspection, normalization, univariate and multivariate analysis, displaying informative visualizations.

## Functions

|  No |  Function | Descriptions |
|---|---|---|
| 1 | banner | function to print tidy labels and separators |
| 2 | boxPlots | function to plot box and violin plot |
| 3 | compareDiamReduction | funcion to compare different diamentionality reduction methods, namely, PCA, OPLS, TSNE and UMAP |
| 4 |correlationPlot | function to plot intragroup correlation plot |
| 5 |diffAbundanceScore | function to analyse differential abundance score of enriched pathways |
| 6 |distributionPlot | function to analyse distribution of different classes of metabolites |
| 7 | enrichedNetwork | function to plot enriched annotated netwroks for the signficantly altered metabolic pathways |
| 8 | enrichedNetwork.KEGG | function to plot influence of metabolic alteration on enriched pathway (direction of change) |
| 9 | exportResults  | function to export results in color coded format in excel table |
| 10 | findMarkers  | function to select and plot most differentiated markers based on roc analysis |
| 11 | ImputeTransformScale | function to impute, transform and scale raw data |
| 12 | installScriptLibs | function to install listed R-packages from CRAN and Bioconductor  |
| 13 | leveneStat   | function to compute Levene's stats for homogeneity of variance across groups |
| 14 | lipidChainLengthCorrelation | function to compute correlation of chain length with abundance |
| 15 | lipidChainLengthDistribution | function to compute distribution of chain length per class |
| 16 | normalizeDat.binary | function to normalize data using fixed effect as well as mixed effects models for comparison of binary groups ie group vs rest |
| 17 | normalizeDat | function to normalize data using fixed effect as well as mixed effects models |
| 18 | piePlot | function to plot pieCharts |
| 19 | plotDiamReduction | function to plot groups and distribution of metsbolites in diamentionality reduction plot |
| 20 | plotMarkerViolin | function to plot violin plots of identified markers |
| 21 | plotMetaboliteAlteration | function to compute atered metabolites distribution per class |
| 22 | plotSplitViolin | function to plot split violin plots |
| 23 | rocPlots | function to compute and plot roc curves |
| 24 | shapiroTest | function to test normality |
| 25 | volcanoPlots | function to plot volcano plots |

<img src="https://github.com/umahajanatlmu/metapacR/blob/master/inst/figures/pipeline.png" alt="alt text" width="900" height="450">

## Package installation

### install bioconductor dependancies before package installation
```{r}
packages <- c("FELLA", "ggtree", "KEGGREST", "ropls")
BiocManager::install(packages)
```
### install package from source  
```{r}
devtools::install_github("umahajanatlmu/metapacR", ref="master", auth_token = "tokenString")
```
