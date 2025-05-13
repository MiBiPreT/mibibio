# `mibibio` Methods Overview

Methods implemented for comprehensive statistical, visual and meta-analysis of large, complex microbiome datasets derived from high-throughput sequencing (16S rRNA amplicon, shotgun metagenomics or metatranscriptomics). Most routines wrap—or reproduce programmatically—the functionality of the MicrobiomeAnalyst web suite and its underlying R packages, so that analyses can be executed reproducibly from the command line or in notebooks. The workflow includes:

## Filtering methods — removing uninformative or misleading features
### Low-information feature removal (`filter_features()`)
This is the first clean-up pass. The routine scans every taxon or gene (a “feature”) across all samples and discards any feature that (a) is zero in every sample, (b) appears in only one sample, or (c) sits below a user-defined minimum abundance. By stripping out such “dead weight” you cut matrix sparsity, speed up calculations, and avoid inflating richness estimates. The result is a minimally-cleaned phyloseq object that still preserves rare—but not useless—features for later α-diversity work.

### Prevalence filtering (`filter_prevalence()`)
Next, you can tighten the focus on the community “core”. For each feature the routine calculates the proportion of samples in which it has a non-zero count; anything below the chosen prevalence threshold (10 % by default) is removed. This retains taxa/genes that occur often enough to carry statistical weight, boosting power for group comparisons, ordination.

### Variance filtering (`filter_variance()`)
Finally, the routine computes a variability score (inter-quartile range or coefficient of variation) for every remaining feature. Those whose abundance scarcely changes between samples are unlikely to discriminate between experimental groups, so they are dropped if their variance falls below the chosen cutoff. The outcome is a lean dataset containing only features with the potential to explain biological differences and keeping the multiple-testing burden in check.

Suggested order: run the three filters sequentially — low-information → prevalence → variance.
Each step returns an updated phyloseq object that you can pass straight to the next filter or into diversity and differential analyses.

## Normalisation methods — putting all libraries on a common scale
### Total-Sum Scaling (TSS) — `normalise_tss()`
Counts in every sample are divided by that sample’s library size and multiplied by a constant (e.g. 1 000 000 or 100 %). The outcome is a compositional matrix of relative abundances whose rows all sum to the same value. This scale is ideal for stacked-bar taxonomic plots and straightforward α-diversity calculations.

### Centered Log-Ratio (CLR) transformation — `normalise_clr()`
After adding a tiny pseudocount to avoid zeros, each count is divided by the sample’s geometric mean and transformed with the natural logarithm. The CLR removes the unit-sum constraint of compositional data and yields continuous, centred values that work seamlessly with Euclidean distance metrics and PCA.

### Rarefaction — `rarefy_even_depth()`
Here every library is randomly subsampled without replacement down to a common depth you specify. Samples that fall short of that depth are removed. Rarefaction equalises library sizes, making richness and α-diversity comparisons fair, but it introduces sampling stochasticity — set a random seed to make runs reproducible

Normalisation routines return a normalised (or transformed) count matrix and an updated phyloseq object on the correct scale for the statistical engine you plan to use next (diversity profiling, ordination, edgeR, metagenomeSeq, …) or for visualisation exploratione by performing rarefaction curve analysis to allow users to visually assess the sequencing depth with regard to the number of OTUs detected.Rarefaction curve `PlotRarefactionCurve()`is usefull to check this trait of you data, it gives you some information on the quality of the sequencing. and could be usefull to determine the normalisation method. It basically gives you perspective on the sampling depth and the number of unique features that are found. Here you can see if the sequencing depth you did was sufficient to 'find all the unique sequences'.


## Community diversity profiling (Marker Data Profiling module)
The analysis can be performed at different taxonomic levels based on the available annotations.
### Alfa diversity metrics `alpha_diversity()`:
The alpha-diversity analyzing the species diversity within a local samples.  The results are plotted across samples `PlotAlphaData()` and are also summarized as box plots for each group `PlotAlphaBoxData()`. The corresponding statistical significance is estimated automatically using either parametric or non-parametric tests based on user selection. Users can also visualize abundance profiles at different taxonomic levels using a stacked area or stacked bar plot `PlotTaxaAundanceBar()` which gives an overall picture with regard to the community. The alpha_diversity analysis function supports following common diversity measures:
`Phylogenetic diversity (faiths)` for which a phyloseq object with a phylogenetic tree is required.`Dominance` alfa diversity(Berger-Parker) is usefull to determine different Dominance parameters of the microbiome indicating the proportion of the most abundant taxon relative to all others.`Dom_taxon`can find the most dominant feature identity. So it can determine which family is present the most in the samples. If this is for example Pseudomonas, Acinetobacter or other known degrading families, then this could be an indication of biodegradation.`informative diversity (Shannon)`reflects how many different microbes are in the sample and how evenly they are distributed within the sample and `Species richness (Observed)` measures how many different ASV/OTUs are present in each samples. This gives information with regard to the species richness.

### Beta_diversity metrics `beta_diversity()`:
The beta-diversity analysis supports five common distance measures (Bray–Curtis, UniFrac variants, Aitchison, etc) `PerformBetaDiversity()`. The results are presented as both 2D and 3D ordination plots based on principal coordinate analysis (PCoA) `PCoA3D.Anal()`or non-metric multidimensional scaling (NMDS). To help identify patterns or gain biological insight, the samples displayed on PCoA or NMDS plots can be colored based on the metadata (default), their alpha diversity measures, or the abundance levels of a particular feature they contain `PerformCategoryComp()`. 

## Predicting metabolic potentials and profiling functional diversity
 From the microbes whose whole genomes have been sequenced and annotated, the 16S rRNA data can be used to infer the metabolic potentials of the corresponding microbial species. In particular, the core microbiome analysis `CoreMicrobeAnalysis()` is usefull for exploratory analysis as some target groups could be identified from literature, it can indicate whether there are aerobic or anaerobic degraders present. The Tax4Fun function used for data annotated by SILVA database. It looks at the specific groups/traitstry and try to link functional traits to different taxa found in the different samples. So the Tax4Fun `Perform16FunAnot_mem` is used to link specific phylogenetic groups to different functions, instead of manually choosing certain groups which are 'known degraders'. The result is a table containing relative KO abundance levels. The KO profiles obtained from predictions or actually measured from shotgun metagenomics or metatranscriptomics can be used for functional diversity profiling `PlotFunAnotSummary()` based on the KEGG (pathways, modules or EC categories) or COG annotation systems.

