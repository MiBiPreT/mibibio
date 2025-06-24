---
title: "Microbiome Data Analysis Pipeline Template"
output: html_document

Originally adapted by **Simon Vandersanden**, simon.vandersanden@uhasselt.be  
Last edits: 24/03/2025

In the following you can read all the steps to create the microbiome data analysis pipeline for your project:

## 1. Setup Working Directory and Environment
This initial step configures your project’s working folder and loads all required R packages.

### Tasks:
- Set the knitr root directory.
- Load necessary libraries.

#### Example:
```r
# Set project root directory for file I/O
knitr::opts_knit$set(
  root.dir = "C:/Amplicon_sequencing_analysis/.../input"
)

# Load libraries
library(phyloseq)
library(dplyr)
library(shiny)
library(crosstalk)
```

## 2. Data Ingestion (Import Raw Microbiome Data)
Gather your data either as a single combined object or as separate tables.

### Tasks:
- Read in a saved phyloseq object, or
- Import OTU/ASV table, taxonomy table, metadata, and phylogenetic tree separately.

#### Example:
```r
# As a phyloseq object
ps <- readRDS("ps_IdTax_noncontam_NoMito_No_Chloro.rds")

# Or as separate files
otu_table   <- read.csv("otu_table.csv")
tax_table   <- read.csv("tax_table.csv")
meta_table  <- read.table("metadata.txt", header=TRUE, sep="\t", row.names=1)
phylo_tree  <- read.tree("phylo_tree.nwk")

phyloseq_obj <- phyloseq(
  otu_table(otu_table, taxa_are_rows=FALSE),
  tax_table(tax_table),
  sample_data(meta_table),
  phy_tree(phylo_tree)
)
```

## 3. Data Filtering (Quality Control)
Filter out unwanted sequences and select sample subsets.

### Tasks:
- Remove contaminants (e.g., mitochondria, chloroplasts).
- Subset samples based on metadata criteria.

#### Example:
```r
# Filter taxa
ps_filt <- subset_taxa(ps, !Family %in% c("Mitochondria", "Chloroplast"))

# Subset samples (e.g., only "soil" samples)
ps_soil <- subset_samples(ps_filt, sample_type == "soil")
```

## 4. Data Normalization
Normalize counts to account for variable sequencing depth and sparsity.

### Tasks:
- Rarefaction or conversion to relative abundance.
- Document normalization choices.

#### Example:
```r
# Rarefy to even depth
ps_rarefied <- rarefy_even_depth(ps_soil, sample.size = min(sample_sums(ps_soil)))

# Transform to relative abundance
ps_relabun <- transform_sample_counts(ps_soil, function(x) x / sum(x))
```

## 5. Exploration and Visualization
Generate initial plots to assess data quality and composition.

### Examples:
- Rarefaction curves  
  ```r
  PlotRarefactionCurve(ps_rarefied)
  ```
- Taxa abundance barplots  
  ```r
  PlotTaxaAundanceBar(ps_relabun, tax_level="Phylum")
  ```

## 6. Alpha Diversity Analysis
Compute within-sample diversity metrics.

### Tasks:
- Shannon diversity, Observed richness.
- Phylogenetic diversity (Faith’s PD).
- Dominance indices (Berger–Parker).

#### Example:
```r
# Shannon index
#A value that reflects how many different microbes are in the sample and how evenly they are distributed within the sample
alpha_diversit_id_tax <-  alpha(Phyloseq_Object, index = c('shannon')
#help(alpha)
write.table(alpha_diversit_id_tax, file = "Shannon_diversit_id_tax.csv", row.names = T, sep = ";", dec = ",")

# Faith's PD
#here you need a phyloseq object with a phylogenetic tree.
#Faiths: the sum of the branch lengths on the minimum spanning tree linking a set of terminal taxa to the root 
#PD=i∈S∑​Li​
#where SSS represents the set of branches in the minimum spanning phylogenetic tree connecting all observed taxa in a sample, and LiL_iLi​ is the branch length of each branch.

library(btools)
phylogenetic_Diversity <-estimate_pd(Mibi_Phylo_Object)
View(phylogenetic_Diversity)
estimate_pd(Phyloseq_Object)
help("estimate_pd")
write.table(phylogenetic_Diversity, file = "phylogenetic_Diversity.csv", row.names = T, sep = ";", dec = ",")

# Berger-Parker dominance
# it is usefull to determine different 'Dominance parameters of your microbiome
#Berger-Parker metric, which indicates the proportion of the most abundant taxon relative to all others
library(microbiome)
dom_idtax <-dominance(Phyloseq_Object_idtax, index = "all", rank = 1, aggregate = T)
dominance
help("dominance")
#Here there are several dominance metrics calculated for your different samples, information can be found in the help tab. I personally use Berger-Parker metric, others are also interesting. 
#Another interesting feature which is usefull is the following, where you can find the most dominant feature identity. So you can see for example at the level of family, which family is present the most in the samples. If this is for example Pseudomonas, Acinetobacter or other known degrading families, then this could be an indication of bio degradation.
Dom_taxon <- dominant(Phyloseq_Object, level = 'family')

# Species richness (Observed)
#A measure of how many different ASV/OTUs are present in each samples. This gives information with regard to the species richness.
alpha_diversit_id_tax <-  alpha(Phyloseq_Object, index = c('observed')
#help(alpha)
write.table(alpha_diversit_id_tax, file = "Shannon_diversit_id_tax.csv", row.names = T, sep = ";", dec = ",")

```

## 7. Beta Diversity Analysis
Assess differences between samples.

### Tasks:
- Calculate distance matrices (Bray-Curtis, UniFrac).
- Perform ordinations (PCoA, NMDS).
- Statistical tests (PERMANOVA/Adonis).

#### Example:
```r
mbSet<-PerformBetaDiversity(mbSet, "beta_diver_0","PCoA","bray","expfac","SC_number","none","OTU","","Shannon", "yes", "adonis", "png", 72, "default", "true");
mbSet<-PCoA3D.Anal(mbSet, "PCoA","bray","OTU","expfac","SC_number","","Shannon","beta_diver3d_0.json")
mbSet<-PerformCategoryComp(mbSet, "OTU", "adonis","bray","SC_number","true");
#Example research question: How much do the different communities look like one another on a phylogenetic level, based on the presence/absence of different features ?
mbSet<-PCoA3D.Anal(mbSet, "PCoA","unifrac","OTU","expfac","SC_number","","Shannon","beta_diver3d_2.json")
mbSet<-PerformCategoryComp(mbSet, "OTU", "adonis","unifrac","SC_number","true");
#Example research question: How much do the different communities look like one another on a phylogenetic level, taking into account the abundance of the different features?
mbSet<-PCoA3D.Anal(mbSet, "PCoA","wunifrac","OTU","expfac","SC_number","","Shannon","beta_diver3d_3.json")
mbSet<-PerformCategoryComp(mbSet, "OTU", "adonis","wunifrac","SC_number","true")
```

## 8. Core Microbiome and Functional Prediction
Identify core taxa and predict functions (e.g., via Tax4Fun).

### Tasks:
- Determine core taxa by prevalence and abundance thresholds.
- Link taxa to functional traits or pathways.

#### Example:
```r
# Core microbiome (>=10% prevalence, >=1% abundance)
core_taxa <- core_members(ps_relabun, detection=0.01, prevalence=0.1)

# Tax4Fun functional prediction
mbSet <- Perform16FunAnot_mem(mbSet, "SILVA", "qi_silva", "gg13")
```

## 9. Merging Phylogenetic Tree (Optional)
Combine or update the phylogenetic tree if it’s maintained separately.

### Example:
```r
merged_tree <- merge_phyloseq(ps_relabun, phy_tree(new_tree))
```

## 10. Output and Reporting
Save results, export tables and plots, and document your findings.

### Tasks:
- Write diversity metrics to CSV.
- Save plots (PNG, PDF).
- Knit final RMarkdown report.

#### Example:
```r
write.csv(alpha_shannon, "alpha_shannon.csv")
ggsave("beta_pcoa_plot.png")
```