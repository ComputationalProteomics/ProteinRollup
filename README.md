# NOTE

This is a fresh software. Initially, I recommend you cross-check the results with InfernoRDN ([homepage](https://omics.pnl.gov/software/infernordn)). This one is a rewrite of the RRollup method found in [InfernoRDN](https://github.com/PNNL-Comp-Mass-Spec/InfernoRDN/blob/master/Rscripts/Rollup/RRollup.R#L55).

Settings to include:

* Presence percentage requirements for peptides (no requirement here, InfernoRDN requires 50% per default)

I wish:

* Optionally write the scaled peptides, and the scaled peptides without outlier matrices
* Optional profile plots of peptides and resulting protein

# Basic run examples from terminal

## With design matrix

If you have multiple annotation columns in your data matrix, you need to use a design matrix
to distinguish these

```
Rscript ProteinRollup.R \
    --rdf_fp testdata/dummy_df.tsv \
    --ddf_fp ProteinRollup/testdata/dummy_df_design.tsv \
    --sample_col sample \
    --protein_col Protein \
    --out_fp test_out.tsv
```

## With single-column annotation

If having an expression dataset with a single annotation column specifying the proteins, the following call can be executed.

```
Rscript ProteinRollup.R \
    --rdf_fp testdata/dummy_df.tsv \
    --out_fp test_out.tsv \
    --one_column_mode TRUE
```

# Basic run example from R

```
library(tidyverse)
source("ProteinRollup.R")

peptide_rdf <- read_tsv("testdata/dummy_df.tsv")
peptide_mat <- peptide_rdf[, -1] %>% as.matrix()
protein_ids <- peptide_mat$Protein

protein_rollup(protein_ids, peptide_mat)
```

# What is RRollup?

* A reference peptide is selected based on (1) fewest missing values and (2) if several peptides ties - the one with highest median intensity
* Ratios to reference are calculated for all. The ratio is simply subtracting the log-scale reference from each peptide
* For each non-reference peptide - the median distance to the reference are calculated
    * If lower than `min_overlap`, then set to zero in the median calculation (this is from Inferno, maybe to NA would make more sense? Now bias towards reference)
* The other peptides are shifted to the reference median. Intrapeptide variances are retained.
* Outlier removal is performed
    * `minPs` specifies how many peptides needs to be present with value in a sample for the outlier calculation to be carried out (default 5)
    * `pvalue` specify how certain the Grubb outlier test should be for removal
* Peptides can optionally be centered (TODO: What consequence?)
* A rollup score is calculated (meaning what? why is not outlier removed peptides used here?)
* The final protein matrix is returned

# Key parameters for rollup

* `min_overlap` - How many values should be in common between reference and other peptides for the other peptides to be considered

