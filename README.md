# Basic run example

```
library(tidyverse)
peptide_rdf <- read_tsv("testdata/dummy_df.tsv")
peptide_mat <- peptide_rdf[, -1] %>% as.matrix()
protein_ids <- peptide_mat$Protein

pr$protein_rollup(protein_ids, peptide_mat)
```

# What is RROllup?

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

