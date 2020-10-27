# cognac: Core Gene Alignement Concatenation

# Description
The cognac function identifies shared genes to be used as phylogenetic markers within the input set of genomes. Marker genes are aligned individually with mafft and concatenated into a single alignment for downstream phylogenetic analysis.

# Install the package
install.packages("devtools")
devtools::install_github("rdcrawford/cognac")
library(cognac)
