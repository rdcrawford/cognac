# cognac: Core Gene Alignment Concatenation

## Description

Cognac is an R package for for generating concatenated gene alignments. The main function identifies shared genes to be used as phylogenetic markers within the input set of genomes. Marker genes are aligned individually with mafft and concatenated into a single alignment for downstream phylogenetic analysis. The algorithm uses mutli-threading and several algorithmic tricks to improve speed and efficiency, which make cognac capable of generating core-gene alignments for thousands of genomes in only a few hours. 

For more information on the algorithm and bench marking, see our preprint posted on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.15.340901v1).

## Install the package

```
install.packages("devtools")
devtools::install_github("rdcrawford/cognac")
library(cognac)
```

[Mafft](https://mafft.cbrc.jp/alignment/software/) and [cd-hit](https://github.com/weizhongli/cdhit) must be in your path. 


## Usage

The most basic command to use cognac is to supply a directories containing whole genome seqences in fasta files and genome annotations in the form of gff. Multithreading is available at multiple steps. The number of threads to be used can be supplied as argument. By default, all available threads are used. 

```
algnEnv = cognac(
  fastaDir   = "path/to/your/fasta/files/",
  featureDir = "path/to/your/gff/files/",
  threadVal  = 4
  )
```

Alternatively you can supply character vectors with the paths to the fasta files and the gff files.

```
fastaFiles   = c( genome1.fasta, genome2.fasta, genome3.fasta )
featureFiles = c( genome1.gff, genome2.gff, genome3.gff )

algnEnv = cognac(
 fastaFiles   = fastaFiles,
 featureFiles = featureFiles,
 threadVal    = 4
 )
```

Optionally, you can generate a codon aware nucleotide alignment from the amino acid alignment by mapping the nucleotide sequence from each gene back to the amino acid alignment position-wise. This may be useful for differentiating closely related genomes by leveraging the degeneracy in the codon code. 


```
algnEnv = cognac(
  fastaDir   = "path/to/fasta/_files/",
  featureDir = "path/to/gff_files/",
  mapNtToAa  = TRUE,
  threadVal  = 4
  )
```


The output for cognac is an [environment](http://adv-r.had.co.nz/Environments.html), where multiple objects can be stored and accessed. By Default, cognac produces two objects: the path to the concatenated gene alignment and the meta-data on the genes included. If the nucleotide alignment was requested the path to the alignment is included in the environment. The meta-data includes columns: the gene description, the comma eliminated gene IDs, the positions of the partitions in the amino alignment, and the positions of the partitions in the nucleotide alignment, if requested. 

```
cat( algnEnv$aaAlgnPath )
cat( algnEnv$ntAlgnPath )
cat( head( algnEnv$geneData ) )
```

We offer the option to create a neighbor joining tree within the cognac function. This a useful method for generating trees based off of the genetic distances, especially for large data sets where other methods may be too computationally intensive. However, this method may not be appropriate for distantly related sequences.

```
algnEnv = cognac(
  fastaDir      = "path/to/your/fasta/files/",
  featureDir    = "path/to/your/gff/files/",
  njTree        = TRUE
  )
  
ape::plot.phylo( algnEnv$njTree )
```

While neighbor joining trees are useful, often times higher resolution methods are required. Maximum likelihood based methods, such as [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/), are more accurate. Alternatively, [FastTree](http://www.microbesonline.org/fasttree/) is an approximate maximum likelihood based method which is easy to use and 100-1000 times faster.

``` 
system( paste( "FastTree <", algnEnv$aaAlgnPath, "> cognac_fastTree.tre" ) )
```


Cognac requires genome annotations in gff format. There are many tools available for generating annotations such as [RAST](https://docs.patricbrc.org/cli_tutorial/rasttk_getting_started.html), [prokka](https://github.com/tseemann/prokka), or [prodigal](https://github.com/hyattpd/Prodigal).

If gff files are not available for your genomes, we have provided a function to generate them with RAST using the command line interface -- installation instructions can be found [here](https://docs.patricbrc.org/cli_tutorial/index.html#installing-the-cli-release). 

```
gffFiles = sapply( fastaFiles, AnnotateGenome, outDir = "path/to/gff_files/" )

```