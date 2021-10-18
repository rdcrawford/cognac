# ------------------------------------------------------------------------------
# Align Gene Sequences
# 2019/12/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Generates an alignment for a single genes. The index of the orthologous gene
# in the list with the cd-hit output, and the environment with the parsed
# gene sequences and meta-data. The gene sequenes are found, identical
# sequences are identified, and only the identical sequences are aligned with
# mafft. The new alignment is then read in and the
# ------------------------------------------------------------------------------

AlgnGeneSeqs = function(
  geneEnv,  # Environment with the parsed gene sequences
  clustIdx, # Vector with unique gene ids to align
  algnDir,  # Path to the fasta file to write gene sequences to
  mafftOpts # Optional. Arguments for mafft. mafft [mafftOpts] in > out
  )
{
  # Set the path to the output files for the genes and the gene alignment
  geneFaPath = paste0( algnDir, "gene_cluster_", clustIdx, ".fasta" )
  algnPath   = paste0( algnDir, "gene_cluster_", clustIdx, "_mafft.fasta" )

  # Find the genes in the corresponding to genomes which are still in the
  # analysis
  algnGenes    = geneEnv$clustList[[ clustIdx ]]
  isAlgnGenome = geneEnv$genomeIdList[[ clustIdx ]] %in% geneEnv$genomeNames
  algnGenes    = algnGenes[ isAlgnGenome ]
  isAlgnGene   = geneEnv$geneIds %in% algnGenes

  # For each unique sequence, find the identical sequences and
  # store them in a list.
  identList = FindIdenticalGenes(
    geneEnv$geneSeqs[ isAlgnGene ], geneEnv$geneIds[ isAlgnGene ]
    )

  # If there is no variation in the gene, it doesnt matter. Just return an
  # empty character vector
  if (length(identList) == 1)
  {
    algn = vector( "character", length(geneEnv$genomeNames) )
    return( list( algn ) )
  }

  # If there are no duplication, align all of the genes
  isGeneRep  = geneEnv$geneIds[ isAlgnGene ] %in% names(identList)
  toAlgnIdxs = which(isAlgnGene)[ isGeneRep ]

  # Write the input fasta file
  sink( geneFaPath )
  for ( i in toAlgnIdxs )
  {
    cat('>', geneEnv$geneIds[i], '\n', geneEnv$geneSeqs[i], '\n', sep = '')
  }
  sink()

  # Generate the mafft command
  if ( missing( mafftOpts ) )
  {
    mafftOpts = "--retree 2 --maxiterate 2 --quiet"
  }
  mafftCmd = paste( "mafft", mafftOpts, geneFaPath, '>', algnPath )

  # Run Mafft
  system( mafftCmd )
  
  # Read in the alignment
  algn = ParseFasta( algnPath )

  # If there were any duplicated genes
  for ( i in 1:length(identList) )
  {
    nGenes = length(identList[[i]])
    if ( nGenes )
    {
      algnVecLen = length( algn ) + 1
      newSeqRange = algnVecLen:(algnVecLen + nGenes - 1)
      algn[ newSeqRange ] = rep( algn[i], nGenes )
      names(algn)[ newSeqRange ] = identList[[i]]
    }
  }

  # Get the genome name of origin for each gene
  algnGenomeIds = sapply( names(algn), GetGenomeId, USE.NAMES = FALSE )

  # Sorted order of the genomes
  genomeIdOrder = vector("integer", length(geneEnv$genomeNames))

  # Sort the alignment by the genome names
  # Iterate over all of the genomes and find the position of the
  # genome in the vector of genome Ids
  for ( j in seq( length(geneEnv$genomeNames) ) )
  {
    # Find the alignment header corresponding to this genome
    isCurGenome = algnGenomeIds == geneEnv$genomeNames[j]

    # If this genome is missing, add a string of '-' the length of the
    # gene to the alignment vector
    if ( !TRUE %in% isCurGenome )
    {
      idx = length(algn) + 1
      algn[ idx ] = paste( rep( '-', nchar(algn[1]) ), collapse = '' )
      algnGenomeIds[ idx ] = geneEnv$genomeNames[ j ]

    } else {

      # Find the position of the genome in the current vector
      idx = which( isCurGenome )
    }

    # Set the order to the position of this genome
    genomeIdOrder[ j ] = idx
  }

  algn = algn[ genomeIdOrder ]
  names(algn) = algnGenomeIds[ genomeIdOrder ]

  # Output the path to mafft alignment
  return( list(algn) )
}

# ------------------------------------------------------------------------------
