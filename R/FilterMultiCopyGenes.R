# ------------------------------------------------------------------------------
# FilterMultiCopyGenes
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function identifies instances where a gene is present at greater than
# one copy and removes it from the dataset
# ------------------------------------------------------------------------------

FilterMultiCopyGenes = function( genePtr, copyNumTresh )
{
  # Set the threshold for defining single copy genes -- default to multi-copy
  # in less than 1 in 1000 genomes
  if ( missing(copyNumTresh) ) copyNumTresh = 0.001

  # Find the fraction of genomes where
  nGenes   = length( genePtr$clustList )
  nGenomes = length( genePtr$genomeNames )
  dupFrac  = sapply( seq(nGenes), function(i)
  {
    numDuplicates = sum( table( genePtr$genomeIdList[[i]] ) != 1 )
    return( numDuplicates / nGenomes )
  })

  # Find any instances where there are genes with multiple copies
  isSingleCopy = dupFrac <= copyNumTresh

  # If there are genes with greater than one copy remove them from the
  # data
  if ( FALSE %in% isSingleCopy )
  {
    cat("  -- Before filtering there were", nGenes, "genes\n")
    genePtr$clustList    = genePtr$clustList[ isSingleCopy ]
    genePtr$genomeIdList = genePtr$genomeIdList[ isSingleCopy ]
    genePtr$geneMat      = genePtr$geneMat[ , isSingleCopy ]

    cat("  -- After filtering there were", sum(isSingleCopy), "genes\n")

    if ( !TRUE %in% isSingleCopy )
      stop("No single copy genes were identified...")

    # Find if there are any genes that have duplications below the
    # threshold
    hasDuplication = dupFrac[ isSingleCopy ] > 0
    if ( TRUE %in% hasDuplication )
    {
      # For each gene that has a gene duplication, remove one copy from the
      # cluster.
      for ( i in which( hasDuplication ) )
      {
        isNotDuplicated = !duplicated( genePtr$genomeIdList[[ i ]] )
        genePtr$clustList[[ i ]] =
          genePtr$clustList[[ i ]][ isNotDuplicated ]
        genePtr$genomeIdList[[ i ]] =
          genePtr$genomeIdList[[ i ]][ isNotDuplicated ]
      }
    }
  }
}

# ------------------------------------------------------------------------------
