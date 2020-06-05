# ------------------------------------------------------------------------------
# FilterMultiCopyGenes
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function identifies instances where a gene is present at greater than
# one copy and removes it from the dataset
# ------------------------------------------------------------------------------

FilterMultiCopyGenes = function( geneEnv, copyNumTresh )
{
  # Set the threshold for defining single copy genes -- default to multi-copy
  # in less than 1 in 1000 genomes
  if ( missing(copyNumTresh) ) copyNumTresh = 0.001

  # Find the fraction of genomes where
  nGenes   = length( geneEnv$clustList )
  nGenomes = length( geneEnv$genomeNames )
  dupFrac  = future.apply::future_sapply( seq(nGenes), function(i)
  {
    numDuplicates = sum( table( geneEnv$genomeIdList[[i]] ) != 1 )
    return( numDuplicates / nGenomes )
  })

  # Find any instances where there are genes with multiple copies
  isSingleCopy = dupFrac <= copyNumTresh

  # If there are genes with greater than one copy remove them from the
  # data
  if ( FALSE %in% isSingleCopy )
  {
    geneEnv$clustList    = geneEnv$clustList[ isSingleCopy ]
    geneEnv$genomeIdList = geneEnv$genomeIdList[ isSingleCopy ]
    geneEnv$geneMat      = geneEnv$geneMat[ , isSingleCopy ]

    cat(
      "  -- ", sum( !isSingleCopy ), " multi-copy genes were identified\n",
      "  -- After filtering there are ", sum(isSingleCopy), 
      " remaining genes\n",
      sep = ''
      )

    if ( !TRUE %in% isSingleCopy )
      stop( "No single copy genes were identified..." )

    # Find if there are any genes that have duplications below the
    # threshold
    hasDuplication = dupFrac[ isSingleCopy ] > 0
    if ( TRUE %in% hasDuplication )
    {
      # For each gene that has a gene duplication, remove one copy from the
      # cluster.
      for ( i in which( hasDuplication ) )
      {
        isNotDuplicated = !duplicated( geneEnv$genomeIdList[[ i ]] )
        geneEnv$clustList[[ i ]] =
          geneEnv$clustList[[ i ]][ isNotDuplicated ]
        geneEnv$genomeIdList[[ i ]] =
          geneEnv$genomeIdList[[ i ]][ isNotDuplicated ]
      }
    }
  }
}

# ------------------------------------------------------------------------------
