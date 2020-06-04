# ------------------------------------------------------------------------------
# Calc Num Core Genes
# 2020/01/31
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function takes the environment with the gene presence absence matrix
# and a threshold for the number of core genes as a fraction of the matrix
# and calcuates the number of genes in the matrix that are present
# at or above that threshold.
# ------------------------------------------------------------------------------

CalcNumCoreGenes = function( genePtr, coreGenomeThresh, isKeeper )
{
  # Define the minimium threshold for number of core genes
  coreCountVal = round(  sum(isKeeper) * coreGenomeThresh, 0 )

  # Recalculate the number of core genes
  isCoreGene = sapply( 1:ncol(genePtr$geneMat),
    function(j) sum( genePtr$geneMat[ isKeeper, j ] > 0 ) >= coreCountVal
    )

  # Returnt the number of core genes
  return( isCoreGene )
}

# ------------------------------------------------------------------------------
