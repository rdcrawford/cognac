# ------------------------------------------------------------------------------
# Remove Outlier Genomes
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Some genes may contain fewer than the expected number of core genes because
# of a poor assembly or misspecification. This function iterately 
# ------------------------------------------------------------------------------

RemoveOutlierGenomes  = function( 
  geneMat, targetGeneNum, coreGeneThresh, isKeeper
  )
{
  # Count the frequency of the genes in the remaining genomes
  geneCount  = sapply( 1:ncol( geneMat ),
    function(j) sum( geneMat[ isKeeper, j ] > 0 )
    )
  
  # Find the top n genes in the dataset 
  topGenes = order( geneCount, decreasing = TRUE )[ 1:targetGeneNum ]
  
  # Of the remaining genomes, find the genome with the fewest number of 
  # core genes and remove it
  toTest = which( isKeeper )
  numMissing = sapply( toTest,
    function(i) sum( geneMat[ i, topGenes ] == 0 )
    )
  
  isMissingGenes = numMissing > 0
  if ( !TRUE %in% isMissingGenes ) return( isKeeper )
  
  isKeeper[ toTest[ which.max(numMissing) ] ] = FALSE
  return( isKeeper )
}

# ------------------------------------------------------------------------------
