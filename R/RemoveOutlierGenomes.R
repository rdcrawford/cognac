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
  counter = 1
  # Iterativly remove genomes until there is the minimium number of core genes
  repeat 
  {
    # Count the frequency of the genes in the remaining genomes
    geneCount  = sapply( 1:ncol(geneMat),
      function(j) sum( geneMat[ isKeeper, j ] > 0 )
      )
    
    # Find the top n genes in the dataset 
    topGenes = order( geneCount, decreasing = TRUE )[ 1:targetGeneNum ]
    
    # Find the minimum frequency to be considered a core gene
    coreCountVal = round( sum(isKeeper) * coreGeneThresh, 0 )
    
    # Get the total count of core genes 
    numCoreGenes = sum( geneCount >= coreCountVal )
    
    # If the target number of genes has been reached, 
    if ( numCoreGenes >= targetGeneNum ) break
    
    # Of the remaining genomes, find the genome with the fewest number of 
    # core genes and remove it
    toTest = which( isKeeper )
    numMissing = sapply( toTest,
      function(i) sum( geneMat[ i, topGenes ] == 0 )
      )
    isKeeper[ toTest[ which.max(numMissing) ] ] = FALSE
    
    if ( sum(isKeeper) == 0 || !TRUE %in% ( geneCount > 0 ) )
      stop("There are no genes in the matrix...")
  }
  
  return( isKeeper )
}

# ------------------------------------------------------------------------------
