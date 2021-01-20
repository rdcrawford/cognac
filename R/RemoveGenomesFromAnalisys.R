# ------------------------------------------------------------------------------
# Remove Outlier Genomes
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Some genes may contain fewer than the expected number of core genes because
# of a poor assembly or misspecification. This function iterately 
# ------------------------------------------------------------------------------

RemoveGenomesFromAnalisys  = function( geneEnv, toRemove )
{
  # Find the genomes to remove
  isOut = geneEnv$genomeNames %in% toRemove
  
  # Update the genomes and genes that are still being used for the analysis
  geneEnv$genomeNames = geneEnv$genomeNames[ !isOut ]
  geneEnv$geneMat     = geneEnv$geneMat[ !isOut, ]
  geneEnv$gfList      = geneEnv$gfList[ !isOut ]
  geneEnv$fastaFiles  = geneEnv$fastaFiles[ !isOut ]
  
  # Initialize a vector of gene IDs to remove
  toRemoveGenes = character()
  
  # For each cluster of genes find the genes that correspond to the genomes
  # that will be removed from the analysis
  for ( i in 1:length( geneEnv$genomeIdList ) )
  {
    # Find the names of the genomes to remove
    isOut = geneEnv$genomeIdList[[ i ]] %in% toRemove
    
    if ( TRUE %in% isOut )
    {
      addStart = length(toRemoveGenes ) + 1
      addRange = addStart:( addStart + ( sum( isOut ) - 1 ) )
      toRemoveGenes[ addRange ] = geneEnv$clustList[[ i ]][ isOut ]
      
      # Remove the gene ids
      geneEnv$genomeIdList[[ i ]] = geneEnv$genomeIdList[[ i ]][ !isOut ]
      geneEnv$clustList[[ i ]]    = geneEnv$clustList[[ i ]][ !isOut ] 
    }
  }
  
  # Remove the genes from the analysis
  isOutGene        = geneEnv$geneIds %in% toRemoveGenes
  geneEnv$geneIds  = geneEnv$geneIds[ !isOutGene ]
  geneEnv$geneSeqs = geneEnv$geneSeqs[ !isOutGene ]
}

# ------------------------------------------------------------------------------
