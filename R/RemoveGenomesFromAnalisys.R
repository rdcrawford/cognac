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
  
  for ( i in 1:length(geneEnv$genomeIdList) )
  {
    isOut = geneEnv$genomeIdList[[ i ]] %in% toRemove
    geneEnv$genomeIdList[[ i ]] = geneEnv$genomeIdList[[ i ]][ !isOut ]
    geneEnv$clustList[[ i ]]    = geneEnv$clustList[[ i ]][ !isOut ]
  }
}

# ------------------------------------------------------------------------------
