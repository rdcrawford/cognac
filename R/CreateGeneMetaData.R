# ------------------------------------------------------------------------------
# Create Gene Meta-Data
# 2019/07/15
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function creates a dataframe with information on the cd-hit clusters
# including, the gene annotation from rast, the length of the gene,
# the cd hit cluster, and whether the gene is single copy in every geneome
# ------------------------------------------------------------------------------

CreateGeneMetaData = function( geneEnv, revTranslate )
{
  # Find a representitive member from each cluster
  nGenes = length( geneEnv$clustList )
  
  # Create a character vector with the gene Ids of all of the genes 
  # in each cluster 
  clGeneIds = sapply( seq(nGenes),
    function(i) paste( geneEnv$clustList[[i]], collapse = ',' )
    )

  # Look up the gene description
  geneDescr = future.apply::future_sapply( seq(nGenes), function(i)
  {
    # Find the genome name of the cluster representitive
    gIdx = which( geneEnv$genomeNames == geneEnv$genomeIdList[[i]][1] )

    # Look up the row index of the gene in the data-frame of
    # parsed genome features
    rIdx = which( geneEnv$gfList[[ gIdx ]]$featId == geneEnv$clustList[[i]][1] )

    # Return the description of the gene
    return( geneEnv$gfList[[ gIdx ]]$description[ rIdx ] )
  })

  # Create a dataframe with the data on all of selected core genes
  geneData = cbind.data.frame(
    geneDescr,
    clGeneIds,
    GetGenePartitions( geneEnv, revTranslate ),
    stringsAsFactors = FALSE
    )
  
  return( geneData )
}

# ------------------------------------------------------------------------------
