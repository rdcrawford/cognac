# ------------------------------------------------------------------------------
# Create Gene Meta-Data
# 2019/07/15
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function creates a dataframe with information on the cd-hit clusters
# including, the gene annotation from rast, the length of the gene,
# the cd hit cluster, and whether the gene is single copy in every geneome
# ------------------------------------------------------------------------------

CreateGeneMetaData = function( genePtr )
{
  # Find a representitive member from each cluster
  nGenes    = length( genePtr$clustList )
  clustReps = future_sapply( seq(nGenes), function(i) genePtr$clustList[[i]][1])

  # Look up the description of each gene
  repIdxs = future_sapply(clustReps, function(x) which(genePtr$geneIds == x))

  # Look up the length of the gene
  geneLens = future_sapply(repIdxs, function(i) nchar(genePtr$geneSeqs[i]))

  # Look up the gene description
  geneDescr = future_sapply( seq(nGenes), function(i)
  {
    # Find the genome name of the cluster representitive
    gIdx = which( genePtr$genomeIds == genePtr$genomeIdList[[i]][1] )

    # Look up the row index of the gene in the data-frame of
    # parsed genome features
    rIdx = which( genePtr$gfList[[i]]$featId == clustReps[i] )

    # Return the description of the gene
    return( genePtr$gfList[[i]]$description[ rIdx ] )
  })

  # Create a dataframe with the data on all of selected core genes
  genePtr$geneData = cbind.data.frame( clustReps, geneLens, geneDescr )
}

# ------------------------------------------------------------------------------
