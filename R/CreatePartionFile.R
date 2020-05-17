# ------------------------------------------------------------------------------
# CreatePartionFile
# 2020/04/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function creates the gene partition file used by RAxML
# ------------------------------------------------------------------------------

CreatePartionFile = function( genePtr, outDir, runId, subModel, seqType )
{
  CLUST_REPS     = 1
  GENE_LENS      = 2
  GENE_DESC      = 3
  AA_PARTITIONS  = 4
  DNA_PARTITIONS = 5
  
  # Set the default substitution model if none was provided
  if ( missing(seqType) )
  {
    colIdx = ncol( genePtr$geneData )
  } else if ( seqType == "DNA" ) {
    colIdx = DNA_PARTITIONS
  } else {
    colIdx = AA_PARTITIONS
  }
  
  # Set the substitution model to be used 
  if ( missing(subModel)  && colIdx == DNA_PARTITIONS )
  {
    subModel = "DNA"
  } else if ( missing(subModel)  && colIdx == AA_PARTITIONS ) {
    subModel = "LG"
  }
  
  # Write the partition file 
  sink( paste( outDir, runId, "concat_algn_gene_partitions.txt" ) )
  for ( i in seq( nGenes ) )
  {
    cat(
      subModel,
      ", ",
      genePtr$geneData[ i, CLUST_REPS ],
      " = ",
      genePtr$geneData[ i, colIdx ],
      '\n',
      sep = ''
      )
  }
  sink( )
}

# ------------------------------------------------------------------------------