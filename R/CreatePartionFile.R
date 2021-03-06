# ------------------------------------------------------------------------------
# CreatePartionFile
# 2020/04/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#' Create Partion File
#' @description 
#' This is a function creates the gene partition file used by RAxML
#' @details 
#' This is a function that writes the alignmnet partition file used by RAxML
#' to apply. Written in the format:
#' DNA, peg.1=1-30
#' DNA, peg.2=31-60
#' @param algnEnv Environment created by the cognac that contains the 
#'   alignment data.
#' @param outDir Directory to write the file to
#' @param runId Run identifier appended to the output file
#' @param subModel substitution model to be used
#' @return void
#' @export
# ------------------------------------------------------------------------------

CreatePartionFile = function( algnEnv, outDir, runId, subModel, seqType )
{
  # Constant declarations: column indicies in the gene metadata
  CLUST_REPS     = 1
  GENE_LENS      = 2
  GENE_DESC      = 3
  AA_PARTITIONS  = 4
  DNA_PARTITIONS = 5
  
  # Set the default substitution model if none was provided
  if ( missing( seqType ) && ncol( algnEnv$geneData ) )
  {
    colIdx = ncol( algnEnv$geneData )
  }
  else if ( seqType == "nt" )
  {
    colIdx = DNA_PARTITIONS
  } 
  else 
  {
    colIdx = AA_PARTITIONS
  }
  
  # Set the substitution model to be used 
  if ( missing( subModel )  && colIdx == DNA_PARTITIONS )
  {
    subModel = "DNA"
  } else if ( missing( subModel )  && colIdx == AA_PARTITIONS ) {
    subModel = "LG"
  }
  
  # Write the partition file 
  sink( paste( outDir, runId, "concat_algn_gene_partitions.txt" ) )
  for ( i in seq( nGenes ) )
  {
    cat(
      subModel,
      ", ",
      algnEnv$geneData[ i, CLUST_REPS ],
      " = ",
      algnEnv$geneData[ i, colIdx ],
      '\n',
      sep = ''
      )
  }
  sink( )
}

# ------------------------------------------------------------------------------