# ------------------------------------------------------------------------------
# Sort Matrix
# 2020/04/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function takes a matrix and another matrix to use as a reference. The 
# rows and columns of the first matrix are sorted by the row names and 
# column names of the second matrix. The sorted matrix is returned.
# ------------------------------------------------------------------------------

CreatePartionFile = function( genePtr, outDir, runId, aaSubModel, dnaSubModel )
{
  # Set the default substitution model if none was provided
  if ( missing(aaSubModel) )
  {
    aaSubModel = "LG"
  }
  if ( missing(dnaSubModel) )
  {
    dnaSubModel = "DNA"
  }
  
  CLUST_REPS     = 1
  GENE_LENS      = 2
  GENE_DESC      = 3
  AA_PARTITIONS  = 4
  DNA_PARTITIONS = 5
  
  nGenes = nrow( genePtr$geneData )
  # Write the DNA partition file 
  sink( paste(outDir, runId, "aa_gene_partitions.txt" ) )
  for ( i in seq( nGenes ) )
  {
    cat(
      aaSubModel,
      ", ",
      genePtr$geneData[ i, CLUST_REPS ],
      " = ",
      genePtr$geneData[ i, AA_PARTITIONS ],
      '\n',
      sep = ''
      )
  }
  sink( )
  
  if ( ncol( genePtr$geneData ) == DNA_PARTITIONS )
  {
    sink( paste(outDir, runId, "nt_gene_partitions.txt" ) )
    for ( i in seq( nGenes ) )
    {
      cat(
        dnaSubModel,
        ", ",
        genePtr$geneData[ i, CLUST_REPS ],
        " = ",
        genePtr$geneData[ i, DNA_PARTITIONS ],
        '\n',
        sep = ''
        )
    }
    sink( )
  }
  
}

# ------------------------------------------------------------------------------