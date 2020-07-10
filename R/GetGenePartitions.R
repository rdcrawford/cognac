# ------------------------------------------------------------------------------
# Update Gene Data With Partitions
# 2020/04/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Add columns to the gene meta-data with the partitions of the genes in the
# amino acid alignment. The columns are in the RAxML format with the start and
# end postions (ex "1-500"). 
# ------------------------------------------------------------------------------

GetGenePartitions = function( genePtr, revTranslate )
{
  # Calculate the boundaries of the genes in the amino acid alignment 
  aaGenePartitions      = character( length( genePtr$genePositions ) )
  aaGenePartitions[ 1 ] = paste0( "1-", genePtr$genePositions[ 1 ] )
  for ( i in 2:length( aaGenePartitions ) )                                           
  {
    aaGenePartitions[ i ] = paste0(
      genePtr$genePositions[ i  - 1 ] + 1, "-", genePtr$genePositions[ i ] 
      )
  }
  
  # If this alignmnt was not reverse translated, return the amino acid 
  # gene partitions.
  if ( !revTranslate ) return( aaGenePartitions )

  # Get the positions of the genes in the nucluetide alignment. This is 
  # the number of amino acids + 1 for the stop codon, and multiplied by
  # 3 for nucleotides per amino acid
  ntGeneEnds = genePtr$genePositions * 3
  
  # Calculate the boundaries of the genes
  ntGenePartitions      = character( length( ntGeneEnds ) )
  ntGenePartitions[ 1 ] = paste0( "1-", ntGeneEnds[ 1 ] )
  for ( i in 2:length( ntGenePartitions ) )                                           
  {
    # Find the partitions of the genes in the nt alignment
    ntGenePartitions[ i ] = paste0( 
      ntGeneEnds[ i  - 1 ] + 1, "-", ntGeneEnds[ i ]
      )
  }

  # Return the gene positions asa character matrix
  return( cbind( aaGenePartitions, ntGenePartitions ) )
}

# ------------------------------------------------------------------------------