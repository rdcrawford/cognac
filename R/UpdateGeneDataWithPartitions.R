# ------------------------------------------------------------------------------
# Update Gene Data With Partitions
# 2020/04/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Add columns to the gene meta-data with the partitions of the genes in the
# amino acid alignment. The columns are in the RAxML format with the start and
# end postions (ex "1-500"). 
# ------------------------------------------------------------------------------

UpdateGeneDataWithPartitions = function( genePtr, revTranslate )
{
  # Calculate the boundaries of the genes in the amino acid alignment 
  aaGenePartitions      = character( length( genePtr$genePositions ) )
  aaGenePartitions[ 1 ] = paste0( "1-", genePtr$genePositions[ 1 ] )
  for ( i in 2:length( aaGenePartitions ) )                                           
  {
    aaGenePartitions[ i ] = paste0( 
      genePtr$genePositions[ i  + 1 ] + 1, "-", genePtr$genePositions[ i ] 
      )
  }
  
  # Add the aa gene partitions to the dataframe with the data on the
  # genes
  genePtr$geneData = cbind.data.frame( 
    genePtr$geneData, aaGenePartitions, stringsAsFactors = FALSE
    )
  
  if ( revTranslate )
  {
    # Get the positions of the genes in the nucluetide alignment. This is 
    # the number of amino acids + 1 for the stop codon, and multiplied by
    # 3 for nucleotides per amino acid
    ntGenePositions = sapply( genePtr$genePositions,
      function(x) ( x + 1 ) * 3
      )
    
    # Calculate the boundaries of the genes
    ntGenePartitions      = character( length( ntGenePositions ) )
    ntGenePartitions[ 1 ] = paste0( "1-", ntGenePositions[ 1 ] )
    for ( i in 2:length( ntGenePartitions ) )                                           
    {
      # Find the partitions of the genes in the nt alignment
      ntGenePartitions[ i ] = paste0( 
        ntGenePositions[ i  + 1 ] + 1, "-", ntGenePositions[ i ] 
        )
      
      # Add the nt gene partitions to the dataframe with the data on the
      # genes
      genePtr$geneData = cbind.data.frame( 
        genePtr$geneData, ntGenePositions, stringsAsFactors = FALSE
        )
    }
  }
}

# ------------------------------------------------------------------------------