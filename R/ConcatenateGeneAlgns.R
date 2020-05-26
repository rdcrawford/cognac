# ------------------------------------------------------------------------------
# CreateGeneData
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function creates a new environment with the the sequences
# of genes used as input to cd-hit. The unique gene ids
# are also stored in a seperate vector
# ------------------------------------------------------------------------------

ConcatenateGeneAlgns = function( geneEnv, outDir, runId )
{
  # Make an output directory to store the mafft alignments
  algnDir = paste0( outDir, "temp_congnac_files", runId, "/mafft_alignments/" )
  if ( !file.exists(algnDir) ) system( paste("mkdir", algnDir) )

  # Generate the mafft alignments using multi-threading via future.apply
  algnList = future.apply::future_sapply( 1:length(geneEnv$clustList) ,
    function(i) AlgnGeneSeqs( geneEnv, i, algnDir )
    )
  
  # Get the length of the alignments
  algnLens = sapply(1:length(algnList), function(i) nchar( algnList[[i]][1] ) )

  # If there is no variation in the final sequence an empty list is returned.
  # Remove any empty genes from the list.
  isEmpty  = algnLens == 0
  algnLens = algnLens[ !isEmpty ]
  algnList = algnList[ !isEmpty ]

  # Remove any genes that had no variation form the data
  geneEnv$clustList = geneEnv$clustList[ !isEmpty ]
  geneEnv$genomeIdList = geneEnv$genomeIdList[ !isEmpty ]

  # Generate a vector with the gene start positions in the alignment
  geneEnv$genePositions = algnLens[ 1 ]

  # If there is more than one gene, create a vector wit the positions
  # of the gene ends in the alignment
  if ( length(algnList) > 1 )
  {
    for ( i in 2:length( algnLens ) )
    {
      geneEnv$genePositions[ i ] =
        geneEnv$genePositions[ i - 1 ] + algnLens[ i ]
    }
  }

  # Create the concatenated alignemtent sequence
  concatAlgn = vector("character", length(geneEnv$genomeNames) )
  for ( i in 1:length(algnList) )
  {
    if ( length(algnList[[i]]) != length(concatAlgn) )
    {
      stop( "Incorrect number of sequences in alignment ", i, "...\n" )
    }

    if ( !identical( geneEnv$genomeNames, names(algnList[[i]]) ) )
    {
      stop( "The alignment for gene ", i, " is in the wrong order...\n" )
    }

    # Append each alignment to the vector of the concatenated
    # alignemnt (passed by reference)
    ConcatenateAlignments( concatAlgn, algnList[[i]] )
  }

  # Write the concatenated gene file
  coGeneFa = paste0(outDir, runId, "concatenated_gene_aa_alignment.fasta" )
  sink( coGeneFa )
  for ( i in 1:length(concatAlgn) )
  {
    cat('>', geneEnv$genomeNames[i], '\n', concatAlgn[i], '\n', sep = '')
  }
  sink()
  
  # Print the statistics on the alignment 
  algnLen = sum( algnLens )
  cat(
    "  -- Number of genes aligned: ", length(algnLens), '\n',
    "  -- The total alignment length is: ", algnLen, '\n',
    sep = ''
    )
  
  return( coGeneFa )
}

# ------------------------------------------------------------------------------
