# ------------------------------------------------------------------------------
# CreateGeneData
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function creates a new environment with the the sequences
# of genes used as input to cd-hit. The unique gene ids
# are also stored in a seperate vector
# ------------------------------------------------------------------------------

ConcatenateGeneAlgns = function( genePtr, outDir, runId )
{
  # Make an output directory to store the mafft alignments
  if ( !file.exists(algnDir) ) system( paste("mkdir", algnDir) )

  # Generate the mafft alignments using multi-threading via future.apply
  algnList = future_sapply( 1:length(genePtr$clustList) ,
    function(i) AlgnGeneSeqs( genePtr, i, algnDir )
    )

  # Get the length of the alignments
  algnLens = sapply(1:length(algnList), function(i) nchar( algnList[[i]][1] ) )

  # save(file = "data/2020_02_14_GenerateGeneAlgns_start.Rdata", list = ls())
  # If there is no variation in the final sequence an empty list is returned.
  # Remove any empty genes from the list.
  isEmpty  = algnLens == 0
  algnLens = algnLens[ !isEmpty ]
  algnList = algnList[ !isEmpty ]

  # Remove any genes that had no variation form the data
  genePtr$clustList = genePtr$clustList[ !isEmpty ]
  genePtr$genomeIdList = genePtr$genomeIdList[ !isEmpty ]

  # Generate a vector with the gene start positions in the alignment
  genePtr$genePositions = algnLens[ 1 ]

  # If there is more than one gene, create a vector wit the positions
  # of the gene ends in the alignment
  if ( length(alignList) > 1 )
  {
    for ( i in 2:length( algnLens ) )
    {
      genePtr$genePositions[ i ] =
        genePtr$genePositions[ i - 1 ] + algnLens[ i ]
    }
  }

  # Create the concatenated alignemtent sequence
  concatAlgn = vector("character", length(genePtr$genomeNames) )
  for ( i in 1:length(algnList) )
  {
    if ( length(algnList[[i]]) != length(concatAlgn) )
    {
      stop( "Incorrect number of sequences in alignment ", i, "...\n" )
    }

    if ( !identical( genePtr$genomeNames, names(algnList[[i]]) ) )
    {
      stop( "The alignment for gene ", i, " is in the wrong order...\n" )
    }

    # Append each alignment to the vector of the concatenated
    # alignemnt (passed by reference)
    ConcatenateAlignments( concatAlgn, algnList[[i]] )
  }

  # Write the concatenated gene file
  coGeneFa = paste0(outDir, runId, "concatenated_gene_aa_alignment.fasta" )
  sink( concatGeneFa )
  for ( i in 1:length(concatAlgn) )
  {
    cat('>', genePtr$genomeNames[i], '\n', concatAlgn[i], '\n', sep = '')
  }
  sink()
  
  # Print the statistics on the alignment 
  algnLen = sum( algnLens )
  cat(
    "  -- Number of genes aligned: ", length(algnLens), '\n',
    "  -- The total alignment length is: ", algnLen, '\n',
    sep = ''
    )
  
  return( concatGeneFa )
}

# ------------------------------------------------------------------------------
