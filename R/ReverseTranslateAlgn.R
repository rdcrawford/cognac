# ------------------------------------------------------------------------------
# ReverseTranslateAlgn
# 2020//13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function inputs a concatenated gene alignment
# ------------------------------------------------------------------------------

ReverseTranslateAlgn = function( genePtr, concatGeneFa, outDir )
{
  # Define the name of the output file
  concatGeneDnaFa = paste0( outDir, "concat_align_nt_seq.fasta" )
  if ( file.exists( concatGeneDnaFa ) ) system( paste("rm", concatGeneDnaFa) )

  # Read in the concatenated gene alignment
  concatGeneSeq = ParseFasta( concatGeneFa )

  # For each observation in the concatenated alignment, look up the sequences
  # in the gff file
  for ( i in 1:length(concatGeneSeq) )
  {
    # Look up the row in the gff file corresponging to each core gene
    gfRowIdxs = sapply( 1:length(genePtr$clustList), function(j)
    {
      # Find position the current genome in the vector of genome names
      isThisGenome = genePtr$genomeIdList[[j]] == genePtr$genomeNames[i]
      if ( !TRUE %in% isThisGenome ) return( NA )

      # Find the gene id and look up the row in the gff file corresponding
      # to this gene
      listIdx = which( isThisGenome )[ 1 ]
      geneId  = genePtr$clustList[[ j ]][ listIdx ]
      return( which( genePtr$gfList[[ i ]][ , GENE_ID ] == geneId ) )
    })

    # If there are any missing genes core genes represented as na in the
    # vecotr, remove them
    if ( TRUE %in% is.na(gfRowIdxs) ) gfRowIdxs = gfRowIdxs[ !is.na(gfRowIdxs) ]

    # Reverse translate the current sequence in the concatenated gene alignemnt
    TranslateAaAlgnToDna(
      genePtr$gfList[[i]][ gfRowIdxs, ],
      genePtr$fastaFiles[i],
      genePtr$genePositions,
      genePtr$genomeNames[i],
      concatGeneSeq[i],
      concatGeneDnaFa
      )
  }

  return( concatGeneDnaFa )
}

# ------------------------------------------------------------------------------
