#  -----------------------------------------------------------------------------
#  ReverseTranslateAlgn
#  2020/05/13
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#' Reverse Translate Alignment
#' @description
#' This function reads in a concatenated gene alignment created by cognac and
#' @param geneEnv, concatGeneFa, outDir
#' @param concatGeneFa Path to the amino acid alignment
#' @param runId Run identifier to append to the alignment file
#' @param algnEnv Environment created by cognac with the data on the alignment
#' @return Path to the reverse translated alignment
#' @export
#  -----------------------------------------------------------------------------

ReverseTranslateAlgn = function( 
  geneEnv, concatGeneFa, outDir, runId, algnEnv
  )
{
  # Constant declarations
  AA_PARTITIONS  = 4 # Column index of the aa gene partitions 
  
  # Get indexes in the alignment where the genome partitions are
  if ( missing(algnEnv) )
  {
    # Assign the indexes of the gene positions
    genePartitions = geneEnv$genePositions
    
  } else {
    
    # Extract the gene positions from the meta-data
    genePartitions = sapply( algnEnv$geneData[ , AA_PARTITIONS ],
      function(x) as.numeric( strsplit(x, ':')[[1]][2] )
      )
  }
  
  # Define the name of the output file
  concatGeneDnaFa = 
    paste0( outDir, runId,  "concatenated_gene_nt_alignment.fasta" )
  if ( file.exists( concatGeneDnaFa ) ) system( paste("rm", concatGeneDnaFa) )

  # Read in the concatenated gene alignment
  concatGeneSeq = ParseFasta( concatGeneFa )

  # For each observation in the concatenated alignment, look up the sequences
  # in the gff file
  for ( i in 1:length(concatGeneSeq) )
  {
    # Look up the row in the gff file corresponging to each core gene
    gfRowIdxs = sapply( 1:length(geneEnv$clustList), function(j)
    {
      # Find position the current genome in the vector of genome names
      isThisGenome = geneEnv$genomeIdList[[j]] == geneEnv$genomeNames[i]
      if ( !TRUE %in% isThisGenome ) return( NA )

      # Find the gene id and look up the row in the gff file corresponding
      # to this gene
      listIdx = which( isThisGenome )[ 1 ]
      geneId  = geneEnv$clustList[[ j ]][ listIdx ]
      return( which( geneEnv$gfList[[ i ]]$featId == geneId ) )
    })
    
    # If there are any missing genes core genes represented as na in the
    # vecotr, remove them
    isMissingGene = is.na(gfRowIdxs)
    if ( TRUE %in% isMissingGene ) gfRowIdxs = gfRowIdxs[ !isMissingGene ]
    
    # Reverse translate the current sequence in the concatenated gene alignemnt
    TranslateAaAlgnToDna(
      geneEnv$gfList[[i]][ gfRowIdxs, ],
      geneEnv$fastaFiles[i],
      genePartitions,
      geneEnv$genomeNames[i],
      concatGeneSeq[i],
      concatGeneDnaFa
      )
  }

  return( concatGeneDnaFa )
}

# ------------------------------------------------------------------------------
