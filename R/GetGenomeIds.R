#  -----------------------------------------------------------------------------
# GetGenomeIds
# 2020/05/21
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function provides error checking on the unique genome identifiers. If
# no genome IDs were supplied, they are created by removing the path and 
# file extension. Genome IDs are then checked for duplicates. 
# ------------------------------------------------------------------------------

GetGenomeIds = function( 
  featureFiles, fastaFiles, fastaExt, featureExt, genomeIds
  )
{
  # Add the names of the genomes included in the analysis  to the environmnet.
  # if they are missing, extract the name from the paths
  if ( missing(genomeIds) )
  {
    # If no extension was supplied, get the extension as the last
    # '.' in the file name
    if ( missing(fastaExt) )
    {
      # Make the genome names by extracting the last element of the path
      # and removing the file extension.
      geneEnv$genomeNames = sapply( 
        fastaFiles, ExtractGenomeNameFromPath, USE.NAMES = FALSE
        )
    
    # Use the file extension to get the name
    } else {
      
      # Make the genome names by extracting the last element of the path
      # and removing the file extension.
      geneEnv$genomeNames = sapply(
        fastaFiles, GetGenomeNameWithExt, fastaExt, USE.NAMES = FALSE
        )

    }
    
    # If no extension was supplied, get the extension as the last
    # '.' in the file name
    if ( missing(featureExt) )
    {
      # Make the genome names by extracting the last element of the path
      # and removing the file extension.
      genomeIds = sapply(
        featureFiles, ExtractGenomeNameFromPath, USE.NAMES = FALSE
        )
      
    # Use the file extension to get the name  
    } else {
      
      # Make the genome names by extracting the last element of the path
      # and removing the file extension.
      genomeIds = sapply( 
        featureFiles, GetGenomeNameWithExt, featureExt, USE.NAMES = FALSE
        )
    }
    
    #If the genome IDs made by using the fasta files and feature files 
    # are not identical, throw warnings
    if ( !identical( geneEnv$genomeNames, genomeIds ) )
    {
      isInFaNames = genomeIds %in% geneEnv$genomeNames
      
      if ( FALSE %in% isInFaNames )
      {
        i = which( !isInFaNames )[1]
        warning(
          "Creating the genome Ids from the fasta files and genome ",
          "feature files produced different results in ",
          sum(isInFaNames), " out of ", length(fastaFiles), 
          " instances\n\n", "For example:\n",
          " -- Fasta file name: ", geneEnv$genomeNames[i], '\n',
          " -- Feature file name: ", genomeIds[i], '\n'
          )
        
      } else {
        
        warning(
          "All genome IDs are present when creating genome IDs from ",
          "fasta files and genome feature files, but they are in a ",
          "different order. The files may not be sorted properly..."
          )
      }
    }
  }

  # Check that all of the genome IDs are unique
  isDuplicated = duplicated( geneEnv$genomeNames )
  if ( TRUE %in% isDuplicated )
  {
    numDuplicates  = sum( isDuplicated )
    duplicatedIdxs = which( isDuplicated )
    MAX_CASES = 5
    # Create the error string stating the duplicated genomes
    if ( numDuplicates > MAX_CASES )
    {
      errStr = paste( 
        c(paste( "  --", geneEnv$genomeNames[isDuplicated][1:MAX_CASES]), 
        "     ...\n"),
        collapse = "\n" 
        )

    } else {
      
      errStr = paste( 
        paste( "  --", geneEnv$genomeNames[isDuplicated]), 
        collapse = "\n" 
        )
    
    }
    
    stop(
      "There are duplicated genome Ids in ", numDuplicates, " cases:\n",
      errStr,
      "\nGenome IDs must be unique\n"
      )
  }
  
  return( genomeIds )
}

# ------------------------------------------------------------------------------