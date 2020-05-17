# ------------------------------------------------------------------------------
# Create Gene Data Env
# 2020/02/05
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function initializes the environment containing data on the genes
# used in the analysis. This function adds the list of parsed gff files,
# fasta file paths, and the genome Ids of the genomes included in the analysis
# ------------------------------------------------------------------------------

CreateGeneDataEnv = function( featureFiles, fastaFiles, genomeIds, outDir )
{
  # Currently fasta files are required. Throw an error if the fasta
  # files are not included
  if ( missing(fastaFiles) )
    stop("Missing required argument with the paths to the fasta files\n")

  # Initalize an environment to store the data on the genes. This will
  # be passed to further functions and iteratively be subset to only
  # include the genes contained in the alignment
  genePtr = new.env()

  # Add the names of the genomes included in the analysis  to the environmnet.
  # if they are missing, extract the name from the paths
  if ( missing(genomeIds) )
  {
    # Make the genome names by extracting the last element of the path
    # and removing the file extension.
    genePtr$genomeNames = sapply( fastaFiles, ExtractGenomeNameFromPath )
    
    # Get counts of the number of input vectors 
    argVecLen = c( length(featureFiles), length(fastaFiles) )

  } else {

    # Add the genome ids to the envriroment that will be returned
    genePtr$genomeNames = genomeIds
    
    # Get counts of the number of input vectors 
    argVecLen =
      c( length(featureFiles), length(fastaFiles), length(genomeIds) )
  }

  # Check that everything has the same length
  if ( length( unique(argVecLen) ) != 1 )
  {
    argNames = c("featureFiles", "fastaFiles", "genomeIds")
    vecLenStr = sapply( 1:length(argVecLen),
      function(i) 
        paste0( "  -- ", argNames[i], ": ", argVecLen[i], " elements\n")
      )
    inputEfforMessage = paste(
      c("\nThe input data vectors are not the same length:\n", vecLenStr), 
      collapse  = ''
      )
    stop( inputEfforMessage )
  }

  # Check that all of the genome IDs are unique
  isDuplicated = duplicated( genePtr$genomeNames )
  if ( TRUE %in% isDuplicated )
  {
    stop(
      "There are duplicated genome ids: \n",
      paste( genePtr$genomeNames[isDuplicated], collapse = ' ' ),
      '\n'
      )
  }

  # Create the path to the cd-hit input file. This is where the amino acid
  # sequences for all of the genome will be written
  faaPath = paste0( outDir, "allGenes.faa" )

  # Parse the data on the input genomes
  CreateCognacRunData( genePtr, featureFiles, fastaFiles, faaPath )

  # Add the fasta files and path to the parsed genome data
  genePtr$fastaFiles = fastaFiles
  genePtr$faaPath    = faaPath
  
  return( genePtr )
}

# ------------------------------------------------------------------------------
