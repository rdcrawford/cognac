# ------------------------------------------------------------------------------
# Get CD-HIT Clusters
# 2018/11/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This script concatenates the given annotated gene files and runs cd-hit
# automatically selecting the appropriate word size for the desired percent
# id threshold for clustering. The CD-Hit results are returned as a list
# ------------------------------------------------------------------------------

FindCogs = function(
  geneEnv,     # Environment with the data on the genes
  outDir,      # Directory to write Cd-Hit output files
  percId,      # Percent ID clustering threashold
  algnCovg,    # Percent coverage for the alignment
  threadVal,   # Number of threads to used
  cdHitFlags,  # Any additional flags to pass to Cd-Hit
  maxMissGenes # Maximium allowable missing genes to be included
  )
{
  # ---- Parse the input arguments ---------------------------------------------

  # Assign missing arguemnts to the default values
  if ( missing( percId ) )       percId     = 0.7
  if ( missing( algnCovg ) )     algnCovg   = 0.80
  if ( missing( outDir ) )       outDir     = paste0(getwd(), '/')
  if ( missing( cdHitFlags ) )   cdHitFlags = "-M 0 -d 0 -g 1"
  if ( missing( threadVal ) )    threadVal  = 1
  if ( missing( maxMissGenes ) )
  {
    minGeneNum = 2
  } else {
    minGeneNum = floor(
      length( geneEnv$genomeNames ) * ( 1 - maxMissGenes )  * 0.85
      )
  }
  
  # The working directory has to end in a forward slash. If it doesn't,
  # add one
  if (!grepl("/$", outDir)) outDir = paste0(outDir, '/')

  # ---- Variable Initializations ----------------------------------------------

  # Create the stings for the cd hit input and output file names
  geneAnnotExt       = ".faa"
  cdHitGenesFileName = paste0( outDir, "cdHitClusters", geneAnnotExt )
  cdHitClstrFileName = paste0( outDir, "cdHitClusters", geneAnnotExt, ".clstr" )
  cdHitLogFile       = paste0( outDir, "cdHit.log" )

  # ---- Select the word size --------------------------------------------------

  if ( percId >= 0.7 ) {
    wordSize = 5
  } else if ( percId >= 0.6 ) {
    wordSize = 4
  } else if ( percId >= 0.5 ) {
    wordSize = 3
  } else {
    wordSize = 2
  }

  # ---- Run CD-Hit ------------------------------------------------------------

  cdHitCmd = paste(
    "cd-hit",
    "-i",  geneEnv$faaPath,    # Input file with the translated aa seqs
    "-o",  cdHitGenesFileName, # Ouput
    "-c",  percId,             # Min percent identity in the gene seqs
    "-aL", algnCovg,           # Minimum disparity in the length of the seqs
    "-T",  threadVal,          # Number of threads
    "-n",  wordSize,           # Size word to fraction sequences into
    trimws( cdHitFlags ),      # Additional flags for the cd-hit run
    ">",   cdHitLogFile        # Write the output to a temp log file
    )
  system( cdHitCmd )
  
  # ---- Parse the CD-Hit Data -------------------------------------------------

  # Assign the clust list and gene matrix to varibles in the
  # environment containing the gene data
  ParseCdHit( cdHitClstrFileName, FALSE, minGeneNum, geneEnv )
  
  # ---- Check if there are genomes that dont share genmes ---------------------
  
  # Find the number of genes present in each genome
  nGenesPerGenome = sapply( 1:nrow( geneEnv$geneMat ),
    function(i) sum( geneEnv$geneMat[ i, ] != 0 )
    )
  
  # Find if there are any genomes with 0 genes
  hasNoGenes =  nGenesPerGenome == 0
  
  # If there are genomes with no genes, remove them from the analysis. 
  # Most likely, these genomes are of poor quality and will have a negative 
  # impact on the down stream analysis.
  if ( TRUE %in% hasNoGenes )
  {
    # Remove any trace that these genomes ever existed
    RemoveGenomesFromAnalisys( geneEnv, geneEnv$genomeNames[ hasNoGenes ] )
    
    nBadGenomes = sum( hasNoGenes )
    cat(
      "  -- There were ", nBadGenomes, 
      " genomes with no genes shared by at least ", minGeneNum, 
      " genomes in this dataset\n",
      "  -- Continuing with the ", length( geneEnv$genomeNames ), 
      " remaining genomes\n",
      sep = ''
      )
  }
  
  # Print out the relevant stats from the cd-hit run
  cat(
    "  -- The genes were classified into ", geneEnv$nCogs, 
    " clusters of orthologous genes\n", 
    "  -- Continuing with ", length( geneEnv$clustList ), " genes present",
    " in at least ", minGeneNum, " occurrences\n",
    sep = ''
    )
}

# ------------------------------------------------------------------------------
