# ------------------------------------------------------------------------------
# cogncac: COre GeNe Alignement Concatenation
# 2019/07/15
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function uses cd-hit data to identify core genes to
# construct a concatenated gene alignment. A distance matrix created using
# the amino acid alignments is then used to create a neighbor joining tree.
# ------------------------------------------------------------------------------

cogncac = function(
  fastaFiles,     # Fasta files for the input genomes
  genomeFeatures, # Gff3 or genbank files for the input genomes
  genomeIds,      # Vector of genomes to
  outDir,         # Optional. Directory to write the output files
  runId,          # Optional. Run ID to appent to output files
  minGeneNum,     # Optional. Minimium number of genes to build the tree
  maxMissGenes,   # Optional. Maximium fraction of missing genes
  coreGeneThresh, # Optional. Fraction of genomes with gene to quality as core
  copyNumTresh,   # Optional. Fraction of genomes for a gene to be single copy
  threadVal,      # Optional. Number of threads availible for mafft
  distMat,        # Optional. Bool to create a distance matrix
  njTree,         # Optional. Bool to create a neighbor joining tree
  revTranslate,   # Optional. Bool to convert the aa algiment to dna
  fastTree,       # Optional. Bool to run fastTree
  outGroup,       # Optional. Vector of genomes to exclude for gene selection
  partitionFile,  # Optional. Bool to write a partition file
  aaSubModel,     # Optional. Type of AA subsitiution model for partition file
  dnaSubModel,    # Optional. Type of DNA subsitiution model for partition file
  keepTempFiles   # Optional. Bool to keep mafft and cd-hit files
  percId,         # Optional. Percent ID for the Cd-hit
  algnCovg,       # Optional. Percent alignment coverage for the Cd-hit
  cdHitFlags      # Optional. Parameters to pass to cd-hit to define clusters
  )
{
  startTime = Sys.time() # Start the timer
  # ---- Parse the input arguments ---------------------------------------------

  # If now output directory was specified, write to the working directory
  if ( missing(outDir) )
  {
    outDir = ''
  }
  # Make sure the output directory ends in a '/'
  else if ( !grepl("/$", outDir) && outDir != '' )
  {
    outDir = paste0( outDir, '/' )
  }

  # Make a temporary directory to store any files made during the run
  tempDir = paste0( outDir, "temp_congnac_files/" )
  if ( !file.exists(tempDir) ) system( "mkdir", tempDir )

  # By default, delete any temporary file that are created
  if ( missing(keepTempFiles) ) keepTempFiles = FALSE

  # Default to one availible thread
  if ( missing(threadVal) ) threadVal = as.numeric(future::availableCores())

  # Create a distance matrix with the pairwise distances between isolates
  if ( missing(distMat) ) distMat = FALSE

  # Set the bool ot create the neighbor joining tree
  if ( missing(njTree) )
  {
    njTree = FALSE
  } else if ( njTree ) {
    distMat = TRUE
  }

  # By defaul do not make the RAxML partition file
  if ( missing(partitionFile) ) partitionFile = FALSE

  # ---- Set up multithreadding ------------------------------------------------

  # It's really dumb that you have to do this, but remove the limitation on the
  # maximum allowable object size by future sapplyf
  options( future.globals.maxSize = Inf )

  # Set up multithreadding via future
  plan( future::tweak( future::multiprocess, workers = threadVal )

  # Set the number of threads fot TBB
  RcppParallel::setThreadOptions( numThreads = threadVal )

  # ---- Find the target genes using the gene ids ------------------------------

  # Parse the input files
  cat("\nStep 1: parsing the data on the input genomes\n")
  cat(
    "  -- Wriing results to: ", outDir, "\n",
    "  -- Creating a tree with ", length(gffFiles), " genomes\n",
    sep = ''
    )
  genePtr = CreateGeneDataEnv( genomeFeatures, fastaFiles, genomeIds, tempDir )

  # Identify orthologous genes with cd-hit
  cat("\nStep 2: identifiy orthologues with cd-hit\n")
  FindCogs( genePtr, tempDir, percId, algnCovg, threadVal, cdHitFlags )

  # Use the gene data in the "genePtr" to subset single copyt genes to
  # those that occur as single copy. This function is void, and the clust list
  # within the "genePtr" leaving only the single copy genes
  cat("\nStep 3: Filter for single copy genes\n")
  FilterMultiCopyGenes( genePtr, copyNumTresh )

  # Use the cd-hit results to identify a set of core genes present in all of
  # the input genomes.
  cat("\nStep 5: Selecting genes to create the tree\n")
  SelectCoreGenes( genePtr, minGeneNum, coreGeneFrac, maxMissGenes, outGroup )

  # Make a data frame with the meta data on the genes still in the
  # cd-hit clust list
  cat("\nStep 4: Getting metadata on the selected genes\n")
  CreateGeneData( genePtr )
  save( file = debuggingData, list = ls() )

  # Individually create a new fasta file for each gene and generate the
  # alignment each gene with mafft
  cat("\nStep 6: Generating gene alignments with mafft\n")
  concatGeneFa = ConcatenateGeneAlgns( genePtr, outDir, threadVal )
  save(file = debuggingData, list = ls())

  # Create the environment with the objects to export
  cat( "\nStep 7: Creating output files\n\n" )
  treeEnv = new.env( )

  # If requested, convert the AA alignment to DNA
  if ( revTranslate )
  {
    concatGeneFaNt = ReverseTranslateAlgn( genePtr, concatGeneFa, outDir )
  }

  # Add columns to the data-frame with the meta-data on the genes that
  # correspond to the partitions of the genes in the alignment
  UpdateGeneDataWithPartitions( genePtr, revTranslate )

  # Run fast tree to make the ML tree for the concatenated genes
  if ( fastTree )
  {
    treeEnv$treeFile = paste0( outDir, runId, "_fastTree.tre" )
    fastTreeCmd = paste( "fasttree <", concatGeneFa, ">", treeEnv$treeFile )
    system( fastTreeCmd )
  }

  # If requested, create a distance matrix with the
  if ( distMat )
  {
    treeEnv$distMat = CreateAlgnDistMat( concatGeneSeq )
  }

  # If requested, make a neighbor joining tree with ape
  if ( njTree ) treeEnv$tree = ape::nj( as.dist(treeEnv$distMat) )


  # Add the the data on the individual genes included in the tree
  # to the output
  treeEnv$geneData       = genePtr$geneData
  treeEnv$concatGeneAaFa = concatGeneFa

  # If requested, write a partition file with the positions of each
  # gene in the a lignment
  if ( partitionFile )
    CreatePartionFile( treeEnv, outDir, aaSubModel, dnaSubModel )

  save(file = debuggingData, list = ls())
  # Stop the timer
  end = Sys.time( )
  treeTime = round( difftime(end, start, units = "mins"), 2 )
  cat("\nTree Finished in: ", treeTime, "\n\n")

  return( treeEnv )
}

# ------------------------------------------------------------------------------
