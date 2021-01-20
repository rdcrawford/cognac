# ------------------------------------------------------------------------------
# Select Core Genes
# 2019/12/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Remove perfectly conserved genes without any variation from the data set and
# remove any genomes that do not have a complete set of core genes
# ------------------------------------------------------------------------------

SelectAlgnGenes = function(
  geneEnv,        # Environment with the data on the genes
  minGeneNum,     # Minimum number of genes to include in the alignment
  coreGeneThresh, # Fraction of genome that a gene must be present in
  maxMissGenes,   # Maximium fraction of genes that can be missing
  outGroup        # Names of out group genomes not used to select core genes
  )
{
  # If missing the fraction of genomes a given gene must be in, set the
  # default to 99%
  if ( missing(coreGeneThresh) ) coreGeneThresh = 0.99

  # Set the minimium number of genes to build the tree
  if ( missing(minGeneNum) ) minGeneNum = 2

  if ( missing(maxMissGenes) ) maxMissGenes = 1
  
  # Save the number of genes and genomes before filtering
  numClusts  = length( geneEnv$clustList )
  numGenomes = length( geneEnv$genomeNames )

  # To select core genes of interest iterate over the data until all
  # of the remaining genes have some degree of variation and
  isKeeper = rep( TRUE, numGenomes )

  # If there is an outgroup selected, remove this genome for the
  # core gene selection
  if ( !missing(outGroup) )
  {
    isOutGroup = geneEnv$genomeNames %in% outGroup
    for ( i in which( isOutGroup ) ) isKeeper[ i ] = FALSE
  }

  repeat
  {
    # Check that there are suffienct gene counts in the data set to make the
    # tree
    if ( is.null(dim(geneEnv$geneMat)) )
      stop( "There are not enough conserved genes to make the tree..." )
    if ( ncol(geneEnv$geneMat) < minGeneNum )
      stop( "There are not enough conserved genes to make the tree..." )

    # Remove any outlier genomes that dont have core genes
    isKeeper = RemoveOutlierGenomes(
      geneEnv$geneMat, minGeneNum, coreGeneThresh, isKeeper
      )

    # Count the number of core genes in the dataset
    isCoreGene = CalcNumCoreGenes( geneEnv, coreGeneThresh, isKeeper )
    coreGeneCount = sum( isCoreGene )

    if ( length( isCoreGene ) == 0 || length( isCoreGene ) < minGeneNum )
      stop( paste( "Unable to find", minGeneNum, "core genes in these data" ) )
    
    # If there is variation in all of the genes after removing outliers
    # that do not share the core set of genes,
    if ( coreGeneCount >= minGeneNum ) break
  }

  cat("  -- Identified ",  sum( isCoreGene ), " core genes\n" )
  
  # Check and see that at least one remaining genome has variation
  # in each gene. If genes are conserved in all of the remaining genomes
  # then there is no point in aligning them
  isNotConserved = sapply(1:ncol(geneEnv$geneMat),
    function(j) FALSE %in% (geneEnv$geneMat[ isKeeper, j ] == 100)
    )
  
  # Subset to only include the core genes
  geneEnv$clustList    = geneEnv$clustList[ isCoreGene ]
  geneEnv$genomeIdList = geneEnv$genomeIdList[ isCoreGene ]
  geneEnv$geneMat      = geneEnv$geneMat[ , isCoreGene ]
  
  # If there are genes without variation, remove them from the dataset
  if ( FALSE %in% isNotConserved )
  {
    geneEnv$geneMat      = geneEnv$geneMat[ , isNotConserved ]
    geneEnv$clustList    = geneEnv$clustList[ isNotConserved ]
    geneEnv$genomeIdList = geneEnv$genomeIdList[ isNotConserved ]
    
    cat( 
      "  -- Removing ", sum( !isNotConserved ), " perfectly conserved",
      "genes from the analysis\n", sep = ''
      )
  }
  
  # Calcuate the fraction of missing genes
  missingGeneFrac = sapply( 1:nrow(geneEnv$geneMat),
    function(i) sum( geneEnv$geneMat[ i, ] == 0 ) / coreGeneCount
    )
  
  if ( 1 %in% missingGeneFrac )
  {
    isMissingAll =  missingGeneFrac == 1
    
    if ( !missing(outGroup) )
    {
      if ( outGroup %in% geneEnv$genomeNames[ isMissingAll ] )
      {
        stop("The supplied outgroup (", outGroup, ") has none of the ",
          "selected core genes. Consider altering the cd-hit parameters ",
          " to be less stringent."
          )
      }
    }
    
    # Update the genomes and genes that are still being used for the analysis
    RemoveGenomesFromAnalisys( geneEnv, geneEnv$genomeNamest[ isMissingAll ] )
    cat(
      "  -- Removing ", sum(isMissingAll), " genomes missing all of the ",
      "selected core genes\n",
      sep = ''
      )
  }
  
  # If there is a threshold on the number of genes which can be missing
  if ( maxMissGenes < 1 )
  {
    # Calcuate the fraction of missing genes
    missingGeneFrac = sapply( 1:nrow(geneEnv$geneMat),
      function(i) sum( geneEnv$geneMat[ i, ] == 0 ) / coreGeneCount
      )
    
    # Remove genomes that are missing too many genes
    isMissingTooMuch = missingGeneFrac > maxMissGenes

    # If there is an outgroup, do not remove it from the analysis
    if ( !missing(outGroup) )
    {
      for ( i in which( isOutGroup ) ) isMissingTooMuch[ i ] = FALSE
    }

    if ( TRUE %in% isMissingTooMuch )
    {
      numStillIn = sum( !isMissingTooMuch )

      # Make sure there is an alignment for atleast two genomes
      if ( numStillIn < 2 )
      {
        stop(
          "There are no genomes with less than",
          maxMissGenes * 100, "% of genes missing\n"
          )
      }

      cat(
        "  -- ",  sum(isMissingTooMuch), " genomes are missing greater than ",
        maxMissGenes * 100, "% of genes\n", "  -- There are ",  numStillIn,
        " genomes remaining\n",
        sep = ''
        )

      # Remove the genomes missing too many genes from the analysis
      RemoveGenomesFromAnalisys( 
        geneEnv, geneEnv$genomeNames[ isMissingTooMuch ] 
        )
    }
  }

  cat(
    "  -- ",  length( geneEnv$clustList ), 
    " genes met the criteria to be included in the alignment\n",
    sep = ''
    )
  
  # Remove any non-core genes from the data
  coreGenes        = unlist( geneEnv$clustList )
  isCoreGene       = geneEnv$geneIds %in% coreGenes
  geneEnv$geneSeqs = geneEnv$geneSeqs[ isCoreGene ]
  geneEnv$geneIds  = geneEnv$geneIds[ isCoreGene ]
}

# ------------------------------------------------------------------------------
