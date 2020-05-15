# ------------------------------------------------------------------------------
# Select Core Genes
# 2019/12/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Remove perfectly conserved genes without any variation from the data set and
# remove any genomes that do not have a complete set of core genes
# ------------------------------------------------------------------------------

SelectCoreGenes = function(
  genePtr,        # Environment with the data on the genes
  minGeneNum,     # Minimum number of genes to include in the alignment
  coreGeneThresh, # Fraction of genome that a gene must be present in
  maxMissGenes,   # Maximium fraction of genes that can be missing
  outGroup        # Names of out group genomes not used to select core genes
  )
{
  # If missing the fraction of genomes a given gene must be in, set the
  # default to 99%
  if ( missing(coreGeneThresh) ) coreGeneThresh = 0.99

  # Create the path to the tree file
  if ( missing(fastTree) ) fastTree = FALSE

  # Set the minimium number of genes to build the tree
  if ( missing(minGeneNum) ) minGeneNum = 2

  # Save the number of genes and genomes before filtering
  numClusts  = length( genePtr$clustList )
  numGenomes = length( genePtr$genomeNames )

  # To select core genes of interest iterate over the data until all
  # of the remaining genes have some degree of variation and
  isKeeper = rep( TRUE, numGenomes )

  # If there is an outgroup selected, remove this genome for the
  # core gene selection
  if ( !missing(outGroup) )
  {
    isOutGroup = genePtr$genomeNames %in% outGroup
    for ( i in which( isOutGroup ) ) isKeeper[ i ] = FALSE
  }

  repeat
  {
    # Check that there are suffienct gene counts in the data set to make the
    # tree
    if ( is.null(dim(genePtr$geneMat)) )
      stop( "There are not enough conserved genes to make the tree..." )
    if ( ncol(genePtr$geneMat) < minGeneNum )
      stop( "There are not enough conserved genes to make the tree..." )

    # Remove any outlier genomes that dont have core genes
    isKeeper = RemoveOutlierGenomes(
      genePtr$geneMat, minGeneNum, coreGeneThresh, isKeeper
      )

    # Check and see that at least one remaining genome has variation
    # in each gene. If genes are conserved in all of the remaining genomes
    # then there is no point in aligning them
    isNotConserved = sapply(1:ncol(genePtr$geneMat),
      function(j) FALSE %in% (genePtr$geneMat[ isKeeper, j ] == 100)
      )

    # If there are genes without variation, remove them from the dataset
    if ( FALSE %in% isNotConserved )
    {
      genePtr$geneMat      = genePtr$geneMat[ , isNotConserved ]
      genePtr$clustList    = genePtr$clustList[ isNotConserved ]
      genePtr$genomeIdList = genePtr$genomeIdList[ isNotConserved]
    }

    # Count the number of core genes in the dataset
    isCoreGene = CalcNumCoreGenes( genePtr, coreGeneThresh )
    coreGeneCount = sum( isCoreGene )

    # If there is variation in all of the genes after removing outliers
    # that do not share the core set of genes,
    if ( coreGeneCount >= minGeneNum ) break
  }

  if ( !missing(maxMissGenes) )
  {
    if ( maxMissGenes == 1 ) break
    cat(
      "  -- Removing genomes with more than ", maxMissGenes * 100,
      "% missing genes\n",
      sep = ''
      )
    # Calcuate the fraction of missing genes
    missingGeneFrac = sapply( 1:nrow(genePtr$geneMat),
      function(i) sum( genePtr$geneMat[i, isCoreGene] == 0 ) / coreGeneCount
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
      # Update the genomes and genes that are still being used for the analysis
      genePtr$gfList      = genePtr$gfList[ !isMissingTooMuch ]
      genePtr$genomeNames = genePtr$genomeNames[ !isMissingTooMuch ]
      genePtr$geneMat     = genePtr$geneMat[ !isMissingTooMuch, ]
      numStillIn          = sum( !isMissingTooMuch )

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

      # Recalculate the number of core genes
      isCoreGene = CalcNumCoreGenes( genePtr, coreGeneThresh )
    }
  }

  # Subset to only include the core genes
  genePtr$clustList    = genePtr$clustList[ isCoreGene ]
  genePtr$genomeIdList = genePtr$genomeIdList[ isSingleCopy ]
  genePtr$geneMat      = genePtr$geneMat[ , isCoreGene ]
  cat(
    "  -- ",  numClusts, " orthologous genes were identified by cd-hit\n",
    "  -- ",  sum( isCoreGene ), " genes were used to create the tree\n",
    sep = ''
    )
}

# ------------------------------------------------------------------------------
