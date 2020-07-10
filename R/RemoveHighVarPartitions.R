# ------------------------------------------------------------------------------
# Benchmark Parse Cd-Hit
# 2020/04/24
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

RemoveHighVarPartitions = function( algnData, algnType, outAlgnPath )
{
  # Set the default to analyze the amino acid alignment 
  if ( missing(algnType) ) algnType = "aa"
  
  # Constants: ro indexes of the gene partition matrix
  START = 1
  END   = 2
  
  # Get the column index of the gene partitions and the path to the alignment
  if ( algnType == "aa" )
  {
    cIdx = 3
    algn = algnData$aaAlgnPath
  } else if ( algnType == "nt" ) {
    cIdx = 4
    algn = algnData$ntAlgnPath
  } else {
    stop("Type must be \"aa\" or \"nt\" ")
  }
  
  # Get the positions of the genes in the alignment
  geneParitions = sapply( algnData$geneData[ , cIdx ], 
    function(x) as.numeric( strsplit(x, '-')[[1]] )
    )
  
  # Calculate the distances for each partition in the ali
  distMatList = CalcAlgnPartitionDists( algn, "raw", geneParitions[ END, ] )
  
  # Get the length of the alignment for each gene
  algnLens =  sapply( 1:ncol(geneParitions),
    function(j) geneParitions[ END, j ] - geneParitions[ START, j ] + 1
    )
  
  # Get the total number of variants in the alignment, normalized 
  # to the length of the gene 
  numVars  = sapply( 1:length(distMatList), 
    function(i) sum( as.dist( distMatList[[i]] ) / algnLens[i] )
    )
  
  # Calculate summary statistics for the distribution
  meanVal = mean( numVars )
  stdDev  = sd( numVars )
  
  # Identify outliers
  cutOffVal = stdDev * 3
  upperVal  = meanVal + cutOffVal
  isOutlier = numVars > upperVal
  
  if ( TRUE %in% isOutlier )
  {
    # Get the begin and end positions of the partitions to erase
    delStartPos = geneParitions[ START, isOutlier ]
    delEndPos   = geneParitions[ END, isOutlier ]
    
    # Delete the outlier partitions from the alignment 
    DeletePartitions( algn, delStartPos, delEndPos, outAlgnPath )
  }
  
  # Calculate the boundaries of the genes
  genePartitions      = character( length( algnLens ) )
  genePartitions[ 1 ] = paste0( "1-", algnLens[1] )
  
  curEnd = algnLens[1]
  for ( i in 2:length( algnLens ) )
  {
    if ( isOutlier[i] )
    {
      genePartitions[ i ] = NA
    } else {
      
      curStart = curEnd + 1
      curEnd   = curStart + algnLens[ i ] - 1
      
      # Find the partitions of the genes in the nt alignment
      genePartitions[ i ] = paste0( curStart, "-", curEnd )
    }
  }
  
  # Add the gene partitions for the updated alignment to the dataframe
  cbind.data.frame( algnData$geneData,  genePartitions )
}

# ------------------------------------------------------------------------------