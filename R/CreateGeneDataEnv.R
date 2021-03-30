# ------------------------------------------------------------------------------
#  CreateGeneDataEnv
#  2020/02/05
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#' Create an environment contining the parsed data on the coding 
#' sequences in the analysis
#' @description
#'   This function initializes the environment containing data on the genes
#'   used in the analysis. This function adds the list of parsed gff files,
#'   fasta file paths, and the genome Ids of the genomes included in the 
#'   analysis to a new environment which is returned. This function uses 
#'   all available threads to parse the genomic data in parallel via the 
#'   RcppParallel package. 
#' @param featureFiles Character vector with the paths to the gff or genbank 
#'   files
#' @param fastaFiles Character vector with the paths to the fasta files
#' @param genomeIds Optional vector with the unique identifiers for the 
#'   genomes included in the analysis
#' @param outDir Directory in which to write the faa file with the translated
#'   gene sequences from all genes.
#' @return The environment containing a list of the genome features ("gfList"), 
#'   vector of gene sequences ("geneSeqs"), unique gene ids ("geneIds"), 
#'   path to the faa file written ("faaPath").
#' @export
#  -----------------------------------------------------------------------------

CreateGeneDataEnv = function( featureFiles, fastaFiles, genomeIds, outDir )
{
  # Initalize an environment to store the data on the genes. This will
  # be passed to further functions and iteratively be subset to only
  # include the genes contained in the alignment
  geneEnv = new.env()
  
  # Add the genome ids to the envrionment 
  geneEnv$genomeNames = genomeIds
  
  # Create the path to the cd-hit input file. This is where the amino acid
  # sequences for all of the genome will be written
  faaPath = paste0( outDir, "allGenes.faa" )

  # Parse the data on the input genomes
  CreateCognacRunData( geneEnv, featureFiles, fastaFiles, faaPath )

  # Add the fasta files and path to the parsed genome data
  geneEnv$fastaFiles = fastaFiles
  geneEnv$faaPath    = faaPath
  
  # Check that parsing the data was sucessful for all genomes
  hasNoGenes = sapply( geneEnv$gfList, nrow ) == 0
  if ( TRUE %in% hasNoGenes )
  {
    sink( paste0( outDir, "removed_genomes.tsv" ) )
    for ( i in which( hasNoGenes ) )
      cat( genomeNames[i], "\t", "failed parsing data", '\n')
    sink()
    
    geneEnv$fastaFiles  = fastaFiles[ !hasNoGenes ] 
    geneEnv$gfList      = geneEnv$gfList[ !hasNoGenes ]
    geneEnv$genomeNames = geneEnv$genomeNames[ !hasNoGenes ]
    
    cat(
      "  -- Failed to parse the data for ", sum( hasNoGenes ), 
      " genomes\n  -- Continuing", length( geneEnv$genomeNames ), "genomes\n",
      sep = ''
      )
    
  } else {
    
    geneEnv$fastaFiles = fastaFiles
  }
  
  cat(
    "  -- Parsed the data for ", length( geneEnv$geneIds ), 
    " coding genes\n", 
    sep = ''
    ) 
  
  return( geneEnv )
}

# ------------------------------------------------------------------------------
