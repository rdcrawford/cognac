# ------------------------------------------------------------------------------
# CreateGeneData
# 2019/12/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function creates a new environment with the the sequences
# of genes used as input to cd-hit. The unique gene ids 
# are also stored in a seperate vector 
# ------------------------------------------------------------------------------

CreateGeneData = function( genePtr )
{
  # Add vectors with the gene sequences and the unique gene ids 
  genePtr$geneSeqs = ParseFasta( genePtr$cdHitInput )
  genePtr$geneIds  = sapply(
    names(genePtr$geneSeqs), ExtractGeneNames, USE.NAMES = FALSE
    )
  
  # Remove any non-core genes from the data
  coreGenes        = unlist( genePtr$clustList )
  isCoreGene       = genePtr$geneIds %in% coreGenes
  genePtr$geneSeqs = genePtr$geneSeqs[ isCoreGene ]
  genePtr$geneIds  = genePtr$geneIds[ isCoreGene ]
  
  # Create meta-data on the genes in the 
  CreateGeneMetaData( genePtr )
}

# ------------------------------------------------------------------------------