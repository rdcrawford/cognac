// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "GenomeData.h"

// -----------------------------------------------------------------------------
// CreateCognacRunData
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------
// This function takes the R enviroment used to store variables containing
// data used in the analysis. This function sets up the congnac run: 1) parsing
// gff files and fasta files, 2) translating the coding sequences and retirves
// the correspondin gene ids, and 3) writes the coding sequences for the
// cd-hit input file
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void CreateCognacRunData(
  Rcpp::Environment                &geneEnv, // R environment to be updated
  const std::vector< std::string > &gfPaths, // Paths to the gff files
  const std::vector< std::string > &faPaths, // Paths to the fasta files
  const std::string                &faaPath  // Paths to the faa to create
  )
{
  Rcpp::Rcout << "CreateCognacRunData step 1: get genome names" << std::endl;
  // Retrieve the genome names from the environment
  std::vector< std::string > genomeIds = geneEnv[ "genomeNames" ];

  // Intialize a vector of genome class objects to parse
  GenomeData genomeData( gfPaths, faPaths, genomeIds );
  // Rcpp::Rcout << "CreateCognacRunData step 3: write the faa" << std::endl;
  // Write the gene sequences to an faa file to be the input for cd-hit
  genomeData.writeFaa( faaPath );
  // Rcpp::Rcout << "CreateCognacRunData step 4: assign the data" << std::endl;
  // Parse the genome features to a data frames
  geneEnv.assign( "gfList", genomeData.createGeneDataFrames() );
  geneEnv.assign( "geneSeqs", genomeData.getAaSeqs() );
  geneEnv.assign( "geneIds", genomeData.getGeneIds() );
}

// -----------------------------------------------------------------------------
