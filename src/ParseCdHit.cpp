// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "CdHitParser.h"

// -----------------------------------------------------------------------------
// Parse CD Hit
// Ryan D. Crawford
// 11/20/2019
// -----------------------------------------------------------------------------
// This function takes the output file for CD-HIT and creates an R list
// class object containing a matrix of presence of the sequences
// and a list with the gene identifiers for each of the sequences
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void ParseCdHit(
  const std::string &cdHitClstFile, // Path to the cd-hit results
  Rcpp::Environment &geneEnv, // Environment with  the data on all of the genes
  bool isBinary, // If true, creates a binary matrix, false the percent ids
  int  minClustCount // Remove genes at low frequency in the cluster
  )
{
  // This creates a cd-hit object with lists for the gene cluster data And
  // a matrix of gene presence absence. "geneEnv" is passed by reference And
  // updated with these attributes
  CdHitParser cdHitParser( cdHitClstFile, isBinary, removeLowFreq, geneEnv );
}

// -----------------------------------------------------------------------------
