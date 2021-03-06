// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
// CreateCoreGenomeDistMat
// Ryan D. Crawford
// 05/15/2020
// -----------------------------------------------------------------------------
// This function reads in the path to a multiple sequence alignment.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix CreateCoreGenomeDistMat( std::string msaPath )
{
  // Create the msa class object
  MultiSeqAlgn multiSeqAlgn( msaPath );

  // Read in the fasta file and check that this is a valid alignmnet
  multiSeqAlgn.parseMsa();

  // Remove any gap positions from the Msa
  multiSeqAlgn.filterMsaColumns( 0, 1 );

  // Create a distance matrix with the raw alignment distances -- number
  // of mutations between each pair of sequences
  return multiSeqAlgn.createDistMat( "shared" );
}

// -----------------------------------------------------------------------------
