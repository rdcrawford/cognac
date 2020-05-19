// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
//  CreateAlgnDistMat
//  Ryan D. Crawford
//  05/15/2020
//  ----------------------------------------------------------------------------
//' @name CreateAlgnDistMat
//' @title  Create Algnment Distance Matrix
//' @description
//'   This function reads in the path to a multiple sequence alignment. The
//'   pairwise alignment distances for each sequence are returned. The 
//'   pairwise distances are calculated in parallel via the RcppParallel 
//'   package.
//' @param Path to the alignment
//' @param Method for calculating distance: "raw", or "shared."
//' @param Logical to specify whether to remove gap positions from the alignment
//'   to create the core genome.
//' @return A numeric matrix
//' @export
//  ----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix CreateAlgnDistMat( 
  std::string msaPath, std::string method, bool isCore
  )
{
  // Create the msa class object
  MultiSeqAlgn multiSeqAlgn( msaPath );

  // Read in the fasta file and check that this is a valid alignmnet
  multiSeqAlgn.parseMsa();

  // If requested, remove the gaps to create the "core genome alignment" 
  if ( isCore ) multiSeqAlgn.removeGaps();
    
  // Create a distance matrix with the raw alignment distances -- number
  // of mutations between each pair of sequences
  return multiSeqAlgn.createDistMat( method );
}

// -----------------------------------------------------------------------------
