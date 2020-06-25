// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
//  Filter Alignment Positions
//  Ryan D. Crawford
//  05/15/2020
//  ----------------------------------------------------------------------------
//' @name FilterAlignmentPositions
//' @title Filter Alignment Positions
//' @description
//'   This function
//' @param msaPath Path to the alignment
//' @param filterMsaPath Path to write the filtered alignment
//' @param double minGapFrac=0.01,
//' @param minSubThresh Minimium number o subsititutions
//' @param genePositions Optional vector of partitions within the alignment
//' @return void
//' @export
//  ----------------------------------------------------------------------------



// [[Rcpp::export]]
void FilterAlgnPositions(
  std::string msaPath, std::string filterMsaPath, double minGapFrac=0.01,
  int minSubThresh=0
  )
{
  // Create the msa class object
  MultiSeqAlgn multiSeqAlgn( msaPath );

  // Parse the msa
  multiSeqAlgn.parseMsa();

  // Remove columns from the msa that dont meet the parameters
  multiSeqAlgn.filterMsaColumns( minGapFrac, minSubThresh );

  // Write the MSA to a new fle
  if ( !multiSeqAlgn.writeSeqs( filterMsaPath ) )
    Rcpp::stop( "Writing to file " + filterMsaPath + "failed..." );
}

// [[Rcpp::export]]
std::vector< int > FilterPartitionedAlgnPositions(
  std::string msaPath, std::string filterMsaPath,
  std::vector<int> genePositions, double minGapFrac=0.01, int minSubThresh=0
  )
{
  // Create the msa class object
  MultiSeqAlgn multiSeqAlgn( msaPath );

  // Parse the msa
  multiSeqAlgn.parseMsa();

  // Remove columns from the msa that dont meet the parameters
  multiSeqAlgn.filterMsaColumns( minGapFrac, minSubThresh, genePositions );

  // Write the MSA to a new fle
  if ( !multiSeqAlgn.writeSeqs( filterMsaPath ) )
    Rcpp::stop( "Writing to file " + filterMsaPath + "failed..." );

  return genePositions;
}

// -----------------------------------------------------------------------------
