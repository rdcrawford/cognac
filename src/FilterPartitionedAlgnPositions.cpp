// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
//  Filter Partitioned Algn Positions
//  Ryan D. Crawford
//  05/15/2020
//  ----------------------------------------------------------------------------
//' @name FilterPartitionedAlgnPositions
//' @title Filter Partitioned Algn Positions
//' @description
//'   This function
//' @param msaPath Path to the alignment
//' @param filterMsaPath Path to write the filtered alignment
//' @param genePositions Vector of gene partitions in the alignment
//' @param minGapFrac Double representing the minimium fraction of gaps to
//'   remain in the alignment.Defaults to 0.01.
//' @param minSubThresh Integer representing the minimum number of
//'   substitutions to remain in the alignment. This is the number of instances
//'   of any minor allele to remain. Defaults to 0.
//' @return Updated vector with the partions in the filtered alignment
//' @export
//  ----------------------------------------------------------------------------

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
