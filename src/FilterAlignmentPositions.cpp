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
//'   This function removes gaps and/or positions without sufficient varition
//'   from the alingment.
//' @param msaPath Path to the alignment
//' @param filterMsaPath Path to write the filtered alignment
//' @param minGapFrac Double representing the minimium fraction of gaps to
//'   remain in the alignment.Defaults to 0.01.
//' @param minSubThresh Integer representing the minimum number of
//'   substitutions to remain in the alignment. This is the number of instances
//'   of any minor allele to remain. Defaults to 0.
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

// -----------------------------------------------------------------------------
