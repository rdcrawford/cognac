// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
// DeletePartitions
// Ryan D. Crawford
// 2020/06/11
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void DeletePartitions( std::string msaPath, std::vector<int> delStart,
  std::vector<int> delEnd, std::string outPath
  )
{
  // Make sure the partitions are the correct size
  if ( delStart.size() != delEnd.size() )
    Rcpp::stop( "The start positions and the end positions are not equal size");

  // Convert the start and stop positions to zero index
  for ( unsigned int i = 0; i < delStart.size(); i++ )
  {
    delStart[i] --;
    delEnd[i] --;
  }

  // Create the msa class object
  MultiSeqAlgn multiSeqAlgn( msaPath );

  // Read in the fasta file and check that this is a valid alignmnet
  multiSeqAlgn.parseMsa();

  // Erase the selected partitions from the msa
  multiSeqAlgn.deletePartitions( delStart, delEnd );

  if ( !multiSeqAlgn.writeSeqs( outPath ) )
    Rcpp::stop( "Unable to write the alignment to " + outPath );
}

// -----------------------------------------------------------------------------
