// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
// CalcDistMatsAtPartions
// Ryan D. Crawford
// 05/15/2020
// -----------------------------------------------------------------------------
// This function reads in the path to a multiple sequence alignment, and Then
// generates a distance matrix at each of the partions specified in the input
// vector, which specifies the end of the partitions.  
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::list < Rcpp::NumericMatrix > CalcAlgnPartitionDists(
  std::string msaPath,  std::string method, std::vector< int > genePartitions
  )
{
  // Create the msa class object
  MultiSeqAlgn multiSeqAlgn( msaPath );

  // Read in the fasta file and check that this is a valid alignmnet
  multiSeqAlgn.parseMsa();

  // Create a distance matrix with the raw alignment distances -- number
  // of mutations between each pair of sequences
  return multiSeqAlgn.calcAlignPartitionDists( method, genePartitions );
}

// -----------------------------------------------------------------------------
