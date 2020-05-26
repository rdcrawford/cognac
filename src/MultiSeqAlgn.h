// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "BioSeq.h"
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MultiSeqAlgn
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------
// This class provides function for parsing and manipulating multiple sequence
// alignments. Functionality is provided to remove gaps from the alignment
// to retrieve the core genome. Then a distance matrix can be calculated
// with the paiwise distances between the sequences in the alignmnet
// -----------------------------------------------------------------------------

#ifndef _MULTI_SEQ_ALGN_
#define _MULTI_SEQ_ALGN_
class MultiSeqAlgn : public BioSeq
{
public:

  // This ctor takes the path to the msa in fasta format and creates the
  // BioSeq class object to parse the
  MultiSeqAlgn( std::string faPath ): BioSeq( faPath )
  { ; }

  // Create a distance matrix from
  Rcpp::NumericMatrix createDistMat( const std::string & distType );

  // This function reads in the fasta file containing the msa and makes sure
  // the file is valid for downstream analysis
  void parseMsa();

  // Iterate over each position in the alignment and remove any position with
  // a gap to generte the core genome alignment
  void removeGaps();

  // Create a list distance matricies for each of the individual partitions.
  // The specified by the input vector of integers.
  std::list< Rcpp::NumericMatrix > calcAlignPartitionDists(
    std::string distType, std::vector< int > genePartitions );

private:

  // Vector of iterators that each point to a sequence in the MSA
  std::vector< std::string::iterator > seqIts;

  // This function iterates over sequence in the  alignment positions until a
  // gap is found. If a gap is present at this position in the alignment,
  // return true. If a gap is found in none of the sequences return false.
  bool findGapPosition();

  // Create a vector with iterators to each sequence in the alignment
  void getSeqIterators();

  // The length of the msa
  unsigned int seqLen;
};
#endif

// -----------------------------------------------------------------------------
