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

// Define a static vector to stand in for the vector of integers
// so that it can be passed as argument by default
static std::vector<int> DEFAULT_VECTOR;

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
  void filterMsaColumns( double minGapFrac,  int minSubThresh,
     std::vector<int> & genePositions=DEFAULT_VECTOR );

  // Create a list distance matricies for each of the individual partitions.
  // The specified by the input vector of integers.
  std::list< Rcpp::NumericMatrix > calcAlignPartitionDists(
    std::string distType, std::vector< int > genePartitions );

  // Delete a selected partitions in the alignment
  void deletePartitions( const std::vector<int> &delStart,
    const std::vector<int> &delEnd );

private:

  // Create a vector with iterators to each sequence in the alignment
  std::vector< std::string::iterator > getSeqIterators();

  // The length of the msa
  unsigned int seqLen;
};
#endif

// -----------------------------------------------------------------------------
