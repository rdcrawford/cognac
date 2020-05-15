// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "BioSeq.h"

// -----------------------------------------------------------------------------
// MultiSeqAlgn
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

class MultiSeqAlgn : public BioSeq
{
public:

  // This ctor takes the path to the msa in fasta format and creates the
  // BioSeq class object to parse the
  MultiSeqAlgn( std::string faPath ):BioSeq( faPath )
  { ; }

  // Create a distance matrix from
  Rcpp::NumericMatrix createDistMat( const &std::string distType );

  // This function reads in the fasta file containing the msa and makes sure
  // the file is valid for downstream analysis
  void parseMsa();

  // Iterate over each position in the alignment and remove any position with
  // a gap to generte the core genome alignment
  void removeGaps();

private:

  // Vector of iterators that each point to a sequence in the MSA
  std::vector< std::string::iterator > strIts;

  // This function iterates over sequence in the  alignment positions until a
  // gap is found. If a gap is present at this position in the alignment,
  // return true. If a gap is found in none of the sequences return false.
  bool findGapPosition();

  // Create a vector with iterators to each sequence in the alignment
  void getSeqIterators();

};

// -----------------------------------------------------------------------------
