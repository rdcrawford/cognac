// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// AlgnSubCalc
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

#ifndef _ALGN_SUB_CALC_
#define _ALGN_SUB_CALC_
class AlgnSubCalc
{
public:

  // Default ctor
  AlgnSubCalc();

  // Update the substitiution matrix with the input alignment
  void updateSubMat( const std::vector< std::string > &seqs );

  // Return the subsitiution matrix with the alignment
  Rcpp::NumericMatrix getSubMat();

  // Divide the counts of
  void normalizeMatrix();

private:

  // Map of amino acids to their frequency
  std::map< char, int > aaCounts;

  // Map to look up the index of the amino acids
  std::map< char, int > aaIdxs;

  // Matrix with the substitution frequencies
  Rcpp::NumericMatrix subMat;

  // Count the chars in a sequence
  void countAaCodes( const std::string &seq );

  // Count the amino acid substitutions between two sequences. The counts
  // for each amino acid are incrementd for each char in the reference
  // sequence.
  void countSubs( const std::string &ref, const std::string &qry );

  //
  bool getAaIdx( char aa, unsigned int &idx );

};
#endif

// -----------------------------------------------------------------------------
