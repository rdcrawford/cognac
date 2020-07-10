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


  // Return the probability of substitutions between two symbols
  double getSubPr( char rCh, char qCh );

  // Calculate the log liklihood of each substitution in the matrix
  void calcLogLikelihoods();

  // Calculate the log normalized substitution probabilities
  void calcNormalizedProbs();


private:

  // Map of amino acids to their frequency
  std::map< char, int > aaCounts;

  // Map to look up the index of the amino acids
  std::map< char, int > aaIdxs;

  // Matrix with the substitution frequencies
  Rcpp::NumericMatrix subMat;

  // Bool indicating if the matrix has been normalized
  bool isNormalized = false;

  // Count the chars in a sequence
  void countAaCodes( const std::string &seq );

  // Count the amino acid substitutions between two sequences. The counts
  // for each amino acid are incrementd for each char in the reference
  // sequence.
  void countSubs( const std::string &ref, const std::string &qry );

  // Look up the row/column index of the input char. If this is a valid aa/nt
  // symbol true is returned and the "idx" variable is updated.
  bool getAaIdx( char aa, unsigned int &idx );

  // Calculate the subsitiution probabilities for each pair of amino acids
  void calcSubProbabilities();

};
#endif

// -----------------------------------------------------------------------------
