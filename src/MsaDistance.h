// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "AlgnSubCalc.h"
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MsaDistance
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------
// This is a functor that provides a parallelized calculation of the
// MSA distances. For each pair of sequences in the alignment, disttance
// is calcuated as the number of mutations between the two sequences.
// The struct is a functor used by tbb via RcppParallel. The name of the
// distancefunction to be used is passed as argument into the ctor and a
// function pointer is used in the call operator to specify the distance
// function to be used.
// -----------------------------------------------------------------------------

#ifndef _MSA_DISTANCE_
#define _MSA_DISTANCE_
struct MsaDistance : public Worker
{
  // Ctor: Initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type). Additionally
  // The integer indicating the type of function to sue
  MsaDistance(
    const std::vector<std::string> &msa, Rcpp::NumericMatrix distMat):
    msa(msa), distMat(distMat)
  { ; }

  // Input multipe sequence alignment to calculate distance from
  std::vector< std::string > msa;

  // Output matrix to write distances to
  RMatrix< double > distMat;

  // Object the calculates the substitution matrix between amino acids or
  // nucleotides in the alignment
  AlgnSubCalc algnSubCalc;

  // Function call operator that work from the range specified by begin and
  // end
  void operator()( std::size_t begin, std::size_t end );

  // Set the function pointer to the type specified by the input
  // argument "distFunType"
  void setDistFunc( std::string distFunType );

  // Set the sequences contained in this alignment.
  void setMsaSeqs( const std::vector< std::string > & msa );

  // This function cacluates the subsiion probabilities between sequences
  // in the alignmet
  void calcSubProbabilities();

  // Returns the raw number of mutations between two sequences
  double calcRawDist( const std::string &ref, const std::string &qry );

  // Returns the sum of the log liklihood of substitutions between two
  // sequences in the alignment
  double calcSharedDist( const std::string &ref, const std::string &qry );

  // Returns the log odds of the substitutions between the two sequences
  double calcNormProbDist( const std::string &ref, const std::string &qry );

  // Calculate the loglilihood of each position in the alignment -- similar
  // to blossum distance
  double calcBlosum( const std::string &ref, const std::string &qry );

  // There are several ways to calculate the distance from an MSA. This
  // provides a function pointer to be called when creating the distance
  // matrx. By default, calculate the raw number of substitutions
  typedef double ( MsaDistance::*DistFunction )( const std::string &refSeq,
    const std::string &qrySeq );
  DistFunction distFunction = &MsaDistance::calcRawDist;
};
#endif

// -----------------------------------------------------------------------------
