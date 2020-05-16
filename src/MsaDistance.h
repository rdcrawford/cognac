// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MsaDistance
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------
// This is a functor that provides a parallelized calculation of the
// MSA distances. For each pair of sequences in the alignment, distcance
// is calcuated as the number of mutations between the two sequences.
// The process is done in parralel using tbb via RcppParallel. The distance
// function to be used is passed as argument into the ctor and a function
// pointer is used in the call operator to specify the function when called.
// -----------------------------------------------------------------------------

#ifndef _MSA_DISTANCE_
#define _MSA_DISTANCE_
struct MsaDistance : public Worker
{
  // Ctor: Initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type). Additionally
  // The integer indicating the type of function to sue
  MsaDistance( const std::vector<std::string> &msa,
    Rcpp::NumericMatrix distMat, std::string distFunType );

  // Input multipe sequence alignment to calculate distance from
  const std::vector< std::string > msa;

  // Output matrix to write distances to
  RMatrix< double > distMat;

  // There are several ways to calculate the distance from an MSA. This
  // provides a function pointer to be called when creating the distance
  // matrx.
  typedef double ( MsaDistance::*DistFunction )( const std::string &refSeq,
    const std::string &qrySeq );
  DistFunction distFunction;

  // Function call operator that work from the range specified by begin and end
  // This enables this object t0 behave like a function when passed to
  // parallel_for.
  void operator()( tbb::blocked_range< int > & range );

  // Retrun the fraction of mutations per shared site. Any position where
  // there is a gap in either sequece is excluded.
  double calcSharedDist( const std::string &ref, const std::string &qry );

  // Return the raw number of mutations between two sequences
  double calcRawDist( const std::string &ref, const std::string &qry );
};
#endif

// -----------------------------------------------------------------------------
