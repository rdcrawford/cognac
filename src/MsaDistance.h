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
// MSA distances
// function. An MSA is input ant for each sequence in the MSA the distcance
// is calcuated between all of the other sequences and sotred in a matrix
// of integers. The process is done in parralel using the RcppParallel
// library, which provides thead safe wrappers for the R data types.
// -----------------------------------------------------------------------------

// ---- Define a worker class for parallel computing ---------------------------

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
  DistFunction distuncPtr;

  // Function call operator that work from the range specified by begin and end
  // This enables this object t0 behave like a function when passed to
  // parallel_for.
  void operator()( tbb::blocked_range< int > & range );

  //
  double calcSharedDist( const std::string &ref, const std::string &qry );
};

// -----------------------------------------------------------------------------
