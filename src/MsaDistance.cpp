// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "MsaDistance.h"
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MsaDistance
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------

MsaDistance::MsaDistance(
  const std::vector<std::string> &msa,
  Rcpp::NumericMatrix distMat, std::string distFunType)
{
  // Assign the inputs to the member variables
  this->msa     = msa;
  this->distMat = distmat;

  switch( distFunType )
  {
    case "raw":
      distFunction = &MsaDistance::calcRawDist;
      break;
    case "shared":
      distFunction = &MsaDistance::calcSharedDist;
      break;
    default:
      Rcpp::stop( "Distance function type: ", distFunType,
        " is not supported\nSupported types are:\n  -- raw\n  -- shared\n",
        )
  }
}

// Function call operator that work from the range specified by begin and end
void MsaDistance::operator()( tbb::blocked_range< int > & range )
{
  double distVal;
  for ( int i = range.begin(); i < range.end(); ++i)
  {
    for ( int j = 0; j < i; ++j )
    {
      // Value to store the count of mismatches between the two sequences
      if ( j == i ) distVal = 0;
      else distVal = ( this->*distFunction )( msa[ i ], msa[ j ] );

      // Assign the position in the distance matrix to the number
      // of mutations
      distMat( i, j ) = distVal;
      distMat( j, i ) = distVal;
    }
  }
}

double MsaDistance::calcRawDist(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a counter for the number of mutations between two sequences
  double numMutations = 0;

  // For each base in the two sequence, see if there is NOT a match at
  // the kth position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for (int i = 0; i < ref.length(); i++)
  {
    // If this position in the reference is not the same as the
    // query..
    if ( ref[i] != qry[i] )
    {
      // And neither sequence has a gap at this position, increment
      // the counter for the numer of numations
      if ( ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
      {
        numMutations ++;
      }
    }
  }
  return numMutations;
}

double MsaDistance::calcSharedDist(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a counters for the number of mutations between two sequences
  // and a counter for the number of sites that the two sequences share
  double numMutations = 0;
  double numSites     = 0;

  // For each base in the two sequence, see if there is NOT a match at
  // the kth position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for (int i = 0; i < ref.length(); i++)
  {
    // If this position in the reference is not the same as the
    // query..
    if ( ref[i] != qry[i] )
    {
      // And neither sequence has a gap at this position, increment
      // the counter for the numer of numations
      if ( ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
      {
        numMutations ++;
        numSites ++;
      }
    }
  }
  return numMutations / numSites;
}
