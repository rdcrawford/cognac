// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "MsaDistance.h"
using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MsaDistance
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------

void MsaDistance::setDistFunc()
{
  if ( distFunType == "raw" )
  {
    distFunction = &MsaDistance::calcRawDist;
  }
  else if ( distFunType ==  "shared" )
  {
    distFunction = &MsaDistance::calcSharedDist;
  }
  else
  {
    std::string errStr = "Distance function type " + distFunType +
      " is not supported\nSupported types are:\n  -- raw\n  -- shared\n";
    Rcpp::stop( errStr );
  }
}

double MsaDistance::calcRawDist(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a counter for the number of mutations between two sequences
  double numMutations = 0;
  setDistFunc();

  // For each base in the two sequence, see if there is NOT a match at
  // the kth position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for ( unsigned int i = 0; i < ref.length(); i++)
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

// void MsaDistance::operator()(std::size_t begin, std::size_t end)
// {
//   for (std::size_t i = begin; i < end; i++)
//   {
//     for (std::size_t j = 0; j < i; j++)
//     {
//       // Strings we will operate on
//       std::string ref = msa[ i ];
//       std::string qry = msa[ j ];
//
//       // Initialize a counters for the number of mutations between two sequences
//       // and a counter for the number of sites that the two sequences share
//       double numMutations = 0;
//       double numSites     = 0;
//
//       // For each base in the two sequence, see if there is NOT a match at
//       // the kth position of of the alignment and there is not an aligned,
//       // base at that position increment the number of mutations.
//       for ( unsigned int i = 0; i < ref.length(); i++ )
//       {
//         // If this position in the reference is not the same as the
//         // query..
//         if ( ref[i] != qry[i] )
//         {
//           // And neither sequence has a gap at this position, increment
//           // the counter for the numer of numations
//           if ( ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&
//            qry[i] != 'N' )
//           {
//             numMutations ++;
//             numSites ++;
//           }
//         }
//       }
//       distMat( i, j ) = numMutations / numSites;
//       distMat( j, i ) = numMutations / numSites;
//     }
//   }
// }

// Function call operator that work from the range specified by begin and end
void MsaDistance::operator()(std::size_t begin, std::size_t end)
{
  setDistFunc();
  double distVal;
  for (std::size_t i = begin; i < end; i++)
  {
    for (std::size_t j = 0; j < i; j++)
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
  for ( unsigned int i = 0; i < ref.length(); i++ )
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
