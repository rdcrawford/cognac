// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include "MsaDistance.h"
#include "AlgnSubCalc.h"
using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MsaDistance
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------

// Set the function pointer to the type specified by the input
// argument "distFunType"
void MsaDistance::setDistFunc( std::string distFunType )
{
  if ( distFunType == "raw" )
  {
    distFunction = &MsaDistance::calcRawDist;
  }
  else if ( distFunType ==  "shared" )
  {
    distFunction = &MsaDistance::calcSharedDist;
  }
  else if ( distFunType ==  "normProb" )
  {
    distFunction = &MsaDistance::calcNormProbDist;

    // Read in the alignment and calcualte the pairwise substitutions in they
    // alignment
    algnSubCalc.updateSubMat( msa );

    // Normalize the matirx to get the log liklihood
    algnSubCalc.calcNormalizedProbs();
  }
  else if ( distFunType ==  "logLike" )
  {
    distFunction = &MsaDistance::calcBlosum;

    // Read in the alignment and calcualte the pairwise substitutions in the
    // alignment
    algnSubCalc.updateSubMat( msa );

    // Normalize the matirx to get the log liklihood
    algnSubCalc.calcLogLikelihoods();
  }
  else
  {
    std::string errStr = "Distance function type " + distFunType +
      " is not supported\nSupported types are:\n  -- raw\n  -- shared\n" +
      "  -- normProb\n  -- logLike";
    Rcpp::stop( errStr );
  }
}

// Set the sequences contained in this alignment.
void MsaDistance::setMsaSeqs( const std::vector< std::string > & msa )
{
  this->msa = msa;
}

// Calculate the raw number of mutations between two sequences
double MsaDistance::calcRawDist(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a counter for the number of mutations between two sequences
  double numMutations = 0;

  // For each base in the two sequence, see if there is NOT a match at
  // the ith position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for ( unsigned int i = 0; i < ref.size(); i++ )
  {
    // If this position in the reference is not the same as the
    // query..
    if ( ref[i] != qry[i] )
    {
      // And neither sequence has a gap at this position, increment
      // the counter for the numer of numations
      if ( ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
      {
        numMutations += 1;
      }
    }
  }
  return numMutations;
}

// Function call operator that work from the range specified by begin and end
void MsaDistance::operator()(std::size_t begin, std::size_t end)
{
  // Initialize the value that stores the distances
  double distVal;

  // FInd the pairwise differences for each sequence
  for ( std::size_t i = begin; i < end; i++ )
  {
    for ( std::size_t j = 0; j < i; j++ )
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

// Returns the number of mutations normalized to the number of shared
// sites (excluding gap potitions)
double MsaDistance::calcSharedDist(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a counters for the number of mutations between two sequences
  // and a counter for the number of sites that the two sequences share
  double numMutations = 0;
  double numSites     = 0;

  // For each base in the two sequence, see if there is NOT a match at
  // the ith position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for ( unsigned int i = 0; i < ref.size(); i++ )
  {
    // If this position in the reference is not the same as the
    // query..
    if (  ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
    {
      // And neither sequence has a gap at this position, increment
      // the counter for the numer of numations
      if ( ref[i] != qry[i] )
      {
        numMutations += 1.0;
      }
      numSites += 1.0;
    }
  }
  return numMutations / numSites;
}

double MsaDistance::calcNormProbDist(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a value to store the calulated distance between the two
  // sequences
  double seqDist = 0;

  // For each base in the two sequence, see if there is NOT a match at
  // the ith position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for ( unsigned int i = 0; i < ref.size(); i++ )
  {
    // If this position in the reference is not the same as the
    // query..
    if (  ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
    {
      // And neither sequence has a gap at this position, increment
      // the counter for the numer of numations
      if ( ref[i] != qry[i] )
      {
        seqDist += algnSubCalc.getSubPr( ref[i], qry[i] );
      }
    }
  }
  return seqDist;
}

double MsaDistance::calcBlosum(
  const std::string &ref, const std::string &qry
  )
{
  // Initialize a value to store the calulated distance between the two
  // sequences
  double seqDist = 0;

  // For each base in the two sequence, see if there is NOT a match at
  // the ith position of of the alignment and there is not an aligned,
  // base at that position increment the number of mutations.
  for ( unsigned int i = 0; i < ref.size(); i++ )
  {
    // If this position in the reference is not the same as the
    // query..
    if (  ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
    {
      seqDist += algnSubCalc.getSubPr( ref[i], qry[i] );
    }
  }
  return seqDist;
}

// -----------------------------------------------------------------------------
