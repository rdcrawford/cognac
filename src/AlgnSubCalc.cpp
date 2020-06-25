// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include "AlgnSubCalc.h"
using namespace Rcpp;

// -----------------------------------------------------------------------------
// CalcSubstitutionMatrix
// Ryan D. Crawford
// 2020/06/11
// -----------------------------------------------------------------------------

AlgnSubCalc::AlgnSubCalc()
{
  // Initialize a vector of the amino acid characters
  std::vector< char > aaCodes
  {
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
  };

  // Create a map to store the frequency of each amino acid.
  for ( auto it = aaCodes.begin(); it != aaCodes.end(); it ++)
    aaCounts.insert( std::pair< char, int >( *it, 0 ) );

  // Create a map to look up the row/column indices of amino acids
  for ( unsigned int i = 0; i < aaCodes.size(); i++ )
    aaIdxs.insert( std::pair< char, int >( aaCodes[i], i ) );

  // Set number of rows and columns to attribute dim
  subMat = NumericMatrix( aaCodes.size(), aaCodes.size() );

  Rcpp::CharacterVector names = Rcpp::wrap( aaCodes );
  rownames( subMat ) = names;
  colnames( subMat ) = names;
}

// Return the subsitiution matrix with the alignment
Rcpp::NumericMatrix AlgnSubCalc::getSubMat()
{
  return subMat;
}

void AlgnSubCalc::countAaCodes( const std::string &seq )
{
  // For each amino acid in the alignment, increment the count
  for ( auto it = seq.begin(); it != seq.end(); it++ )
  {
    auto aaCount = aaCounts.find( *it );
    if ( aaCount != aaCounts.end() ) aaCount->second ++;
  }
}

bool AlgnSubCalc::getAaIdx( char aa, unsigned int &idx )
{
  auto it = aaIdxs.find( aa );
  if ( it == aaIdxs.end() ) return false;
  idx = it->second;
  return true;
}

void AlgnSubCalc::updateSubMat(
  const std::vector< std::string > &seqs
  )
{
  // For each sequence get the frequency of each amino acid and the
  // frequeny at which an amino acid is mutated
  unsigned int numSeqs = seqs.size();
  for ( unsigned int i = 0; i < numSeqs; i++ )
  {
    // Count the frequency of each amino acid in this sequence
    countAaCodes( seqs[i] );

    // For each Additional sequence,
    for ( unsigned int j = i + 1; j < numSeqs; j++ )
      countSubs( seqs[ i ], seqs[ j ] );
  }
}

void AlgnSubCalc::normalizeMatrix()
{
  // Transform the matrix. Each non-zero entry is the log of the
  // frequency of the mutation divided by the probability of the amino
  // acid ( log Mij / pi )
  auto aaCount = aaCounts.begin();
  for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  {
    // For each Additional sequence,
    for ( unsigned int j = 0; j < subMat.ncol(); j++ )
    {
      if ( subMat( i, j ) )
        subMat( i, j ) = -log( subMat( i, j ) / aaCount->second );
    }
  }
}

void AlgnSubCalc::countSubs( const std::string &ref, const std::string &qry )
{
  unsigned int rIdx;
  unsigned int qIdx;

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
        if ( getAaIdx( ref[i], rIdx ) && getAaIdx( qry[i], qIdx ) )
        {
          subMat( rIdx, qIdx ) ++;
          subMat( qIdx, rIdx ) ++;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
