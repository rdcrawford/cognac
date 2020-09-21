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

void AlgnSubCalc::calcSubProbabilities()
{
  // for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  // {
  //   for ( unsigned int j = 0; j < subMat.ncol(); j++ )
  //   {
  //     Rcout << subMat( i, j ) << ' ';
  //   }
  //   Rcout << std::endl;
  // }

  double numSubs;
  for ( unsigned int i = 0; i < subMat.nrow(); i++ )
    for ( unsigned int j = 0; j < subMat.ncol(); j++ )
      numSubs += subMat( i, j );
  Rcout << std::endl << "numSubs: " << numSubs
        << std::endl << std::endl << std::endl;

  for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  {
    for ( unsigned int j = 0; j < subMat.ncol(); j++ )
    {
      if ( subMat( i, j ) )
      {
        subMat( i, j ) = subMat( i, j ) / numSubs;
      }
    }
  }
  // for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  // {
  //   for ( unsigned int j = 0; j < subMat.ncol(); j++ )
  //   {
  //     Rcout << subMat( i, j ) << ' ';
  //   }
  //   Rcout << std::endl;
  // }
}

// Calculate the log liklihood of each substitution in the matrix
void AlgnSubCalc::calcLogLikelihoods()
{
  // Calculate the subsitiution probabilities for each pair of amino acids
  calcSubProbabilities();

  // Get the total counts of the amino acids
  int aaCount = 0;
  for ( auto it = aaCounts.begin(); it != aaCounts.end(); it ++ )
    aaCount += it->second;

  for ( auto it = aaCounts.begin(); it != aaCounts.end(); it ++ )
  {
    Rcout << it->first << ": "<< it->second << std::endl;
  }
  Rcout << std::endl;
  // Initialize an iterator for the first amino acid count
  std::vector< double > aaProbs;
  for ( auto it = aaCounts.begin(); it != aaCounts.end(); it ++ )
    aaProbs.push_back( ( double ) it->second / aaCount );

  for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  {
    // For each Additional sequence,
    for ( unsigned int j = 0; j < subMat.ncol(); j++ )
    {
      if ( subMat( i, j ) )
      {
        subMat( i, j ) =
          log( ( subMat( i, j ) / ( aaProbs[i] * aaProbs[j] ) ) );
      }
    }
  }
  // Rcout << std::endl<< std::endl<< std::endl<< std::endl;
  // for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  // {
  //   for ( unsigned int j = 0; j < subMat.ncol(); j++ )
  //   {
  //     Rcout << subMat( i, j ) << ' ';
  //   }
  // }
}

// Calculate the log normalized substitution probabilities
void AlgnSubCalc::calcNormalizedProbs()
{
  // Calculate the subsitiution probabilities for each pair of amino acids
  calcSubProbabilities();

  for ( unsigned int i = 0; i < subMat.nrow(); i++ )
  {
    // For each Additional sequence,
    for ( unsigned int j = 0; j < subMat.ncol(); j++ )
    {
      if ( subMat( i, j ) )
      {
        subMat( i, j ) = -log( subMat( i, j ) );
      }
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
    if ( ref[i] != '-' && qry[i] != '-' && ref[i] != 'N' &&  qry[i] != 'N' )
    {
      if ( getAaIdx( ref[i], rIdx ) && getAaIdx( qry[i], qIdx ) )
      {
        subMat( rIdx, qIdx ) ++;
        subMat( qIdx, rIdx ) ++;
      }
    }
  }
}

double AlgnSubCalc::getSubPr( char rCh, char qCh )
{
  // Initialize the row and column indicies
  unsigned int rIdx;
  unsigned int qIdx;

  // Look up the row and column indicies in the matrix
  if ( getAaIdx( rCh, rIdx ) && getAaIdx( qCh, qIdx ) )
  {
    return subMat( rIdx, qIdx );
  }
  return 0;
}


// -----------------------------------------------------------------------------
