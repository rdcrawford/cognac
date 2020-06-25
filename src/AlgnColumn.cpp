// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "AlgnColumn.h"
using namespace Rcpp;

// -----------------------------------------------------------------------------
// AlgnColumn
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------

// This function counts the number of gaps in a given column and
void AlgnColumn::calcColStats( double minGapFrac, unsigned int minSubThresh )
{
  // Get the number of sequences in the alignment
  unsigned int numSeqs = seqRef->size();

  // Count the character at each position in the alignment
  for ( unsigned int i = 0; i < numSeqs; i++ )
    updatCharCounts( ( *seqRef )[ i ][ cIdx ] );

  // Check if there are too many gaps at this position
  auto it = charCounts.find( '-' );
  if ( it == charCounts.end() )
  {
    isNotGappy = true;
  }
  else // If there are gaps, check that they are below the threshold
  {
    isNotGappy = ( it->second / (double) numSeqs ) <= minGapFrac;
  }

  // Check that there is the correct count of symbols at this position
  if ( charCounts.size() >= minSubThresh ) isHighDiversity = true;
  else isHighDiversity = false;
}

// Update the counts for a character in the alignemnt. If there isn't an
// entry for this character then one is added to the map
void AlgnColumn::updatCharCounts( char algnChar )
{
  // Search if this char has been observed
  auto it = charCounts.find( algnChar );

  // If this char is not in the map add it
  if ( it == charCounts.end() )
  {
    charCounts.insert( std::pair< char, int >( algnChar, 1 ) );
  }
  else // If already observed increment the count
  {
    it->second ++;
  }
}

// Return a bool to check if this postion meeths the qualifications to
// keep in the alignment
bool AlgnColumn::getColStatus()
{
  return isNotGappy && isHighDiversity;
}

// -----------------------------------------------------------------------------
