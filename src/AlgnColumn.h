// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// AlgnColumn
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------
// This class defines the column of an alignment and specifys
// -----------------------------------------------------------------------------

#ifndef _ALGN_COLUMN_
#define _ALGN_COLUMN_
class AlgnColumn
{
public:

  // Value ctor: takes the reference to an MSA and
  AlgnColumn( std::vector< std::string > *seqRef, int cIdx ):
    seqRef( seqRef ), cIdx( cIdx )
  { ; }

  // This function counts the number of gaps in a given column and
  void calcColStats( double minGapFrac, unsigned int minSubThresh );

  // Return a bool to check if this postion meeths the qualifications to
  // keep in the alignment. Returns true if this position is
  bool getColStatus();

private:

  // Reference to the alignment being analyzed
  std::vector< std::string > *seqRef;

  // Index of this column in the alignment
  int cIdx;

  // Map with the counts of all of the chars included in the alignments
  std::map< char, int > charCounts;

  // Bool indicating that
  bool isNotGappy;

  // Bool indicating that there is insufficient variation in this column
  bool isHighDiversity;

  // Update the counts for a character in the alignemnt. If there isn't an
  // entry for this character then one is added to the map
  void updatCharCounts( char algnChar );

  // Calcuate the number of sequences with the minor allele
  int getNumMinorAlleles();
};
#endif

// -----------------------------------------------------------------------------
