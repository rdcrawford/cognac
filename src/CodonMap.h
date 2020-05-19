#include <Rcpp.h>
#include <fstream>
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------
// Codon Map
// Ryan D. Crawford
// 02/09/2019
// -----------------------------------------------------------------------------
// This object represents a map for codons to amino acids. Includeds a
// translate an input string with the nuclueotide sequence, and outputs a
// string containing the amino acid sequence.
// -----------------------------------------------------------------------------

//Define an object to store the map to look up the amino acid from codon
//sequence
#ifndef _CODON_MAP_
#define _CODON_MAP_
class CodonMap
{
public:

  //This is the only Ctor. Initializes a map of all the codons
  CodonMap();

  //Take a string and translate the sequence. Returns true if the sequences was
  //able to be translated. The input string is updated from nuclueotide to The
  //amino acid sequence.
  bool Translate( string &seq );
  
  //Take a vector of strings and translate each of them. Returns a vector
  //of the same length
  vector<string> TranslateVec( const vector<string> &ntSeqs);

  //Looks up the corresponding amino acid for a codon
  char LookUpAa( const string &codon );
  
private:

  //Map containing the codon values
  std::map< std::string, char> codonToAaMap;
};
#endif

// -----------------------------------------------------------------------------
