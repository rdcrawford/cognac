#include <Rcpp.h>
#include "CodonMap.h"
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// -------------------------------------------------------------------------- //
// Codon Map                                                                  //
// Ryan D. Crawford                                                           //
// 02/09/2019                                                                 //
// -------------------------------------------------------------------------- //
// This object represents a map for codons to amino acids. Includeds a        //
// translate an input string with the nuclueotide sequence, and outputs a     //
// string containing the amino acid sequence.                                 //
// -------------------------------------------------------------------------- //

//Ctor
CodonMap::CodonMap()
{
  //Initialize the codon map to the
  codonToAaMap =
  {
    //Isoleucine
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'},

    //Leucine
    {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'}, {"TTA", 'L'},
    {"TTG", 'L'},

    //Valine
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},

    //Phenylalanine
    {"TTT", 'F'}, {"TTC", 'F'},

    //Methionine
    {"ATG", 'M'},

    //Cysteine
    {"TGT", 'C'}, {"TGC", 'C'},

    //Alanine
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},

    //Glycine
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'},

    //Proline
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},

    //Threonine
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},

    //Serine
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'}, {"AGT", 'S'},
    {"AGC", 'S'},

    //Tyrosine
    {"TAT", 'Y'}, {"TAC", 'Y'},

    //Tryptophan
    {"TGG", 'W'},

    //Glutamine
    {"CAA", 'Q'}, {"CAG", 'Q'},

    //Asparagine
    {"AAT", 'N'}, {"AAC", 'N'},

    //Histidine
    {"CAT", 'H'}, {"CAC", 'H'},

    //Glutamic acid
    {"GAA", 'E'}, {"GAG", 'E'},

    //Aspartic acid
    {"GAT", 'D'}, {"GAC", 'D'},

    //Lysine
    {"AAA", 'K'}, {"AAG", 'K'},

    //Arginine
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'}, {"AGA", 'R'},
    {"AGG", 'R'},

    //Stop codons
    {"TGA", 'Z'}, {"TAA", 'Z'}, {"TAG", 'Z'}
    };

}

//Look up values in the map
char CodonMap::LookUpAa( const std::string &codon )
{
  char aa = codonToAaMap[ codon ];
  if (aa == '\0') return 'X';
  return aa;
}

//Take a string and translate the sequence. Returns true if the sequences was
//able to be translated. The input string is updated from nuclueotide to The
//amino acid sequence.
bool CodonMap::Translate( std::string &seq )
{
  //If the input sequence is not a multiple of 3 return an NA
  if (!seq.length() / 3) return false;

  // Initalize a string to store the translated sequence
  // Look up the stop codon. If there is not a cannonical stop
  // codon, translate the entire sequence.
  char stopCodon = LookUpAa( seq.substr( seq.size() - 4, 3) ) ;
  unsigned int aaSeqLen;
  if ( stopCodon == 'Z' )
  {
    aaSeqLen = seq.length() / 3 - 1;
  }
  else
  {
    aaSeqLen = seq.length() / 3;
  }
  std::string aaSeq(aaSeqLen, 'X');

  //Don't look up the first amino acid, it always codes for Met
  aaSeq[0] = 'M';

  //For the second codon throught the length of the
  //gene codon look up the amino acid
  for( unsigned int i = 1; i < aaSeqLen; i++ )
  {
    int startPos = i * 3;
    aaSeq[i] = LookUpAa( seq.substr(startPos, 3) );
  }

  //Update the seq variable (passed by reference)
  seq = aaSeq;

  return true;
}

vector<string> CodonMap::TranslateVec( const vector<string> &ntSeqs )
{
  //Initialize the vector to return
  vector<string> aaSeqs( ntSeqs.size() );

  //Translate each nt seq
  for ( unsigned int i = 0; i < ntSeqs.size(); i++ )
  {
    string seq = ntSeqs[ i ];

    //Seq is passed by reference to the translate function. If it is able
    //to be translated, add it to the vector
    if ( Translate(seq) )
    {
      aaSeqs[ i ] = seq;
    }
  }

  return aaSeqs;
}


// -------------------------------------------------------------------------- //
