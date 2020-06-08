// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <fstream>
#include <algorithm>
#include "GenomeFeatures.h"
#include "CodonMap.h"
#include "BioSeq.h"
#include "Genome.h"
using namespace std;


// -----------------------------------------------------------------------------
// Genome
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------

void Genome::parseGenome()
{
  // Read in the fasta file
  if ( !parseFasta() )
    Rcpp::stop( "Failed parsing fasta file: " + faPath );

  // Parse the relevant attributes of the gff file to the component vectors
  if ( !parseGfs( this ) )
    Rcpp::stop( "Failed parsing gff file: " + gfPath );

}

bool Genome::parseGenome(
  const std::string &faPath,
  const std::string &gffPath,
  const std::string &genomeId
  )
{
  if ( !setFasta( faPath ) ) return false;
  return setGff( gffPath, genomeId, this );
}

// Update the string passed by reference to that of the current gene.
// If the sequnce was updated, return true.
bool Genome::getGeneSeq( std::string &seq )
{
  // Make sure that a gene at the current index exists
  if ( gIdx > startPos.size() ) return false;

  // Get the sequence at the coordinates of this entry in the gff file
  bool isUpdated =
    getSeqAtCoord( contig[ gIdx ], startPos[ gIdx ], endPos[ gIdx ], seq );

  // If these were not valid coordinates, return false indicating that
  // the sequnce was not updated
  if ( !isUpdated )
  {
    // Increment the index for the next gene
    gIdx ++;
    return false;
  }

  // If this is the reverse strand, get the reverse compliment
  if ( strand[ gIdx ].compare("-") == 0 ) getReverseCompliment( seq );
  // Rcpp::Rcout << seq << std::endl;
  // Increment the index for the next gene
  gIdx ++;

  // Retur true to indicate that the gene was updated
  return true;
}

void Genome::getReverseCompliment( std::string &sequence )
{
  // Reverse the iput string
  reverse( sequence.begin(), sequence.end() );

  // Iterate over the sequence and substiture the complementary base
  for ( unsigned int i = 0; i < sequence.length(); i ++ )
  {
    if      ( sequence[i] == 'A' ) sequence[i] = 'T';
    else if ( sequence[i] == 'T' ) sequence[i] = 'A';
    else if ( sequence[i] == 'C' ) sequence[i] = 'G';
    else if ( sequence[i] == 'G' ) sequence[i] = 'C';
    else                           sequence[i] = 'N';
  }
}

// Translate the coding sequences for this genome. The internal variable
// "aaSeqs" is updated. Returns true if all of the genes were translated
bool Genome::translateSeqs()
{
  // Initialize the codon map object to look up the corresponding amino acid
  // for each codon
  CodonMap codonMap;

  std::string seq; // Empty string to populate with the gene seqences
  auto numGenes = featId.size();

  // Allocate a vector with the number of genes
  aaSeqs.reserve( numGenes );

  // Iterate over all of the genes
  while ( gIdx < numGenes )
  {
    // Get the sequnece of the current gene
    if ( !getGeneSeq( seq ) ) return false;

    // Translate the nucleotide sequence
    if ( !codonMap.Translate( seq ) ) return false;

    // If the gene was updated and translated, add it to the vector of
    // amino acid sequences
     aaSeqs.push_back( seq );
  }

  return true;
}

// Returns the vector of translated gene sequences
std::vector< std::string > Genome::getAaSeqs()
{
  return aaSeqs;
}

// Returns a vector to the translate gene sequences
std::vector< std::string > *Genome::getAaSeqRef()
{
  return & aaSeqs;
}
