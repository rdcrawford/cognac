// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "GenomeFeatures.h"
#include "CodonMap.h"
#include "BioSeq.h"

// -----------------------------------------------------------------------------
// Genome
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------
// This class inherits the propertes of the BioSeq class, representing the
// whole genome sequence, and the GenomeFeatures class which is the parsed
// gff file. This allows easy parsing of the fasta file to retrieve the
// gen sequences and translate them to return the amino acid sequences.
// -----------------------------------------------------------------------------

#ifndef _GENOME_
#define _GENOME_
class Genome: public BioSeq, public GenomeFeatures
{
public:

  // Default Ctor
  Genome()
  { ; }

  // Value ctor: takes the
  Genome( const std::string &faPath, const std::string &gffPath,
    const std::string &genomeId ):
    BioSeq( faPath ), GenomeFeatures( gffPath, genomeId )
  { ; }

  Genome( const std::string &faPath, vector<char> strand,
    vector<int> startPos, vector<int> endPos, vector<int> contig ):
    BioSeq( faPath ), GenomeFeatures( strand, startPos, endPos, contig )
  { ; }

  // This functions takes and integer correponding the the current gene,
  // and a string corresponding to a gene sequence to be updated, both passed
  // by reference. The current gene index is incremented with each function
  // call and the seq variable is updated to
  bool getGeneSeq( string &seq );

  // Parse the inut files
  void parseGenome();

  // If the default ctor was used this function inputs the vaiables
  // for
  bool parseGenome( const std::string &faPath, const std::string &gffPath,
    const std::string &genomeId );

  // Translate the coding sequences for this genome. The internal variable
  // "aaSeqs" is updated
  bool translateSeqs();

  // Returns the vector of translated gene sequences
  std::vector< std::string > getAaSeqs();

  // Returns a vector to the translate gene sequences
  std::vector< std::string > *getAaSeqRef();

private:

  // Integer to keep track of the current gene
  int gIdx = 0;

  // Vector of amino acid sequences
  std::vector< std::string> aaSeqs;

  // Pass a gene nt sequence by reference and update the sequence
  // to the reverse complement
  void getReverseCompliment( std::string &sequence );

};
#endif

// -----------------------------------------------------------------------------
