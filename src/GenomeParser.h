// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "Genome.h"
using namespace std;

// -----------------------------------------------------------------------------
// GenomeParser
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------
// This is the functor (call operator enabling this struct to behave like a
// function), which Intel TBB will use inside the tbb::parallel_for
// function. Provides the functionally to parse
// -----------------------------------------------------------------------------

struct GenomeParser : public worker
{
  // Ctor. Gives the reference to the vector of genome class objects
  // to be parsed.
  GenomeParser( std::vector< Genome > &genomeData ): genomeData( genomeData )
  { ; }

  // Vector of genome data class objects to be parsed
  std::vector< Genome > genomeData;

  // This function will be called inside tbb::parallel_for.
  // Iterates over the input range and parses the genomes
  void operator()( tbb::blocked_range< int > & range )
  {
    // Loop over the input range and parse the genome class objects
    for ( int i = range.begin(); i < range.end(); ++i )
    {
      // Read in and parse the fasta and gff files
      genomeData[ i ].parseGenome();

      // Translate the amino acid sequences
      if ( !genomeData[ i ].translateSeqs() )
        Rcpp::stop(
          "Translating genes for ", genomeData[i].getGenomeId(), " failed"
          );
    }
  }
};

// -----------------------------------------------------------------------------
