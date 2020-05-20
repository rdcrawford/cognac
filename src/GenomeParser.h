// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include "Genome.h"
#include <RcppParallel.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// GenomeParser
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------
// This is the functor (call operator enabling this struct to behave like a
// function), which Intel TBB will use inside the tbb::parallel_for
// function. Provides the functionally to parse
// -----------------------------------------------------------------------------

#ifndef _GENOME_PARSER_
#define _GENOME_PARSER_
struct GenomeParser : public Worker
{
  // Ctor. Gives the reference to the vector of genome class objects
  // to be parsed.
  GenomeParser( std::vector< Genome > genomeData ): genomeData( genomeData )
  { ; }

  // Vector of genome data class objects to be parsed
  std::vector< Genome > genomeData;

  // This function will be called inside tbb::parallel_for.
  // Iterates over the input range and parses the genomes
  void operator()( std::size_t begin, std::size_t end )
  {
    // Loop over the input range and parse the genome class objects
    for ( auto i = begin; i < end; ++i )
    {
      // Read in and parse the fasta and gff files
      genomeData[ i ].parseGenome();

      // Translate the amino acid sequences
      if ( !genomeData[ i ].translateSeqs() )
        Rcpp::stop( "Translate failed" );
      Rcpp::Rcout << genomeData[ i ].getNumGenes() << endl;
    }
  }
};
#endif

// -----------------------------------------------------------------------------
