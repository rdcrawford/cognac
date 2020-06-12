// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "BioSeq.h"

// -----------------------------------------------------------------------------
// ParseFasta
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------
// This functions uses the "BioSeq" class object to read in the fasta file.
// The data in the fasta file is returned into R as a character vector
// with the contig names
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::CharacterVector ParseFasta( const std::string &faPath )
{
  // Create the "wholeGenomeSeq" class object of ject for he inut genome
  BioSeq bioSeq( faPath );

  // If the fasta file was unable to be parsed, throw an error
  if ( !bioSeq.parseFasta() ) Rcpp::stop( "Unable to read: ", faPath );

  // Convert the contigs to a character vector
  Rcpp::CharacterVector seqs = Rcpp::wrap( bioSeq.getSeqs() );

  seqs.attr("names") = Rcpp::wrap( bioSeq.getSeqNames() );

  // Return the contigs as a vector
  return seqs;
}

// -----------------------------------------------------------------------------
