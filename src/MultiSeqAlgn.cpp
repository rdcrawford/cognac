// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "BioSeq.h"
#include "MultiSeqAlgn.h"
#include "MsaDistance.h"

// -----------------------------------------------------------------------------
// MultiSeqAlgn
// Ryan D. Crawford
// 05/14/2020
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

Rcpp::NumericMatrix MultiSeqAlgn::createDistMat( const &std::string distType )
{
  // Allocate the matrix to store the differences between aligned
  // sequences we will return
  Rcpp::NumericMatrix distMat( seqs.size(), seqs.size() );

  // Create the functor
  auto func =  msaDistance( seqs, distMat, distType );

  // To use Intel TBB we first have to create a range object. In this case
  // we want to parse all genomes in the vector
  tbb::blocked_range< int > range( 0, seqs.size() );

  // Call tbb::parallel_for, the code in the functor will be executed
  // on the availible number of threads.
  tbb::parallel_for( range, func );

  CharacterVector names = Rcpp::wrap( seqNames );
  rownames( distMat )   = names;
  colnames( distMat )   = names;

  // Set the row and column names of the matirx and Return
  return distMat;
}

// This function reads in the fasta file containing the msa and makes sure
// the file is valid for downstream analysis
void MultiSeqAlgn::parseMsa()
{
  if ( ! parseFasta() ) Rcpp::stop("Unable to open the alignment for reading\n")

  // Find the size of the alignment
  seqLen = seqs[0].size();

  // Make sure that this alignment contains sequences
  if ( !seqLen )
  {
    Rcpp::stop( "There are no sequences in this alignment...\n" );
  }

  // Iterate over the alignment and ensure that all of the sequences are
  // the same length. If not, throw and error.
  seqLen = msa[0].size();
  for ( auto it = seqs.begin(); it != seqs.end(); it++ )
  {
    // If this is not the same size as the first alignment, throw an error
    if ( it->size() != seqLen )
      Rcpp::stop( "The alignments must be the same length...\n" );
  }
}

// This function iterates over sequence in the  alignment positions until a
// gap is found. If a gap is present at this position in the alignment,
// return true. If a gap is found in none of the sequences return false.
bool MultiSeqAlgn::findGapPosition( )
{
  // Set the iterator to the first sequence
  auto it = seqIts.begin();

  // Iterate over the sequences until a gap is found
  while ( it != seqIts.end() )
  {
    // If this position in the sequnce is a gap, return true
    if ( **it == '-') return true;

    // Advance to the next reside in the sequence
    ++it;
  }
  // If we have iterated over every sequence in the alignment and not
  // found a gap, return false -- no gap found
  return false;
}

// Create a vector with iterators to each sequence in the alignment
void MultiSeqAlgn::getSeqIterators()
{
  // Allocate sufficient space for each sequence in the alignment
  strIts.reserve( seqs.size() );

  // Iterate over the alignment and get a pointer to the start of each
  // sequence in the alginment
  for ( auto it = seqs.begin(); it != seqs.end(); it++ )
    strIts.push_back( it->begin() );
}

// Iterate over each position in the alignment and remove any position with
// a gap to generte the core genome alignment
void MultiSeqAlgn::removeGaps()
{
  // Create iterators for each sequence in the alignment
  getSeqIterators();

  // Get pointers to the start and end of the alignment to keep track
  // of the current position in the alignment
  auto algnPos = strIts[0].begin();
  auto algnEnd = strIts[0].end();

  // Iterate over each position in the alignment
  while ( algnPos != algnEnd )
  {
    // If there is a gap at the curret position in the algnment, remove it
    if ( findGapPosition() )
    {
      for ( int i = 0; i < seqs.size(); i ++ )
        seqs[ i ].erase( strIts[i] );
    }
    // If there is not a gap advance to the next position
    else
    {
      for ( auto & it : strIts ) ++it;
    }
    algnPos ++;
    R_CheckUserInterrupt();
  }
}

// -----------------------------------------------------------------------------
