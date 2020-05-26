// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "BioSeq.h"
#include <RcppParallel.h>
#include "MultiSeqAlgn.h"
#include "MsaDistance.h"
using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// MultiSeqAlgn
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------

Rcpp::NumericMatrix MultiSeqAlgn::createDistMat( const std::string & distType )
{
  // Allocate the matrix to store the differences between aligned
  // sequences we will return
  Rcpp::NumericMatrix distMat( getNumSeqs(), getNumSeqs() );

  // Create the functor
  MsaDistance msaDistance( seqs, distMat, distType );

  // Call tbb::parallel_for, the code in the functor will be executed
  // on the availible number of threads.
  parallelFor( 0, getNumSeqs(), msaDistance );

  Rcpp::CharacterVector names = Rcpp::wrap( seqNames );
  rownames( distMat )         = names;
  colnames( distMat )         = names;

  // Set the row and column names of the matirx and Return
  return distMat;
}

// This function reads in the fasta file containing the msa and makes sure
// the file is valid for downstream analysis
void MultiSeqAlgn::parseMsa()
{
  if ( ! parseFasta() )
   Rcpp::stop("Unable to open the alignment for reading\n");

  // Find the size of the alignment
  this->seqLen = seqs[0].size();

  // Make sure that this alignment contains sequences
  if ( !seqLen )
  {
    Rcpp::stop( "There are no sequences in this alignment...\n" );
  }

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
    if ( **it == '-' ) return true;

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
  seqIts.reserve( seqs.size() );

  // Iterate over the alignment and get a pointer to the start of each
  // sequence in the alginment
  for ( auto it = seqs.begin(); it != seqs.end(); it++ )
    seqIts.push_back( it->begin() );
}

// Iterate over each position in the alignment and remove any position with
// a gap to generte the core genome alignment
void MultiSeqAlgn::removeGaps()
{
  // Create iterators for each sequence in the alignment
  getSeqIterators();

  // Iterate over each position in the alignment
  for ( unsigned int algnPos = 0; algnPos < seqs[0].size(); algnPos ++ )
  {
    // If there is a gap at the curret position in the algnment, remove it
    if ( findGapPosition() )
    {
      for ( unsigned int i = 0; i < seqs.size(); i ++ )
        seqs[ i ].erase( seqIts[i] );
    }
    // If there is not a gap advance to the next position
    else
    {
      for ( auto & it : seqIts ) ++it;
    }

    R_CheckUserInterrupt();
  }
}


std::list< Rcpp::NumericMatrix > MultiSeqAlgn::calcAlignPartitionDists(
  std::string distType, std::vector< int > genePartitions
  )
{
  // Initialize the list, reserving the number of matricies to output
  std::list< Rcpp::NumericMatrix > distMatList;

  // Initialize the start of the gene partition
  unsigned int gStart;
  unsigned int numSeqs = getNumSeqs();
  unsigned int len;

  // Convert the row names to an Rcpp character vector to be used as the
  // row and column names
  Rcpp::CharacterVector names = Rcpp::wrap( seqNames );

  for ( unsigned int i = 0; i < genePartitions.size(); i++ )
  {
    if ( i == 0 ) gStart = 0;
    else gStart = genePartitions[ i - 1 ];
    len = genePartitions[ i ] - gStart;

    // Allocate the matrix to store the differences between aligned
    // sequences we will return
    Rcpp::NumericMatrix distMat( numSeqs, numSeqs );

    // Create an alignment with only the sequences in the partitions
    std::vector< std::string > subSeqs( numSeqs );
    for ( unsigned int j = 0; j < numSeqs; j++ )
      subSeqs[ j ] = seqs[ j ].substr( gStart, len );

    // Create the functor
    MsaDistance msaDistance( subSeqs, distMat, distType );

    // Call tbb::parallel_for, the code in the functor will be executed
    // on the availible number of threads.
    parallelFor( 0, numSeqs, msaDistance );

    // Set the row and column names
    rownames( distMat ) = names;
    colnames( distMat ) = names;

    // Set the row and column names of the matirx and Return
    distMatList.push_back( distMat );
  }

  return distMatList;
}

// -----------------------------------------------------------------------------
