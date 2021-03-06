// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "BioSeq.h"
#include "MultiSeqAlgn.h"
#include "MsaDistance.h"
#include "AlgnColumn.h"
#include "AlgnSubCalc.h"
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
  MsaDistance msaDistance( seqs, distMat );

  // Set the function pointer to the requested distance function
  msaDistance.setDistFunc(distType );

  // Call tbb::parallel_for, the code in the functor will be executed
  // on the availible number of threads.
  parallelFor( 0, seqs.size(), msaDistance );

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
    Rcpp::stop( "There are no sequences in this alignment...\n" );

  for ( auto it = seqs.begin(); it != seqs.end(); it++ )
  {
    // If this is not the same size as the first alignment, throw an error
    if ( it->size() != seqLen )
      Rcpp::stop( "The alignments must be the same length...\n" );
  }
}

// Create a vector with iterators to each sequence in the alignment
std::vector< std::string::iterator > MultiSeqAlgn::getSeqIterators()
{
  // Allocate sufficient space for each sequence in the alignment
  std::vector< std::string::iterator > seqIts;
  seqIts.reserve( seqs.size() );

  // Iterate over the alignment and get a pointer to the start of each
  // sequence in the alginment
  for ( auto it = seqs.begin(); it != seqs.end(); it++ )
    seqIts.push_back( it->begin() );

  return seqIts;
}

// Iterate over each position in the alignment and remove any position with
// a gap to generte the core genome alignment

void MultiSeqAlgn::filterMsaColumns( double minGapFrac,
  int minSubThresh, std::vector<int> & genePositions
  )
{
  // Create iterators for each sequence in the alignment
  auto seqIts = getSeqIterators();

  // Check each column in the alignment there there is at least one
  // subsitiution and there there is less than 50% gaps
  std::vector< AlgnColumn > algnCols;
  algnCols.reserve( seqLen );
  for ( unsigned int i = 0; i < seqLen; i++ )
    algnCols.push_back( AlgnColumn(& seqs, i ) );

  // Use tbb parallel for to iterate over each column in the alingmnet and
  // check that each position is of sufficient quality
  tbb::parallel_for_each(
    algnCols.begin(),
    algnCols.end(),
    [&] ( AlgnColumn & algnColumn )
  {
    algnColumn.calcColStats( minGapFrac, minSubThresh );
  });

  // Allocate a new vector with the new alignment sequnces
  std::vector< std::string > filterSeqs( seqs.size() );

  // Reserve enough space for each string to encompass the entire alignments
  for ( auto & seq : filterSeqs ) seq.reserve( seqLen );

  // Iterate over each position in the alignment and copy over each char
  // that meets the criteria for inclusion in the alignment
  for ( auto it = algnCols.begin(); it != algnCols.end(); it++ )
  {
    // If there are too many gaps or no variants at this position in the
    // alignment, advance the iterators
    if ( it->getColStatus() )
    {
      // Get the current column
      auto colIdx = std::distance( algnCols.begin(), it );

      // Add this position to the filtered alignment
      for ( unsigned int i = 0; i < seqs.size(); i ++ )
        filterSeqs[ i ].push_back( seqs[ i ][ colIdx ] );
    }

    R_CheckUserInterrupt();
  }

  // Get rid of any unused space
  for ( auto & seq : filterSeqs ) seq.shrink_to_fit();

  // Reseq the sequences in the msa to the filtered sequences
  seqs = filterSeqs;

  // If a vector of gene positions was input, update the vector so it now
  // reflects the gene partitions that were erased
  if ( genePositions.size() )
  {
    int algnIdx   = 0; // Position in the alignment
    int numErased = 0; // Nnumber of positions in the alignemt that were erased

    // Iterate over each position in the alignment
    for ( auto & pos : genePositions )
    {
      while ( algnIdx <= pos )
      {
        algnIdx ++;
        if ( !algnCols[ algnIdx ].getColStatus() ) numErased ++;
      }

      pos -= numErased;

      R_CheckUserInterrupt();
    }
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

  // Allocate the matrix to store the differences between aligned
  // sequences we will return
  Rcpp::NumericMatrix distMat( numSeqs, numSeqs );

  // Initialize the functor with the entire alignment
  MsaDistance msaDistance( seqs, distMat );

  // Set the function pointer to the requested distance function. If the
  // function type requires, the alignment substitution probabilities are
  // calculated once on the entire alignment
  msaDistance.setDistFunc( distType );

  // For each set of patitions
  for ( unsigned  int i = 0; i < genePartitions.size(); i++ )
  {
    if ( i == 0 ) gStart = 0;
    else gStart = genePartitions[ i - 1 ];
    len = genePartitions[ i ] - gStart;

    // Create an alignment with only the sequences in the partitions
    std::vector< std::string > subSeqs( numSeqs );
    for ( unsigned int j = 0; j < numSeqs; j++ )
      subSeqs[ j ] = seqs[ j ].substr( gStart, len );

    // Set the sequences to calculate the distance for to the subsequence
    // corresponding to this partition
    msaDistance.setMsaSeqs( subSeqs );

    // Call tbb::parallel_for, the code in the functor will be executed
    // on the availible number of threads.
    parallelFor( 0, numSeqs, msaDistance );

    // Set the row and column names
    rownames( distMat ) = names;
    colnames( distMat ) = names;

    // Set the row and column names of the matirx and Return
    distMatList.push_back( Rcpp::clone( distMat ) );
  }

  return distMatList;
}

// Delete a selected partitions in the alignment
void MultiSeqAlgn::deletePartitions( const std::vector<int> &delStart,
  const std::vector<int> &delEnd )
{
  // Get an iterator with the beginning of each string
  auto seqIts = getSeqIterators();

  // Iterate over the partitions to be delected and remove them from each
  // sequence in the alignement
  int delLen = 0;

  // Iterate over each partions to delete
  for ( unsigned int i = 0; i < delStart.size(); i++ )
  {
    for ( unsigned int j = 0; j < seqs.size(); j++ )
    {
      // Find the position to erase and delete the sequence from the alignmnt
      auto curStart = seqIts[ j ] + ( delStart[ i ] - delLen );
      auto curEnd = seqIts[ j ] + ( delEnd[ i ] - delLen + 1 );
      seqs[ j ].erase( curStart, curEnd );
    }

    // Keep track of how much of the sequence has been deleted to adjust the
    // start positions
    delLen += delEnd[i] - delStart[i] + 1;
  }
}

std::vector< double > MultiSeqAlgn::calcAlgnQualScores(
  std::string method, int stepVal, int windowSize
  )
{
  // Initialize the start of the gene partition
  unsigned int numSeqs = getNumSeqs();

  // Allocate the matrix to store the differences between aligned
  // sequences we will return
  Rcpp::NumericMatrix distMat( numSeqs, numSeqs );

  // Initialize the functor with the entire alignment
  MsaDistance msaDistance( seqs, distMat );

  // Set the function pointer to the requested distance function. If the
  // function type requires, the alignment substitution probabilities are
  // calculated once on the entire alignment
  msaDistance.setDistFunc( method );

  //
  unsigned int startPos = 0;
  unsigned int endPos   = windowSize;

  //
  std::vector< double > minDistVals;

  while ( endPos < seqLen )
  {
    // Create an alignment with only the sequences in the partitions
    std::vector< std::string > subSeqs( numSeqs );
    for ( unsigned int j = 0; j < numSeqs; j++ )
      subSeqs[ j ] = seqs[ j ].substr( startPos, windowSize );

    // Set the sequences to calculate the distance for to the subsequence
    // corresponding to this partition
    msaDistance.setMsaSeqs( subSeqs );

    // Call tbb::parallel_for, the code in the functor will be executed
    // on the availible number of threads.
    parallelFor( 0, numSeqs, msaDistance );


    double minVal = distMat( 0 , 1 );
    for ( int i = 0; i < distMat.nrow(); i++ )
     for ( int j = i + 1; j < distMat.ncol(); j++ )
       if ( distMat( i , j ) < minVal ) minVal = distMat( i , j );


    // Set the row and column names of the matirx and Return
    minDistVals.push_back( minVal );

    // Increment he start position by the window size
    startPos += stepVal;
    endPos   = startPos + windowSize;
  }

  return minDistVals;
}

// -----------------------------------------------------------------------------
