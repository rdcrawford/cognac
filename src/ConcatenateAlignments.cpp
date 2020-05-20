#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Concatenate Alignments
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------
// This function takes two sorted alignments and concatenates them. 
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void ConcatenateAlignments( 
  Rcpp::StringVector &concatAlgn,
  const Rcpp::StringVector &algn
  )
{
  // Check that the input alignment are the same length 
  if (concatAlgn.size() != algn.size() ) 
  {
    std::cout << "Alignments must be the same length!\n" 
              << "Length of the concatenated gene sequence: "
              << concatAlgn.size() 
              << std::endl
              << "Length of the alignment: "
              << algn.size()
              << std::endl;
    exit( 1 );
  }
  
  // Iterate over the concatented gene alignment and append the single 
  // alignment from the single gene alignment
  for ( int i = 0; i < concatAlgn.size( );  i++ ) concatAlgn[ i ] += algn[ i ];
}

// -----------------------------------------------------------------------------