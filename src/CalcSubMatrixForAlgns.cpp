// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "AlgnSubCalc.h"
#include "MultiSeqAlgn.h"

// -----------------------------------------------------------------------------
// CalcSubMatrixForAlgns
// Ryan D. Crawford
// 2020/06/11
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

Rcpp::NumericMatrix CalcSubMatrixForAlgns( std::vector< std::string > msaPaths )
{
  // Initialize the object to
  AlgnSubCalc algnSubCalc;

  for ( auto & msa : msaPaths )
  {
    // Create the msa class object
    MultiSeqAlgn multiSeqAlgn( msa );

    // Read in the fasta file and check that this is a valid alignmnet
    multiSeqAlgn.parseMsa();

    // Update the object with the
    algnSubCalc.updateSubMat( multiSeqAlgn.getSeqs() );
  }

  // Normalize the matrix by the probability of each mutation occuring
  algnSubCalc.normalizeMatrix();

  return algnSubCalc.getSubMat();
}

// -----------------------------------------------------------------------------
