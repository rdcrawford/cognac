// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "AlgnSubCalc.h"

// -----------------------------------------------------------------------------
// CalcSubstitutionMatrix
// Ryan D. Crawford
// 2020/06/11
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix CalcAlgnSubMatrix( std::vector< std::string > seqs )
{
  // Initialize the object to
  AlgnSubCalc algnSubCalc;

  // Update the object with the
  algnSubCalc.updateSubMat( seqs );

  // Normalize the matrix by the probability of each mutation occuring
  algnSubCalc.normalizeMatrix();

  // Calculate the distance matrix and return
  return algnSubCalc.getSubMat();
}

// -----------------------------------------------------------------------------
