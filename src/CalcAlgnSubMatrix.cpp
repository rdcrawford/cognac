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

  return algnSubCalc.getSubMat();
}

// -----------------------------------------------------------------------------
