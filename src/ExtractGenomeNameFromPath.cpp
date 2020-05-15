# include <Rcpp.h>

// -----------------------------------------------------------------------------
// ExtractGenomeNameFromPath
// 2019/07/10
// Ryan D. Crawford
// -----------------------------------------------------------------------------
// The cd-hit output file follows the rast syntax for naming genes:
// "fig|genomeName.peg". This function extracts the genome name from the gene
// name.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::string ExtractGenomeNameFromPath( std::string path )
{
  int start = path.find_last_of( "/" ) + 1;
  int len = path.find_last_of("." ) - start;
  return path.substr( start, len );
}

// [[Rcpp::export]]
std::string GetGenomeNameWithExt( std::string path, std::string ext )
{
  int start = path.find_last_of( "/" ) + 1;
  int len = path.find( ext ) - start;
  return path.substr( start, len );
}

// ----------------------------------------------------------------------------
