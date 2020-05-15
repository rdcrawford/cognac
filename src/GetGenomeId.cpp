#include <Rcpp.h>
using namespace std;

// [[Rcpp::interfaces(r, cpp)]]

// -------------------------------------------------------------------------- //
// Parse Fasta Files By Position                                              //
// Ryan D. Crawford                                                           //
// 02/09/2019                                                                 //
// -------------------------------------------------------------------------- //
// This functions updates an input sequence to the reverse compliment.        //
// The case of the string is maintained.                                      //
// -------------------------------------------------------------------------- //

// [[Rcpp::export]]
std::string GetGenomeId( std::string inStr )
{
  int start = inStr.find( "fig|" ) + 4;
  int len  = inStr.find( ".peg" ) - start;
  return inStr.substr( start, len );
}

// -------------------------------------------------------------------------- //
