// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <fstream>
#include <algorithm>
#include "CdHitParser.h"
using namespace std;
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Parse CD Hit
// Ryan D. Crawford
// 11/20/2019
// -----------------------------------------------------------------------------

// Ctor
CdHitParser::CdHitParser(
  const std::string &cdHitClstFile, // Path to the cd-hit results
  bool isBinary,                    // Return binary matrix
  int clSizeThesh,                  // Remove low frequency genes
  Rcpp::Environment &geneEnv        // Environment to append results to
  )
{
  // Set the class member variables
  genomeIds = Rcpp::as<std::vector<std::string> >( geneEnv[ "genomeNames" ] );
  this->clSizeThesh = clSizeThesh; // Set the minimium cluster size to process

  ifstream    ifs;         // Input file stream for the cd-hit output
  std::string line;        // Currnt line in the text file
  int         lineNum = 0; // Counter for the current line

  // Open the cd-hit results for reading
  ifs.open( cdHitClstFile.c_str() );
  if ( !ifs.is_open() )
    Rcpp::stop("Cannot open the cd-hit results...");

  // Read in the results  line by line
  while ( getline( ifs, line ) )
  {
    // If this is the start of a new cluster, mark this as a header
    if ( line[0] == '>' )
    {
      headers.push_back( lineNum );
    }
    // Append the line number -- the cluster will start on the line
    // after a header
    else
    {
      // Append the line to a vector
      cdHitResults.push_back( line );

      // Increment the current postion in the cd hit results
      lineNum ++;
    }
  }

  // Find the positions of all of the fasta headers in the results
  for ( auto it = headers.begin() + 1; it != headers.end(); ++it )
  {
    // The end of the cluster is the position before the next header
    tailers.push_back( *it - 1 );
  }
  // The final end position is the last string in the vector
  tailers.push_back( cdHitResults.size() - 1 );

  // Assign the number of cd-hit clusters to the environment
  geneEnv.assign( "nCogs", headers.size() );

  // If requested remoeve any low frequency clusters
  if ( clSizeThesh > 1 )
  {
    if (  !RemoveLowFreqClusts() )
      Rcpp::stop(
        "No clusters with sufficient numbers of genes were identified..."
        );
  }

  // Transform the raw results to a list of vectors containing the names
  // of the genes
  CreateClustList();

  if ( isBinary )
  {
    // Free the memory associated with the cd-hit results
    cdHitResults.clear();

    // Create a genome x gene matrix with the presene or absene of each cluster
    geneEnv.assign( "geneMat", CreateBinaryMat() );

  } else {

    // Create the list of gene identities to the reference sequence
    CreateItentList();

    // Free the memory associated with the cd-hit results
    cdHitResults.clear( );

    // Create a genome x gene matrix with the presene or absene of each cluster
    geneEnv.assign( "geneMat", CreateIdentMat() );
  }

  geneEnv.assign( "clustList", clustList );
  geneEnv.assign( "genomeIdList", gIdList );
}

bool CdHitParser::RemoveLowFreqClusts( )
{
  auto headerIt  = headers.begin();      // Iterator for the cluster start
  auto tailerIt  = tailers.begin();      // Iterator for the cluster ends
  int  clustSize;                        // Number of genes in the cluster

  // Remove any low frequency genes from the data
  while ( headerIt != headers.end() )
  {
    // Calculate the number of genes in the custer
    clustSize = ( *tailerIt - *headerIt ) + 1;

    //   If there are fewer than the required number of genes remove them
    if ( clustSize < clSizeThesh )
    {
      headers.erase( headerIt );
      tailers.erase( tailerIt );

    } else {

      // Move to the next cluster in the cd-hit resuls
      headerIt ++;
      tailerIt ++;
    }
    R_CheckUserInterrupt();
  }

  if ( cdHitResults.size() ) return true;
  return false;
}

std::string CdHitParser::GetGeneName( const std::string &inStr )
{
  int start = inStr.find( "fig|" );
  int len  = inStr.find( "..." ) - start;
  return inStr.substr( start, len );
}

std::string CdHitParser::GetGenomeId( const std::string &inStr )
{
  int start = inStr.find( "fig|" ) + 4;
  int len  = inStr.find( ".peg" ) - start;
  return inStr.substr( start, len );
}

void CdHitParser::CreateClustList( )
{
  // Set the size of the list to the number of clusters
  clustList.resize( headers.size() );
  gIdList.resize( headers.size() );

  auto headerIt = headers.begin();      // Iterator for the cluster start
  auto tailerIt = tailers.begin();      // Iterator for the cluster ends
  auto clustIt  = clustList.begin();    // It for the list of gene IDs
  auto gIdIt    = gIdList.begin();      // It for the list of genomes IDs

  // Iterate over each header
  while ( headerIt != headers.end() )
  {
    // For each line in the cluster fill in the vectors for the
    // gene IDs and genome IDs for the respective list element
    for ( int i = *headerIt; i <= *tailerIt; i++ )
    {
      // Extract the name of the gene from the line in the results
      std::string geneId = GetGeneName( cdHitResults[i] );

      // Update the vectors in the gene id and genome ID
      clustIt->push_back( geneId );
      gIdIt->push_back( GetGenomeId( geneId ) );
    }

    // Move to the next cluster in the cd-hit resuls
    headerIt ++;
    tailerIt ++;
    clustIt ++;
    gIdIt ++;

    R_CheckUserInterrupt();
  }
}

double CdHitParser::GetGeneIdent( const std::string &line )
{
  if ( line.find("... *") != string::npos ) return( 100.00 );
  int start = line.find_last_of("at") + 2;
  int len   = line.find_last_of('%') - start;
  return stod( line.substr( start, len ) );
}

void CdHitParser::CreateItentList( )
{
  // Set the size of the list to the number of clusters
  clustIdentList.resize( headers.size() );

  auto headerIt = headers.begin();        // Iterator for the cluster start
  auto tailerIt = tailers.begin();        // Iterator for the cluster ends
  auto itentIt  = clustIdentList.begin(); // It for the list of gene IDs

  // Iterate over each header
  while ( headerIt != headers.end() )
  {
    // For each line in the cluster fill in the vectors for the
    // gene IDs and genome IDs for the respective list element
    for ( int i = *headerIt; i <= *tailerIt; i++ )
    {
      // Extract the identity of the gene from the line in the results
      itentIt->push_back( GetGeneIdent( cdHitResults[i] ) );
    }

    // Move to the next cluster in the cd-hit resuls
    headerIt ++;
    tailerIt ++;
    itentIt ++;

    R_CheckUserInterrupt();
  }
}

int CdHitParser::GetGenomeIdx( const string & name )
{
  auto it = std::find( genomeIds.begin(), genomeIds.end(), name );
  if ( it == genomeIds.end() ) return -1;
  return it - genomeIds.begin();
}

Rcpp::NumericMatrix CdHitParser::CreateBinaryMat(  )
{
  // Initialize the matirx to output
  Rcpp::NumericMatrix geneMat( genomeIds.size(), gIdList.size() );

  // Iterate over each header
  for ( auto gIdIt = gIdList.begin(); gIdIt != gIdList.end(); gIdIt ++ )
  {
    // Look up the position in the list for the column index
    int colIdx = std::distance( gIdList.begin(), gIdIt );

    // For each line in the cluster fill in the vectors for the
    // gene IDs and genome IDs for the respective list element
    for ( auto it = gIdIt->begin(); it < gIdIt->end(); it++ )
    {
      int rIdx = GetGenomeIdx( *it );
      if ( rIdx != -1 ) geneMat( rIdx, colIdx ) = 1;
    }

    // Move to the next list element
    R_CheckUserInterrupt();
  }

  rownames( geneMat ) = Rcpp::wrap( genomeIds );
  return geneMat;
}

Rcpp::NumericMatrix CdHitParser::CreateIdentMat(  )
{
  // Initialize the matirx to output
  Rcpp::NumericMatrix geneMat( genomeIds.size(), gIdList.size() );

  // Create an iterator for the list of identities
  auto itentIt = clustIdentList.begin();

  // Iterate over each header
  for ( auto gIdIt = gIdList.begin(); gIdIt != gIdList.end(); gIdIt++ )
  {
    // Calculate the index in the list
    int colIdx = std::distance( gIdList.begin(), gIdIt );

    // Create an iterator for the vecotor of percents
    auto percIt = itentIt->begin();

    // For each line in the cluster fill in the vectors for the
    // gene IDs and genome IDs for the respective list element
    for ( auto it = gIdIt->begin(); it < gIdIt->end(); it++ )
    {
      int rIdx = GetGenomeIdx( *it );
      if ( rIdx != -1 )
      {
        if ( !geneMat( rIdx, colIdx ) ) geneMat( rIdx, colIdx ) = *percIt;
      }
      percIt ++;
    }

    // Move to the next list element
    itentIt ++;

    R_CheckUserInterrupt();
  }

  rownames( geneMat ) = Rcpp::wrap( genomeIds );
  return geneMat;
}

// -----------------------------------------------------------------------------
