#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// -----------------------------------------------------------------------------
// Find Identical Genes
// Ryan D. Crawford
// 2020/01/23
// -----------------------------------------------------------------------------
// This function writes a fasta file including only the identical genes. The 
// gene id for any duplicate genes is stored in a 
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List FindIdenticalGenes(
  const Rcpp::StringVector &genes,
  const Rcpp::StringVector &geneIds
  )
{
  // This liststores the genes which are redundant. where each
  // element is the names of the identical genes 
  Rcpp::List identGeneList( 0 );
  
  // Logical vector indicating whether a gene has been classified as
  // a representitive sequence or 
  std::vector<bool> isClassified( genes.size() );
  
  // Calculate the length of the genes
  std::vector<int> geneLens( genes.size() );
  for (int i = 0; i < genes.size(); i++ ) geneLens[ i ] = genes(i).size();
  
  // Start at the first gene to enter the loop
  int geneIdx = 0;
  Rcpp::StringVector geneReps; 
  // cout << "Total gene count: " << genes.size() << endl;
  
  // Iterate over all of the genes and find any identical sequences 
  while ( geneIdx < genes.size() )
  {
    // cout << "geneIdx: " << geneIdx
    //      << " geneRep count: " << geneReps.size()
    //      << endl;
    // Initialize a vector to store the identical genes 
    Rcpp::StringVector identGeneIds( 0 );
    isClassified[ geneIdx ] = true;
    geneReps.push_back( geneIds[ geneIdx ] );
      
    // For each unclassified gene, find if this gene is identical to the
    // current subsect 
    for ( int i = geneIdx + 1; i < genes.size(); i++ )
    {
      if ( !isClassified[ i ] )
      {
        if ( geneLens[i] == geneLens[geneIdx])
        {
          if ( genes[i] == genes[geneIdx] )
          {
            identGeneIds.push_back( geneIds[i] );
            isClassified[ i ] = true;
          }
        }
      }
    }
    
    // Add the vector with the 
    identGeneList.push_back( identGeneIds );
    
    // Find the next unclassified gene
    do {
      geneIdx ++;
    } while( isClassified[ geneIdx ] && geneIdx < genes.size() );
  }
  identGeneList.names() = geneReps;
  return identGeneList;
}

// -----------------------------------------------------------------------------