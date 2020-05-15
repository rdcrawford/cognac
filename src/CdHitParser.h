// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <fstream>
#include <algorithm>

class CdHitParser
{
public:

  // Ctor
  CdHitParser(
    const std::string &cdHitClstFile, // Path to the cd-hit results
    bool isBinary,                    // Return binary matrix
    bool removeLowFreq,               // Remove low frequency genes
    Rcpp::Environment &genePtr        // Environment to append results to
    );

  // Create a list with the gene presence matrix and the list of
  // genes to output to R
  Rcpp::List CreateOutputList( );

private:

  // Lists containing the gene and genome Ids of each gene in the cluster
  std::list< std::vector< std::string > > clustList;
  std::list< std::vector< std::string > > gIdList;
  std::list< std::vector< double > >      clustIdentList;

  std::vector<std::string> cdHitResults; // Vector to store the lines results
  std::vector<int>         headers;      // Position of te start clusters
  std::vector<int>         tailers;      // Position of the cluster end

  // Minimium number of genes in a cluster to keep
  int minGeneNum

  // Vecotr containing the names of all of the genomes in the analysis
  std::vector< std::string > genomeIds;

  // Remove genes beneath the threshold -- 2 if < 1000 genomes
  // leff than 80%  oft the number of input genomes otherwise
  bool RemoveLowFreqClusts( );

  // Create the lists of genes in each cluster
  void CreateClustList( );

  // Create a list of identites of the genes within each cluster
  void CreateItentList( );

  // Create a genome x gene matrix with the presene or absene of each cluster
  Rcpp::NumericMatrix CreateBinaryMat(  );

  // Create a genome x gene matrix with the presene or absene of each cluster
  Rcpp::NumericMatrix CreateIdentMat(  );

  // Extract the name of the gene from the cd-hit entry
  std::string GetGeneName( const std::string &inStr);

  // Extract the name of the genome from the gene identifier
  std::string GetGenomeId( const std::string &inStr );

  // Extract the percent id of the gene to the reference
  double GetGeneIdent( const std::string &inStr );

  // Input a vector of genomes and a genome name and return the
  // index of the genome name in the vector. Returns -1 if the
  // genome name is not found
  int GetGenomeIdx( const string &name );
};
