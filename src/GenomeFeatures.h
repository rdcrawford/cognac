#include "BioSeq.h"
using namespace std;

// -----------------------------------------------------------------------------
// GenomeFeatures
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------
// This class reads in the relevant columns from a gff file ad parses them
// into vectors. A unique gene ID is assigned to each entry in the gff file.
// The contig names are converted to the index of the contig in the
// fasta file.
// -----------------------------------------------------------------------------

// This class is composed a fasta file and the components of a gff file.
// One member function parses the fasta file and returns the neucleotide
// sequence of the gene
#ifndef _GENOME_FEATURES_
#define _GENOME_FEATURES_
class GenomeFeatures
{
public:

  // Default ctor
  GenomeFeatures( )
  { ; }

  GenomeFeatures( const std::string &gffPath, const std::string &genomeId ):
    gffPath( gffPath ), genomeId( genomeId )
  { ; }

  // Ctor: takes the relevant columns from the gff for parsing fasta files.
  // Other attributes remain empty. Might be a bad idea
  GenomeFeatures( vector<char> strand, vector<int> startPos,
    vector<int> endPos, vector<int> contig ):
    strand( strand ), startPos( startPos ),
    endPos( endPos ), contig( contig )
  { ; }

  // Dtor: does nothing
  ~GenomeFeatures( )
  { ; }

  // Return the number of genes
  int getNumGenes();

  // Parse the gff file
  bool parseGff( BioSeq *wgs );

  // Create the data-frame to be read into R
  Rcpp::DataFrame createGeneData();

  // Pass in the name of the gff file and parse the data
  bool setGff( std::string gffPath, const std::string &genomeId, BioSeq *wgs );

  // Return the identifier for this genome
  std::string getGenomeId();

  // Get the reference to the gene ids
  std::vector< std::string >* getGeneIdRef();

private:

  // Allow the genome class acess to the features to facilitate parsing genes
  friend class Genome;

  // The name assigned to this genome
  std::string genomeId;

  // Path to the gff file
  std::string gffPath;

  // Vectors containing the relevant columns of the gff file
  std::vector< std::string > featId;
  std::vector< std::string > description;
  std::vector< char >        strand;
  std::vector< int >         startPos;
  std::vector< int >         endPos;
  std::vector< int >         contig;

  // Takes a line from a gff file and assigns the attributes of this
  // gene to the component vectors
  bool parseGffEntry( const std::string &line, BioSeq *wgs );

  // Parse the gene annotation from the gff file to get the relevant
  // attributes of the gene
  std::string getDescription( std::string annotation );

  int numGenes = 0;
};
#endif
