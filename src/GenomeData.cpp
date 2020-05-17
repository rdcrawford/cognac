// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include "Genome.h"
#include "GenomeParser.h"
#include <RcppParallel.h>
#include <Rcpp.h>
#include <fstream>
#include <algorithm>
using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// GenomeData
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------

// Value ctor for inputs of gff files, fasta files, and the corresponding
// genome names
GenomeData::GenomeData(const std::vector< std::string > &gffPaths,
  const std::vector< std::string > &faPaths,
  const std::vector< std::string > &genomeIds)
{
  // Initialize the vector of genome class objects
  for ( int i = 0; i < gffPaths.size(); i++ )
    genomeData.push_back( Genome( gffPaths[i], faPaths[i], genomeIds[i] ) );

  // Parse the data
  parseGenomeData();

  // Count the total number of genes that are in this dataset
  countNumGenes();
}

// Value for ctor for gff files with the
GenomeData::GenomeData( const std::vector< std::string > &gffPaths,
  const std::vector< std::string > &genomeIds)
{
  // Initialize the vector of genome class objects
  for ( int i = 0; i < gffPaths.size(); i++ )
    genomeData.push_back( Genome( gffPaths[i], faPaths[i], genomeIds[i] ) );

  // Parse the data
  parseGenomeData();

  // Count the total number of genes that are in this dataset
  countNumGenes();
}

// Value ctor for gbk files
GenomeData::GenomeData(const std::vector< std::string > &gbkPaths,
  const std::vector< std::string > &faPaths,
  const std::vector< std::string > &genomeIds)
{
  // Parse the data: read in the genome features and genome sequences
  parseGenomeData();

  // Count the total number of genes that are in this dataset
  countNumGenes();

  // Concatenate the gene sequences
  concatenateGeneSeqs()
}

// Parse the genome data using tbb parallel_for
void GenomeData::parseGenomeData()
{
  // To use Intel TBB we first have to create a range object. In this case
  // we want to parse all genomes in the vector
  tbb::blocked_range< int > range( 0, genomeData.size() );

  // We have to create a temporary instance of the functor class for the
  // Genome parser. The 'func' object is a struct that behaves like a function
  auto func = GenomeParser( genomeData );

  // Call tbb::parallel_for, the code in the functor will be executed
  // on the availible number of threads.
  tbb::parallel_for( range, func );
}


// Return the anino acid sequences for all of the input genomes
// Create the concatenated gene sequences with the gene ids
void GenomeData::concatenateGenneVecs();
{
  // Allocate space for the total number of genes. This is important
  // because the vector will be very large
  aaSeqs.reserve( numGenes );

  // For each genome Retrive the coding sequences
  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
  {
    // Get the reference to the vector of amino acid sequences
    auto seqRef = it->getAaSeqRef();

    // Concatenate the vectors
    aaSeqs.insert( aaSeqs.end(), seqRef->begin(), seqRef->end() );
  }

  // Allocate space for the total number of genes. This is important
  // because the vector will be very large
  geneIds.reserve( numGenes );

  // For each genome Retrive the coding sequences
  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
  {
    // Get the reference to the vector of amino acid sequences
    auto geneIdRef = it->getAaSeqRef();

    // Concatenate the vectors
    geneIds.insert( aaSeqs.end(), geneIdRef->begin(), geneIdRef->end() );
  }
}

// Return the anino acid sequences for all of the input genomes
std::vector< std::string > GenomeData::getAaSeqs()
{
  return aaSeqs;
}

// Return the gene ids for all of the input genomes
std::vector< std::string > GenomeData::getGeneIds()
{
  return geneIds;
}

// Write the amino acid sequences to an faa file
void GenomeData::writeFaa( std::string faaPath )
{
  ofstresm ofs.open( faaPath.c_str() );

  for ( int i = 0; i < aaSeqs.size(); i++ )
  {
    ofs >> '>' >> geneIds[i] >> endl >> aaSeqs[i] >> endl;
  }
}

// Count the total number of genes that are in this dataset
void GenomeData::countNumGenes()
{
  numGenes = 0;

  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
    numGenes += it->getNumGenes();
}

// -----------------------------------------------------------------------------
