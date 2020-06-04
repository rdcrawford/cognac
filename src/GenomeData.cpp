// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include "GenomeData.h"
#include "Genome.h"
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
GenomeData::GenomeData(
  const std::vector< std::string > &gffPaths,
  const std::vector< std::string > &faPaths,
  const std::vector< std::string > &genomeIds
  )
{
  // Initialize the vector of genome class objects
  genomeData.reserve( gffPaths.size() );
  for ( unsigned int i = 0; i < gffPaths.size(); i++ )
    genomeData.push_back( Genome( faPaths[i], gffPaths[i], genomeIds[i] ) );

  // Parse the data
  parseGenomeData();

  // Count the total number of genes that are in this dataset
  countNumGenes();

  // Concatenate the amino acid sequences into a single vector
  concatenateGeneVecs();
}

// Parse the genome data using tbb parallel_for
void GenomeData::parseGenomeData()
{
  // Use paralle for each to iterate over the vector of genomes and
  // parse the
  tbb::parallel_for_each(
    genomeData.begin(),
    genomeData.end(),
    [&] ( Genome &g )
  {
    // Read in and parse the fasta and gff files
    g.parseGenome();

    // Translate the amino acid sequences
    if ( !g.translateSeqs() )
      Rcpp::stop( "Translating genes failed for " + g.getGenomeId() );
  });
}


// Return the anino acid sequences for all of the input genomes
// Create the concatenated gene sequences with the gene ids
void GenomeData::concatenateGeneVecs()
{
  // Allocate space for the total number of genes. This is important
  // because the vector will be very large
  aaSeqs.reserve( numGenes );
  geneIds.reserve( numGenes );

  // For each genome Retrive the coding sequences
  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
  {
    // Get the reference to the vector of amino acid sequences
    auto seqRef   = it->getAaSeqRef();
    auto seqIdRef = it->getGeneIdRef();

    // Concatenate the vectors
    aaSeqs.insert( aaSeqs.end(), seqRef->begin(), seqRef->end() );
    geneIds.insert( geneIds.end(), seqIdRef->begin(), seqIdRef->end() );
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
  // Initialize the output file stream and open for writing
  std::ofstream ofs;
  ofs.open( faaPath.c_str() );
  for ( unsigned int i = 0; i < aaSeqs.size(); i++ )
  {
    ofs << ">" << geneIds[i] << endl << aaSeqs[i] << endl;
  }
}

// Count the total number of genes that are in this dataset
void GenomeData::countNumGenes()
{
  // Initialize the number of genes to zero
  this->numGenes = 0;

  // Iterate over the vector of genomes and get the count of the
  // number of genes in each genome
  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
    numGenes += it->getNumGenes();
}

// Create a list of dataframes containing the parsed genome
// features
Rcpp::List GenomeData::createGeneDataFrames()
{
  // Initialize the list to output
  Rcpp::List gfList( genomeData.size() );

  // For each genome, create a data frame with the features
  for ( unsigned int i = 0; i < genomeData.size(); i++ )
    gfList[i] = genomeData[i].createGeneData();

  return gfList;
}

// -----------------------------------------------------------------------------
