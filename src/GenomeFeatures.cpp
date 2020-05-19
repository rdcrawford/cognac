// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <fstream>
#include <algorithm>
#include "GenomeFeatures.h"
#include "BioSeq.h"
using namespace std;

// -----------------------------------------------------------------------------
// GenomeFeatures
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------

bool GenomeFeatures::parseGfs( BioSeq *wgs )
{
  ifstream    gff;
  std::string line;

  // Open the input stream
  gff.open( gfPath.c_str() );
  if ( !gff.is_open() ) return false;

  // Read in the file line by line
  while ( getline( gff, line ) )
  {
    // Some gff files have the wgs appended to the end. If this has
    // the wgs, there are no more genes to parse.
    if ( line[0]  == '#' )
    {
      if ( line.find( "FASTA" ) != string::npos ) break;
    }
    else
    {
      // Add the gene described in this line to the
      parseGffEntry( line, wgs );
    }
  }
  gff.close();
  return true;
}

int GenomeFeatures::getNumGenes()
{
  return numGenes;
}
std::string GenomeFeatures::getDescription( std::string attributes )
{
  // Find the position of the delimiter in the target string and subset the
  // string
  unsigned int start = attributes.find( "Name=" );

  // If the name was not found, look for a product ID
  if ( start == string::npos )
  {
    start = attributes.find( "product=" );
    if ( start == string::npos ) return "";
    start += 8;
  }
  else
  {
    start += 5;
  }

  // Find the end of the description
  unsigned int end = attributes.find( ";", start );
  if ( end == string::npos ) end = attributes.length();

  // Get the substring corresponding to the descrition
  string description = attributes.substr( start, end - start );

  // Check if there is a note in the attributes
  unsigned int noteStart = attributes.find( "Note=" );
  if ( noteStart == string::npos ) return description;
  else noteStart += 5;

  // Find the end of the note
  end = attributes.find( ";", noteStart );
  if ( end == string::npos ) end = attributes.length();

  return description + ' ' + attributes.substr( noteStart, end - noteStart );
}


bool GenomeFeatures::parseGffEntry( const std::string &line, BioSeq *wgs )
{
  stringstream ss( line );
  string       fContig;
  string       method;
  string       type;
  string       attributes;
  string         fStrand;
  int          fStart;
  int          fEnd;
  char         whoKnows;
  int          noIdea;
  int          contIdx;

  // Get the attributes of the line
  ss >> fContig;
  ss >> method;
  ss >> type;
  ss >> fStart;
  ss >> fEnd;
  ss >> whoKnows;
  ss >> fStrand;
  ss >> noIdea;

  // The final entry is the remainder of the line
  attributes = ss.str();

  // If this contig was not able to be assigned, dont continue
  if ( !wgs->getSeqIndex( fContig, contIdx ) ) return false;

  // IF this is not a coding sequence, I dont care about it right now
  if ( type != "CDS" ) return false;

  // Increment the number of genes to keep track of how many genes have
  // been added
  numGenes ++;

  // Add the atributes of the genes to the vectors
  featId.push_back( "fig|" + genomeId + ".peg." + to_string( numGenes ) );
  description.push_back( getDescription( attributes ) );
  contig.push_back( contIdx + 1 );
  startPos.push_back( fStart );
  endPos.push_back( fEnd );
  strand.push_back( fStrand );

  // Finished parsing this entry. Return true to indcate that this gene
  // was added sucessfully
  return true;
}
  // Get the reference to the gene ids
std::vector< std::string >* GenomeFeatures::getGeneIdRef()
{
  return & featId;
}

Rcpp::DataFrame GenomeFeatures::createGeneData()
{
  Rcpp::DataFrame geneData = Rcpp::DataFrame::create(
    Rcpp::_["featId"]           = featId,
    Rcpp::_["description"]      = description,
    Rcpp::_["contig"]           = contig,
    Rcpp::_["startPos"]         = startPos,
    Rcpp::_["endPos"]           = endPos,
    Rcpp::_["strand"]           = strand,
    Rcpp::_["stringsAsFactors"] = false
    );

  return geneData;
}

// Return the identifier for this genome
std::string GenomeFeatures::getGenomeId()
{
  return genomeId;
}

bool GenomeFeatures::setGff(
  std::string gffPath, const std::string &genomeId, BioSeq *wgs
  )
{
  this->gfPath   = gfPath;
  this->genomeId = genomeId;
  return parseGfs( wgs );
}

// -----------------------------------------------------------------------------
