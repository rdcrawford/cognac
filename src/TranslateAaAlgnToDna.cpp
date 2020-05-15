#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------
// Parse Fasta Files By Position
// Ryan D. Crawford
// 04/10/2020
// -----------------------------------------------------------------------------
// This function prefoms reverse translation of an amino acid alignment to
// a DNA alignment. A parsed gff fle is input with only the genes contained
// in the alignment. The fasta file contining the whole genome sequence is
// parsed to extract the gene sequences, and the amino acid alignment is
// used to guide the placement of .
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void TranslateAaAlgnToDna(
  const Rcpp::DataFrame &gffData,
  const vector<string>  &faPath,
  const vector<int>     &genePositions,
  const string          &genomeName,
  const string          &aaAlgn,
  const string          &outputFile
  )
{
  // ---- Constant Declarations ------------------------------------------------

  const int ANNOT_ID  = 0; // Col index of annot IDs
  const int GENE_NAME = 1; // Col index of the gene names
  const int CONTIG    = 2; // Col index of the contig index
  const int LEFT_POS  = 3; // Col index of Gene start positons
  const int RIGHT_POS = 4; // Col index of Gene end positions
  const int STRAND    = 5; // Col index indicating is on the antisense strand

  // ---- Variable Initializations ---------------------------------------------

  // Initailize a genome class objects with the genom features necessary
  // to parse the fasta file.
  Genome genome( faPath, gffData[ STRAND ], gffData[ LEFT_POS ],
    gffData[ RIGHT_POS ], gffData[ CONTIG ] );

  // Read in th fasta file
  genome.parseFasta();

  ofstream ofs;           // Output file stream to write the nt alignment
  string   seq;           // String to store the current nucleotide sequence
  int      algnPos   = 0; // Current position in the aa alignment
  int      genePos   = 0; // Position in the nucleotide sequence
  int      alGeneIdx = 0; // Gene position in the amino acid alignment

  // ---- Create  vector with the end positions for each gene ------------------

  // Initialize the vector to the number of genes in the alignment
  vector< int > geneEnd( genePositions.size() );

  // For each gene find the right bound in the alignment. We need to keep track
  // of wherethe genes fall in order to correctly place the stop codons
  for (int i = 0; i < genePositions.size(); i ++)
  {
    geneEnd[ i ] = genePositions[ i ] - 1;
  }

  // ---- Reverse translate the alignment --------------------------------------

  // Open the output file stream
  ofs.open( outputFile.c_str(), ios::out | ios::app );
  if ( !fout.fail() ) Rcpp::stop( "Cannot write to file: ", outputFile );

  // Write the header
  ofs << ">" << genomeName << endl;

  // Get the first gene to enter the loop
  genome.getGeneSeq( seq );

  // Iterate along the alignment and place the codons and gaps in the nt
  // alignment in accordance with the aa alignment
  while ( algnPos < aaAlgn.size() )
  {
    // Write the codon correpsonding to thecurrent position in the alignment
    if ( aaAlgn[ algnPos] == '-' )
    {
      ofs << "---";
    }
    else
    {
      // Write the curren cdon and move the position in the gene to the
      // next codon
      ofs << seq.substr( genePos, 3 );
      genePos += 3;
    }

    // If this is the end of the gene, write the stop codon if the
    // gene is present in this genome
    if ( algnPos == geneEnd[ alGeneIdx ] )
    {
      // Increment the current gene in the alignment to advance to the
      // next end position
      alGeneIdx ++ ;

      // If we have written the entire sequece of this gene, get the next gene
      if ( genePos != 0 )
      {
        // ofs << seq.substr( genePos, 3 );
        bool isUpdated = genome.getGeneSeq( seq );
        genePos = 0;

        gCount ++;

      } else {

        // Insead of a stop codon, write a gap
        ofs << "---";
      }

    }
    algnPos ++; // Move to the aa residue in the alignment

    // Allow user to prematurely break out of the loop
    // if ( alGeneIdx == 6 ) break;
    R_CheckUserInterrupt();
  }

  // Write a new line character for the end of the alignment
  ofs << endl;
  ofs.close(); // Close the output file stream
}

// -----------------------------------------------------------------------------
