% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotateGenome.R
\name{AnnotateGenome}
\alias{AnnotateGenome}
\title{Create a gff file for a bacterial genome}
\usage{
AnnotateGenome(faPath, outDir, genomeName, keepGto)
}
\arguments{
\item{faPath}{Path to the fasta file containing the whole genome sequence
to be annotated.}

\item{outDir}{Directory where the the gff and gto are to be written.}

\item{genomeName}{Optional name for the gff file. By default, the extension
is removed from the fasta file and the remainder is used.}

\item{keepGto}{Optional bool to keep the "genome type object" created
by RAST for further processing (see the rast documentation). False by 
default}
}
\value{
The path to the gff file
}
\description{
This function uses the RAST command line interface to annotate the 
  genome and generate a GFF-3 file. First, a genome type object is
  generated with RAST. Then the genome is processed, and the 
  features are exported to a GFF file.
}
