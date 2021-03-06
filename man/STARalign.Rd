% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STARalign.R
\name{STARalign}
\alias{STARalign}
\alias{STAR}
\title{Align fastq.gz files with STARlong}
\usage{
STARalign('FastqFiles/','~/STARlong','~/starIndex_gencode10_sequins/',32)
}
\arguments{
\item{fqFolder}{fastq.gz files from a single sample or replicate,
or barcoded output}

\item{starProgPath}{path to STARlong aligner}

\item{refGenome}{location of reference genome. Defaults to mouse mm10}

\item{numThreads}{number of parallel threads to use, Defaults to 8}
}
\value{
STARoutput Folder containing star output files and reports
}
\description{
This function is a wrapper for STARlong
with stringent parameters optimized for PacBio sequenced
data. It performs sam to bam conversion to reduce
space, and .bam output is necessary for downstream analysis
}
\seealso{
\code{\link{MapAndFilter}}
}
