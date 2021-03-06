% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift_footprints.R
\name{shift_footprints}
\alias{shift_footprints}
\title{Shifts footprints}
\usage{
shift_footprints(bam_input, selected_lengths, selected_shifts)
}
\arguments{
\item{bam_input}{Path to bam file. Bam file should have .bai file generated.}

\item{selected_lengths}{Numeric vector of lengths of footprints you select for shifting.}

\item{selected_shifts}{Numeric vector of shifts for coresponding selected_lengths.
eg. c(10, -10) with selected_lengths of c(31, 32) means length of 31 will be shifted left by 10.
Footprints of length 32 will be shifted right by 10.}
}
\value{
A GRanges object of shifted footprints.
}
\description{
Function shifts footprints from bam file and loads them to memory as GRanges.
Resizes them to single base in 5' end fashion,
treated as p site. Takes account for junctions in cigars. Length of footprint is saved in
'size' parameter of GRanges output. Footprints are also sorted according to their genomic
position, ready for saving e.g. as bed file.
}
\examples{
#shift_footprints()
}
