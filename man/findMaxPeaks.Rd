% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CageData_Integration.R
\name{findMaxPeaks}
\alias{findMaxPeaks}
\title{Find max peak for each transcript,
returns as data.table, without names, but with index}
\usage{
findMaxPeaks(cageOverlaps, filteredrawCageData)
}
\arguments{
\item{cageOverlaps}{The cageOverlaps between cage and extended 5' leaders}

\item{filteredrawCageData}{The filtered raw cage-data
used to reassign 5' leaders}
}
\value{
a data.table of max peaks
}
\description{
Find max peak for each transcript,
returns as data.table, without names, but with index
}
