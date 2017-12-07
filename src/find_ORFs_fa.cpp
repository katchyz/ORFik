// Find many orfs faster by putting all the data into Cpp directly
// Cpp format: Webkit: 80 character max line

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Rcpp.h>
#include "findORFsHelpers.h"

using vi = std::vector<int>;
using namespace Rcpp;

Function GRangesC("GRanges", Environment::namespace_env("GenomicRanges"));
Function IRangesC("IRanges", Environment::namespace_env("IRanges"));

char complement(char n)
{
  switch (n) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  case 'N':
    return 'N';
  case 'a':
    return 't';
  case 't':
    return 'a';
  case 'g':
    return 'c';
  case 'c':
    return 'g';
  }
  assert(false);
  return ' ';
}

// ! Only for procaryotic circular genomes
// Takes a fasta file and finds all orfs,
// including orfs that span from before stop to over start position of genome
// [[Rcpp::export]]
S4 findORFs_procaryote(std::string file,
                       std::string startCodon,
                       std::string stopCodon,
                       bool longestORF,
                       int minimumLength)
{
  std::vector<int> all_orfs;
  std::vector<int> Seqnames;
  std::vector<std::string> strands;
  std::ifstream in(file.c_str());
  in.get(); // remove first '>'
  std::string rec;
  int n = 0;
  while (getline(in, rec, '>')) { // For each chromosome
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0, newLineLoc);
    std::string fastaSeq = rec.substr(newLineLoc + 1,
                                      rec.length() - newLineLoc - 2);
    // get all orfs for start to stop
    std::vector<int> ORFdef = orfs_as_vector(fastaSeq, startCodon,
                                             stopCodon, longestORF,
                                             minimumLength);
    all_orfs.insert(all_orfs.end(), ORFdef.begin(), ORFdef.end());
    Seqnames.insert(Seqnames.end(), ORFdef.size() / 2, n++);
    strands.insert(strands.end(), ORFdef.size() / 2, "+");

    // Now do stop/start boundary, +/-3000, only keep the ones
    // who are overlapping start/stop boundary
    std::string startStopBoundary = fastaSeq.substr(
      fastaSeq.length() - 3000, 3000);
    startStopBoundary.append(fastaSeq.substr(0, 3000));
    std::vector<int> ORFdefBoundary = orfs_as_vector(startStopBoundary,
                                                     startCodon,
                                                     stopCodon, longestORF,
                                                     minimumLength);
    // now filter out wrong ones, the ones that does not contain point
    // 3000 in range start/stop
    std::vector<int> ORFdefOverlap; // <- vector only for valid ones
    int rescaleStart = fastaSeq.length() - 3000;
    int rescaleStop = 3000;
    for (size_t i = 0; i < ORFdefBoundary.size() / 2; i++) {
      if (ORFdefBoundary[2 * i] < 3000) { // start
        if (ORFdefBoundary[2 * i + 1] > 3000) { // stop
          ORFdefOverlap.push_back(ORFdefBoundary[2 * i] + rescaleStart);
          ORFdefOverlap.push_back(ORFdefBoundary[2 * i + 1] - rescaleStop);
        }
      }
    }
    all_orfs.insert(all_orfs.end(), ORFdefOverlap.begin(),
                    ORFdefOverlap.end());
    Seqnames.insert(Seqnames.end(), ORFdefOverlap.size() / 2, n);
    strands.insert(strands.end(), ORFdefOverlap.size() / 2, "+");

    // now do on reverse compliment
    std::reverse(fastaSeq.begin(), fastaSeq.end());
    std::transform(fastaSeq.begin(), fastaSeq.end(), fastaSeq.begin(),
                   complement);

    ORFdef = orfs_as_vector(fastaSeq, startCodon,
                            stopCodon, longestORF,
                            minimumLength); // <-get all orfs for start to stop
    all_orfs.insert(all_orfs.end(), ORFdef.begin(), ORFdef.end());
    Seqnames.insert(Seqnames.end(), ORFdef.size() / 2, n);
    strands.insert(strands.end(), ORFdef.size() / 2, "-");

    // Now do stop/start boundary, +/-3000, only keep the ones
    // who are overlapping start/stop boundary
    startStopBoundary = fastaSeq.substr(fastaSeq.length() - 3000, 3000);
    startStopBoundary.append(fastaSeq.substr(0, 3000));
    ORFdefBoundary = orfs_as_vector(startStopBoundary, startCodon,
                                    stopCodon, longestORF,
                                    minimumLength);

    // now filter out wrong ones, the ones that does
    // not contain point 3000 in range start/stop
    std::vector<int> ORFdefOverlapMin; // <- vector only for valid ones

    for (size_t i = 0; i < ORFdefBoundary.size() / 2; i++) {
      if (ORFdefBoundary[2 * i] < 3000) { // start
        if (ORFdefBoundary[2 * i + 1] > 3000) { // stop
          ORFdefOverlapMin.push_back(ORFdefBoundary[2 * i] + rescaleStart);
          ORFdefOverlapMin.push_back(ORFdefBoundary[2 * i + 1] - rescaleStop);
        }
      }
    }
    all_orfs.insert(all_orfs.end(),
                    ORFdefOverlapMin.begin(), ORFdefOverlapMin.end());
    Seqnames.insert(Seqnames.end(), ORFdefOverlapMin.size() / 2, n);
    strands.insert(strands.end(), ORFdefOverlapMin.size() / 2, "-");
  }

  // all_orfs is an interlaced vector. We de-interlace it into two vectors.
  std::vector<vi> result_value(2);
  result_value[0].resize(all_orfs.size() / 2);
  result_value[1].resize(all_orfs.size() / 2);
  for (size_t i = 0; i < all_orfs.size() / 2; i++) {
    result_value[0][i] = all_orfs[2 * i];
    result_value[1][i] = all_orfs[2 * i + 1];
  }

  return (GRangesC(wrap(Seqnames),
                   IRangesC(wrap(result_value[0]),
                            wrap(result_value[1])),
                            wrap(strands)));
}
