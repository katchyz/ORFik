#' Converts different type of files to Granges
#'
#' column 5 will be set to score
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x An data.frame from imported bed-file,
#'  to convert to GRanges
#' @param bed6 If bed6, no meta column is added
#' @return a GRanges object from bed
bedToGR <- function(x, bed6 = TRUE){

  if (!bed6) {
    gr <- GRanges(x[, 1], IRanges(x[, 2] + 1, x[, 3]))
    return(gr)
  }
  starts <- x[, 2] + 1
  ends <- x[, 3]
  gr <- GRanges(x[, 1], IRanges(starts, ends),
                strand = x[, 6])
  score(gr) <- x[, 5]
  if (ncol(x) > 6) mcols(gr) <- x[, 7:ncol(x)]
  return(gr)
}


#' Get Cage-Data From a file-path
#' @param filePath The location of the cage-file
#' @importFrom data.table fread setDF
#' @return a GRanges object
cageFromFile <- function(filePath){

  if (.Platform$OS.type == "unix") {
    if (file.exists(filePath)) {
      if (any(gsub(pattern = ".*\\.", "",
                   filePath) == c("gzip", "gz", "bgz"))) {
        rawCage <- bedToGR(setDF(
          fread(paste("gunzip -c", filePath), sep = "\t")))
      } else if (gsub(pattern = ".*\\.", "", filePath) == "bed"){
        rawCage <- bedToGR(setDF(fread(filePath, sep = "\t")))
      } else {stop("Only bed and gzip formats are supported for filePath")}
    } else {stop("Filepath specified does not name existing file.")}
  } else {
    stop("Only unix operating-systems currently support filePath,
         use cage as GRanges argument instead.")}

  message("Loaded cage-file successfully")
  return(rawCage)
}

#' Filter peak of cage-data by value
#' @param rawCage The raw cage-data
#' @param filterValue The number of counts(score) to filter on
#'  for a tss to pass as hit
#' @return the filtered Granges object
filterCage <- function(rawCage, filterValue = 1){

  if (is.null(rawCage$score)) stop("Found no score column in the cageData",
                                   "-file, bed standard is column 5.")
  filteredCage <- rawCage[rawCage$score > filterValue, ] #filter on score
  return(filteredCage)
}

#' Match seqnames
#'
#' Check that seqnames of fiveUTRs and cage uses same standard, i.g chr1 vs 1.
#' @param filteredCage Cage-data to check seqnames in
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @return filteredCage with matched seqnames convention
matchSeqnames <- function(filteredCage, fiveUTRs){
  fiveSeqlevels <- seqlevels(unlist(fiveUTRs, use.names = FALSE))
  cageSeqlevels <- seqlevels(filteredCage)
  if (length(grep(pattern = "chr", fiveSeqlevels)) > 0 &&
     length(grep(pattern = "chr", cageSeqlevels)) == 0) {
    message("seqnames use different chromosome naming conventions,",
            " trying to fix them")
    # chr1, chr2, not chrX, chrY etc. ->
    regexNormalChr <- '(^[a-zA-Z])*([0-9]+)'
    normalChr <- paste0("chr", grep(regexNormalChr,
                                    cageSeqlevels, value = TRUE))
    normalChrInd <- grep(regexNormalChr, cageSeqlevels)

    for(i in normalChrInd){
      seqlevels(filteredCage)[i] <- sub(regexNormalChr,
                                        normalChr[i], cageSeqlevels[i])
    }
  }
  if (length(grep("chrY", fiveSeqlevels)) == 0 &&
        length(grep("chrY", cageSeqlevels)) != 0)
    seqlevels(filteredCage) <- sub("chrY", "Y", seqlevels(filteredCage))
  if (length(grep("chrX", fiveSeqlevels)) == 0 &&
        length(grep("chrX", cageSeqlevels)) != 0)
    seqlevels(filteredCage) <- sub("chrX", "X", seqlevels(filteredCage))
  if (length(grep("chrM", fiveSeqlevels)) == 0 &&
        length(grep("chrM", cageSeqlevels)) != 0)
    seqlevels(filteredCage) <- sub("chrM", "MT", seqlevels(filteredCage))

  return(filteredCage)
}

#' Extends leaders downstream
#'
#' When reassigning Transcript start sites,
#'  often you want to add downstream too.
#' This is a simple way to do that
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cds If you want to extend 5' leaders downstream,
#'  to catch uorfs going into cds, include it.
#' @importFrom S4Vectors pc
#' @return a GRangesList of cds exons added to ends
addFirstCdsOnLeaderEnds <- function(fiveUTRs, cds){

  if (length(cds) == 0) {
    warning("cds is empty, returning without using it.")
    return(fiveUTRs)
  }
  if (is.null(names(cds))) {
    warning("cds have no names, returning without using it.")
    return(fiveUTRs)
  }
  matchingNames <- names(fiveUTRs) %in% names(cds)
  areValidNames <- (sum(matchingNames) - length(names(fiveUTRs))) != 0
  if (areValidNames) {
    warning("not all cds names matches fiveUTRs names,
            returning without using cds.")
    return(fiveUTRs)
  }
  # get only the ones we need
  # select first in every, they must be sorted!
  firstExons <- firstExonPerGroup(cds[names(fiveUTRs)])
  gr <- unlist(firstExons, use.names = FALSE)
  # fix mcols of cds, so that pc() will work
  mcols(gr) <- as.data.frame(mcols(unlist(fiveUTRs,
                use.names = FALSE)))[1:length(gr),]

  grl <- relist(gr, firstExons)
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl)

  return(reduce(fiveUTRsWithCdsExons))
}

#' Extend first exon of each transcript with length specified
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param extension The number of basses upstream to add on transcripts
#' @return granges object of first exons
extendsTSSexons <- function(fiveUTRs, extension = 1000){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  if (is.null(fiveAsgr$exon_rank))
    stop("fiveUTRs need column called exon_rank,",
         " see ?makeTranscriptDbFromGFF")
  ##TODO: I should make this optional,
  #so that we dont need the exon_rank column..

  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  posIDs <- strandBool(firstExons)
  promo <- promoters(firstExons, upstream = extension)

  start(firstExons[posIDs]) <- start(promo[posIDs])
  end(firstExons[!posIDs]) <- end(promo[!posIDs])

  return(firstExons)
}

#' Find max peak for each transcript,
#' returns as data.table, without names, but with index
#' @param cageOverlaps The cageOverlaps between cage and extended 5' leaders
#' @param filteredrawCageData The filtered raw cage-data
#'  used to reassign 5' leaders
#' @importFrom data.table as.data.table
#' @return a data.table of max peaks
findMaxPeaks <- function(cageOverlaps, filteredrawCageData){

  dt <- as.data.table(filteredrawCageData)
  dt <- dt[from(cageOverlaps)]
  dt$to <- to(cageOverlaps)

  maxPeaks <- dt[, max(score), by = to]
  names(maxPeaks) <- c("to", "score")
  maxPeaks <-  merge(maxPeaks, dt)

  return(maxPeaks[!duplicated(maxPeaks$to)])
}



#' Finds max peaks per trancsript from reads in the cagefile
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cageData The location of the cage-file
#' @param extension The number of basses upstream to add on transcripts
#' @return a Hits object
findNewTSS <- function(fiveUTRs, cageData, extension){

  shiftedfiveUTRs <- extendsTSSexons(fiveUTRs, extension)
  cageOverlaps <- findOverlaps(query = cageData, subject = shiftedfiveUTRs)
  maxPeakPosition <- findMaxPeaks(cageOverlaps, cageData)
  return(maxPeakPosition)
}

#' After all transcript start sites have been updated from cage,
#'  put grangeslist back together
#' @param firstExons The first exon of every transcript from 5' leaders
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @return a GRangesList
assignFirstExons <- function(firstExons, fiveUTRs){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  fiveAsgr[fiveAsgr$exon_rank == 1] <- firstExons
  return(relist(fiveAsgr, fiveUTRs))
}

#' add cage max peaks as new transcript start sites for each 5' leader
#' (*) strands are not supported, since direction must be known.
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param maxPeakPosition The max peak for each 5' leader found by cage
#' @return a GRanges object of first exons
addNewTSSOnLeaders <- function(fiveUTRs, maxPeakPosition){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  maxPeakPosition$names <- names(firstExons[maxPeakPosition$to])
  posIDs <- maxPeakPosition$to[maxPeakPosition$strand == "+"]
  minIDs <- maxPeakPosition$to[maxPeakPosition$strand == "-"]

  firstExons[posIDs] <- resize(
    firstExons[posIDs],
    width = end(firstExons[posIDs]) - maxPeakPosition$start[
      maxPeakPosition$strand == "+"] + 1,
    fix = "end")
  firstExons[minIDs] <- resize(
    firstExons[minIDs],
    width = maxPeakPosition$end[
      maxPeakPosition$strand == "-"] - start(firstExons[minIDs]) + 1,
    fix = "end")
  # Might need an chromosome boundary here? current test show no need.

  return(firstExons)
}


#' Reassign all Transcript start sites(tss')
#'
#' Given a GRangesList of 5' utrs or transcripts, reassign the start
#'  postitions. A max peak is defined as new tss if it is within boundary
#'  of 5'leader range, specified by extension
#'  A max peak must also be greater that the cage peak cutoff specified
#'  in filterValue
#'  The new TSS will then be the position of the cage read,
#'  with highest read count in the interval.
#'
#' This method can also be used with other ngs methods, if i.g. if you have a
#' list of possible cds starts (with scores) and cds'. It would reassign the
#' cds starts.
#' @param fiveUTRs The 5' leader or transcript
#'  sequences as GRangesList
#' @param cage Either a  filePath for cage-file, or already loaded
#'  R-object as GRanges
#' @param extension The maximum number of basses upstream the
#'  cage-peak can be from original tss
#' @param filterValue The number of counts(score) to filter on,
#'  for a tss to pass as hit
#' @param cds If you want to extend 5' leaders downstream,
#'  to catch upstream ORFs going into cds, include it.
#' @examples
#' library(GenomicFeatures)
#' samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                          package = "GenomicFeatures")
#' txdb <- loadDb(samplefile)
#' fiveUTRs <- fiveUTRsByTranscript(txdb) # <- extract only 5' leaders
#' cds <- cdsBy(txdb,"tx",use.names = TRUE)[1:length(fiveUTRs)]
#' names(cds) <- names(fiveUTRs)
#' # make cage from promoter of leaders (just for example)
#' cage <- GRanges(seqnames = as.character(seqnames(fiveUTRs)[1:2]),
#'    ranges =  IRanges(as.integer(start(fiveUTRs)[1 : 2] - 500) ,
#'                                as.integer(start(fiveUTRs)[1 : 2])),
#'              strand = as.character(strand(fiveUTRs)[1 : 2]), score = c(5, 10))
#' # now run it:
#' test_result <- reassignTSSbyCage(fiveUTRs[1], cage = cage,
#'  cds = cds)
#' @export
#' @return a GRangesList
reassignTSSbyCage <- function(fiveUTRs, cage, extension = 1000,
                              filterValue = 1, cds = NULL){

  ###Read in cage files
  validGRL(class(fiveUTRs), "fiveUTRs")
  if (!is.null(cds) & class(cds) != "GRangesList")
    stop("cds must be type GRangesList!")
  if (!is.numeric(extension)) stop("extension must be numeric!")
  if (!is.numeric(filterValue)) stop("filterValue must be numeric!")
  if (is.null(cage)) stop("Cage can not be NULL")

  if (class(cage) == "character") { # <- get cage file
    filteredCage <- filterCage(cageFromFile(cage),
                               filterValue) # get the cage data
  } else if (class(cage) == "GRanges") {
    filteredCage <- filterCage(cage, filterValue)
  } else {
    stop("Cage-file must be either a valid character filepath or GRanges object")
  }
  # check that seqnames match
  filteredCage <- matchSeqnames(filteredCage, fiveUTRs)

  maxPeakPosition <- findNewTSS(fiveUTRs, filteredCage, extension)
  fiveUTRs <- assignFirstExons(addNewTSSOnLeaders(fiveUTRs,
                                                  maxPeakPosition), fiveUTRs)
  if(!is.null(cds)) fiveUTRs <- addFirstCdsOnLeaderEnds(fiveUTRs, cds)

  return(fiveUTRs)
}
