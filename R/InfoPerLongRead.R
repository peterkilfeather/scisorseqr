#' Combine all information per read to a flat file
#' @aliases LongReadInfo
#' @description Function to concatenate all the information per read, i.e
#' gene-name, cellular barcode, UMI, cell-type information, and isoform
#' information into one file. Also outputs basic stats such as number of
#' reads/ genes/ UMIs per cellular barcode
#' @seealso \code{\link{MapAndFilter}}
#' @seealso \code{\link{GetBarcodes}}
#' @param barcodeOutputFile .csv file containing barcode and cell-type information
#' per read from the output of \code{\link{GetBarcodes}}
#' @param mapAndFilterOut output directory of the mapping function. If full-length
#' reads have been filtered using CAGE and PolyA site peaks, then it defaults
#' to that output, else it uses the canonically spliced full-length reads
#' @param minTimesIsoObserve minimum number of times an isoform is observed in the
#' dataset. Defaults to 5
#' @export
InfoPerLongRead <- function(barcodeOutputFile, mapAndFilterOut, minTimesIsoObserve = 5) {
  longReadInfo_sh <- system.file("bash", "longReadInfo.sh", package = "scisorseqr")

  longReadInfoFolder <- "LongReadInfo/"
  dir.create(longReadInfoFolder)

  if(file.exists(file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz'))){
    stretchesFile <- file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz')
    incompFile <- file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/incompleteStretches.gz')
  } else {
    stretchesFile <- file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/stretches.gz')
    incompFile <- ""
  }

  geneFile <- file.path(mapAndFilterOut,'mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz')
  exonFile <- file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/exonStretches.gz')

  longReadInfoComm <- paste(longReadInfo_sh, barcodeOutputFile, geneFile, stretchesFile, 
          exonFile, minTimesIsoObserve, incompFile)
  system(longReadInfoComm)

}
