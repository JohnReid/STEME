#' Read, reshape and analyse the results of STEME PWM scans.
#'
#' The STEME package is designed to analyse the results of motif scans
#' performed by the STEME python software.
#'
#' @docType package
#' @name STEME
NULL

# library(stringr)
# library(reshape2)
# library(dplyr)
# library(ggplot2)
# library(rtracklayer)
# library(preprocessCore)
# library(org.Mm.eg.db)


#' Example STEME scan results
#'
#' A dataset containing the results of a STEME scan of 90 PWMs over
#' >1000 sequences.
#' The variables are as follows:
#'
#' \itemize{
#'   \item hits. A data frame with a row for each hit.
#'   \item seqs. A data frame with a row for each sequence.
#' }
#'
#' @format A list containing two data frames with information about the
#'   hits and the sequences
#' @name example.steme.scan
NULL


#' Read the results of a STEME scan from a directory
#'
#' @param scan.dir The directory to read the scan from. The directory should
#'   contain the files 'steme-pwm-scan.out' and 'steme-pwm-scan.seqs'
#' @export
read.steme.scan <- function(scan.dir) {
    hits.filename <- paste(scan.dir, 'steme-pwm-scan.out', sep='/')
    seqs.filename <- paste(scan.dir, 'steme-pwm-scan.seqs', sep='/')
    hits <- read.csv(
        hits.filename,
        skip=1,
        header=FALSE,
        col.names=c("motif", "Wmer", "seqidx", "position", "strand", "Z", "score", "pvalue"),
    )
    hits$seqidx <- hits$seqidx + 1  # Convert to 1-based index
    # dim(hits)
    # sapply(hits, class)
    # head(hits)
    # sample_n(hits, 10)
    seqs <- read.csv(
        seqs.filename,
        skip=1,
        header=FALSE,
        col.names=c("Length", "ID"),
        colClasses=c("integer", "factor"))
    # dim(seqs)
    # sapply(seqs, class)
    # head(seqs)
    # sample_n(seqs, 10)
    # print(nrow(seqs))
    stopifnot(max(hits$seqidx) <= nrow(seqs))  # Make sure seqidx are legal
    result <- list(hits=hits, seqs=seqs)
    class(result) <- c(class(result), "steme.scan")
    return(result)
}


#' Print the results of a STEME scan
#'
#' @param steme.scan The scan to print
print.steme.scan <- function(steme.scan) {
    print(steme.scan$hits)
    print(steme.scan$seqs)
}


#' Summarise the results of a STEME scan
#'
#' @param object The scan to summarise
#' @param ... Other arguments (ignored)
#' @export
## @examples
## data(example.steme.scan, package="STEME")
## summary(example.steme.scan)
summary.steme.scan <- function(object, ...) {
    cat("STEME scan:",
        nrow(object$hits), "hits",
        "for", length(levels(object$hits$motif)), "motifs",
        "in", nrow(object$seqs), "sequences,",
        "comprising", sum(object$seqs$Length), "base pairs,",
        sprintf('Z range: [%.3f, %.3f]\n', min(object$hits$Z), max(object$hits$Z)))
}


#' Calculate sequence-centric statistics from a steme.scan hits object
#'
#' @param hits The hits member of a steme.scan
#' @return A data frame containing sequence centric statistics
#' @export
calc.seq.centric <- function(hits) {
    return(
        hits
        %>% group_by(seqidx, motif)
        %>% summarise(
                total=length(Z),
                expected=sum(Z),
                best=max(Z),
                anywhere=1-Reduce("*", 1-Z, 1.)
            )
    )
}


#' Widen the results of a sequence-centric set of hits
#'
#' @param seq.centric The sequence centric statistics
#' @export
widen.seq.centric <- function(seq.centric) {
    seq.centric.molten <- melt(
        seq.centric,
        id.vars=c("seqidx", "motif"),
        variable.name="stat")
    return(dcast(seq.centric.molten, seqidx ~ motif + stat))
}


#' Find hits of two motifs a certain distance apart
#'
#' @param hits The hits member of a steme.scan
#' @param motif.x The first motif to compare scores for
#' @param motif.y The second motif to compare scores for
#' @param spacing Distance(s) between start of hits
#' @param same.strand Match hits on the same strand?
find.spaced.pairs <- function(hits, motif.x, motif.y, spacing=0, same.strand=TRUE) {
    stopifnot(motif.x %in% hits$motif)
    stopifnot(motif.y %in% hits$motif)
    hits.x <- hits %>% filter(motif == motif.x)
    W.x <- str_length(hits.x$Wmer[1])
    # print(head(hits.x))
    # print(dim(hits.x))
    hits.y <- hits %>% filter(motif == motif.y)
    W.y <- str_length(hits.y$Wmer[1])
    # print(dim(hits.y))
    hits.both <- (
        inner_join(hits.x, hits.y, by="seqidx")
        %>% mutate(offset.y=ifelse(strand.x == '+',
                                        position.y - position.x,
                                        position.x + W.x - position.y - W.y))
        %>% filter(xor(! same.strand, strand.x == strand.y),
                   offset.y %in% spacing)
    )
    return(list(
        hits.x = hits.x,
        hits.y = hits.y,
        W.x = W.x,
        W.y = W.y,
        hits.both = hits.both))
}
