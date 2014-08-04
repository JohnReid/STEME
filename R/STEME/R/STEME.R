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


#' Scatter plot the scores for two motifs on the same sites
#'
#' @param hits The hits member of a steme.scan
#' @param motif.1 The first motif to compare scores for
#' @param motif.2 The second motif to compare scores for
#' @param offset.1 Offset to move position of first motif hits
#' @param same.strand Match hits on the same strand?
plot.site.scores <- function(hits, motif.x, motif.y, offset.x=0, same.strand=TRUE) {
    stopifnot(motif.x %in% hits$motif)
    stopifnot(motif.y %in% hits$motif)
    hits.x <- hits %>% filter(motif == motif.x)
    # print(head(hits.x))
    # print(dim(hits.x))
    hits.y <- hits %>% filter(motif == motif.y)
    # print(dim(hits.y))
    hits.both <- (
        inner_join(hits.x, hits.y, by="seqidx")
        %>% mutate(position.diff=ifelse(strand.x == '+',
                                        position.y - position.x,
                                        position.x - position.y))
        %>% filter(xor(! same.strand, strand.x == strand.y),
                   position.diff == offset.x)
    )
    # colnames(hits.both)
    gp <- ggplot(hits.both, aes(x=Z.x, y=Z.y)) + geom_point()
    # sample_n(hits.both, 30)
    return(list(
        hits.both=hits.both,
        hits.x=hits.x,
        hits.y=hits.y,
        gp=gp))
}
