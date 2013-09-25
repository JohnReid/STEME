/** Copyright John Reid 2011
 *
 * \file
 * \brief A convenient place to keep our global TODO list so that doxygen picks it up.
 *
 * \todo Check large objects are not copied around
 *
 * \todo Check performance after PSP has been fully implemented
 *
 * \todo Make as easy to install as possible. Provide Windows installer. Make 3rd party dependencies minimal
 *
 * \todo Test to see if storing positions as (seq, offset) pairs is more efficient than as global positions
 *
 * \todo Profile STEME without any inlining to try to catch any unnecessary copies
 *
 * \todo Check if finding best w-mers can be speeded up with heap storage
 *
 * \todo Check if find_best_wmers is quicker with list storage instead of vector.
 *
 * \todo When searching for starts for 2nd, 3rd, 4th, motif etc... see if start information from previous
 * motifs can be used.
 *
 * \todo Parallelise start finding
 *
 * \todo Perhaps use LLR to score starts and initialise p-value tables concurrently in separate thread?
 *
 * \todo Check E-values are calculated correctly. Large data sets (e.g. MEIS1) have given both +inf and -inf
 * E-values.
 *
 * \todo Check E-values are sorted correctly.
 *
 * \todo Add test with unknown values in sequences.
 *
 * \todo Ability to save and load Markov models to/from files.
 *
 * \todo Add tests to make sure we are discretising correctly. I.e. picking correct numbers of sites.
 *
 * \todo Add tests to make sure we are comparing the results of EM correctly.
 *
 * \todo Add tests to make sure we are calculating E-values correctly.
 *
 * \todo Make sure we do a fair comparison between MEME and FAST p-value calculations.
 *
 * \todo Investigate whether sparse data structure would speed up shifted Hirji implementation.
 *
 * \todo We can often find the same motif more than once. Perhaps we should erase more occurrences of each motif
 * than we currently do. This is especially true when the number of sites is maximal, i.e. the algorithm could not
 * have chosen more even if it would have been more significant.
 *
 * \todo Revamp start finding for multiple motifs. It is surely very wasteful to go through the whole procedure
 * every time.
 *
 * \todo Add mosaic background models.
 *
 * \todo Test if high-order background models are good at discriminative motif finding.
 *
 * \todo Compare to other motif finders that are suitable for large data sets:
 * ChipMunk, DREME, MoAn, Trawler, DRIM, Amadeus and others?
 *
 * \todo If using a separate fasta file for the background model, still use the frequencies from the input
 * sequences to calculate the significance.
 *
 * \todo Fix memory_error on large input sequences.
 */


/**
 * Dummy struct to store already done TODO items.
 *
 * \todo Change to use number of W-mers without 'N's for significance calculations.
 *
 * \todo Investigate valgrind errors in FAST calculations
 *
 * \todo Check memory usage. In particular check that EMDescender does not hang on to large Z array when
 * no longer needed.
 *
 * \todo Reorganize python code so that each motif has its own output directory to avoid clashes
 *
 * \todo Refactor EM_descender, find_best_w_mers, start_finder out of stem class. Use TSpec template parameters
 * as in seqan to parameterise.
 *
 * \todo Add support for finding more than one motif.
 *
 * \todo Add test(s) for finding more than one motif.
 *
 * \todo Copy all input data into one string so we can use global positions to index the W-mers. This
 * should be quicker.
 *
 * \todo Improve find best w-mer efficiency by preferentially descending tree to find best w-mers first
 *
 * \todo Rewrite python stempy module. There is way too much code that could be refactored...
 *
 * \todo Remove old unused code from project. There is a lot of clutter.
 *
 * \todo Improve how background models are calculated. Do not need to scan all sequences when we already have
 * suffix array.
 *
 */
struct AlreadyDone { };
