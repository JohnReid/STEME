
STEME Options
=============




.. _option-group-multiple-motifs:

Multiple motifs
---------------


Control how to find more than one motif.


.. _option-num-motifs:

``--num-motifs`` *(default=1)*: Number of motifs to look for.

.. _option-prediction-z-threshold:

``--prediction-Z-threshold`` *(default=0.3)*: The threshold on Z used to erase instances of motifs. The lower this is, the more instances will be erased.




.. _option-group-output:

Output
------


Control the output location, format, writing logos, etc...


.. _option-output-dir:

``--output-dir`` *(default=output)*: Output directory.

.. _option-meme-like-output:

``--meme-like-output`` *(default=meme.txt)*: Produce MEME-like output so that it can be parsed by downstream tools.

.. _option-html-output:

``--html-output`` *(default=STEME.html)*: Produce HTML output.

.. _option-print-sites:

``--print-sites`` *(default=False)*: Write a file containing the sites that were used to make the motif.

.. _option-dont-write-logos:

``--dont-write-logos`` *(default=False)*: Don't write logos for motifs.

.. _option-write-em-logos:

``--write-em-logos`` *(default=False)*: Write logos for motifs during EM algorithm.

.. _option-write-em-stats:

``--write-em-stats`` *(default=False)*: Write statistics for EM algorithm.

.. _option-tomtom:

``--tomtom`` *(default=[])*: Run TOMTOM tool from the the MEME suite on the motifs using the specified motif databases.

.. _option-max-seqs-to-write:

``--max-seqs-to-write`` *(default=1000)*: Maximum number of sequences to write information about in output.

.. _option-max-sites-to-write:

``--max-sites-to-write`` *(default=1000)*: Maximum number of sites to write information about in output.




.. _option-group-background-model:

Background model
----------------


Control the background model.


.. _option-bg-model-order:

``--bg-model-order`` *(default=2)*: Order of the background Markov model.

.. _option-bg-fasta-file:

``--bg-fasta-file`` : If specified, STEME builds its background model from the sequences in this file rather than from the input sequences.

.. _option-back-dist-prior:

``--back-dist-prior`` *(default=1.0)*: Pseudo-counts for Markov background model.




.. _option-group-start-finding:

Start finding
-------------


Control how starts are found.


.. _option-max-start-finding-time:

``--max-start-finding-time`` *(default=0.0)*: How many seconds to dedicate to finding starts for each motif. If not given, STEME will look at each possible start (can be slow).

.. _option-min-num-sites:

``--min-num-sites`` : Minimum number of sites. Defaults to # sequences / 10.

.. _option-max-num-sites:

``--max-num-sites`` : Maximum number of sites. Defaults to 50% more than # sequences.

.. _option-width:

``--width`` *(default=[])*: If specified, search for motifs of this width (can specify more than one).

.. _option-min-w:

``--min-w`` *(default=6)*: Minimum width of motif to find.

.. _option-max-w:

``--max-w`` *(default=14)*: Maximum width of motif to find.

.. _option-starts-per-motif:

``--starts-per-motif`` *(default=4)*: Number of starts to find per motif.

.. _option-use-seed:

``--use-seed`` : If specified, only use this seed as a start.

.. _option-starts-seed-pseudo-counts:

``--starts-seed-pseudo-counts`` *(default=0.5)*: Pseudo counts with which to smooth possible starting seeds.

.. _option-starts-speed-up:

``--starts-speed-up`` *(default=0)*: Speed up the start finding by ignoring so many potential starting points.

.. _option-candidate-starts-factor:

``--candidate-starts-factor`` *(default=1.41421356237)*: The factor for the geometric progression that determines which numbers of sites to try when start finding.




.. _option-group-em:

EM
--


Control the behaviour of the Expectation Maximization algorithm.


.. _option-max-iters:

``--max-iters`` *(default=1000)*: Maximum number of iterations for EM algorithm.

.. _option-dont-discretize:

``--dont-discretize`` *(default=False)*: Don't run discretisation after EM.

.. _option-convergence-distance:

``--convergence-distance`` *(default=1e-05)*: Threshold between successive iterations at which to stop EM.

.. _option-wnsites:

``--wnsites`` *(default=0.8)*: Weight on number of sites. Used when updating lambda in EM algorithm.

.. _option-em-seed-pseudo-counts:

``--em-seed-pseudo-counts`` *(default=0.01)*: Pseudo counts for motif model in EM algorithm.

.. _option-epsilon:

``--epsilon`` *(default=0.4)*: Allowed error in motif probabilities for EM algorithm.

