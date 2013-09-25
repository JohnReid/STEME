.. index:: using

Using STEME
===========

If you have any queries regarding how to use STEME please don't hesitate to get in touch. My email address is given with the original STEME
paper_. My home page is here_.

.. _paper: http://nar.oxfordjournals.org/content/early/2011/07/23/nar.gkr574.long
.. _here: http://sysbio.mrc-bsu.cam.ac.uk/group/index.php/John_Reid



.. index:: motif finding

Quick start
-----------

To use STEME as a motif finder you can execute the following command::

  steme <input fasta sequence file>
  
If you need a FASTA file to experiment with, there is one_ included in the STEME distribution ``python/stempy/test/fasta/yst09m.fasta``. 

By default STEME will place its output in a directory named ``output``. Here you will find logos for the motifs STEME found,
a file called ``meme.txt`` that contains the details of the motifs in MEME format 
and a log file, ``STEME.log``, where information about the significance of the motifs can be found. The ``meme.txt`` 
file is compatible with the MEME format and can therefore be used with the downstream
tools in the MEME suite. 

By default, STEME will choose sensible defaults for its parameters. However 
there are probably some parameters you will want to set yourself. We just mention two here. 
You can see all the available options with ``--help`` option::

  steme --help

Choosing the background model is very important. Recently there has been a shift towards discriminative motif finders, that
is motif finders that learn a motif that can discriminate between a positive and a negative set of sequences. In a motif
finder like STEME, the background model plays the part of the negative set. We have found that using a high order Markov
model on the negative sequences can be a powerful technique for discriminative learning, for example::

  steme --bg-model-order=5 --bg-fasta-file=<negative sequences> <positive sequences>
  
In addition to choosing the background model, if your file is large you probably won't want STEME to examine every 
possible starting seed as this would take too long. You can control STEME's running time through the options 
in :ref:`controlling-the-running-time`. For example::

  steme --max-start-finding-time=<seconds> <input fasta sequence file>
  
Obviously the longer STEME is allowed to spend searching for starting seeds the better the results but in practice on
large data sets this option can be very useful.

.. _one: ../../../python/stempy/test/fasta/yst09m.fasta







.. _number-of-motifs:

Multiple motifs
---------------

STEME can find more than one motif via the ``--num-motifs`` option. After each motif is found, STEME erases the
instances of the motif from its model. You can control the threshold used to find these instances via the
``--prediction-Z-threshold`` option.





.. index:: Controlling the running time
.. _controlling-the-running-time:

Controlling the running time
----------------------------

By default STEME examines every possible word as a starting point for the EM algorithm. As the word length and the size of
the data set grow this can take some time. STEME offers options to limit the length of time taken searching for these starting
seeds. Currently STEME searches for starts for each motif it is looking for, so if you set the option 
``--max-start-finding-time=3600`` it will look for one hour for each motif.





.. index:: Characteristics of the motifs
.. _motif-characteristics:

Characteristics of the motifs
-----------------------------

By default STEME will search for motifs between width 6 and 14. This can be changed via the ``--minw`` and ``--maxw`` options.
Additionally the ``--min-sites`` and ``--max-sites`` concentrate STEME's search to motifs that have a certain number of instances.
STEME chooses defaults for these settings based on the number of sequences but if STEME is finding motifs with too
few instances you can raise the lower limit. Raising the maximum number of sites can also be useful but this can slow down the
significance calculations. If not specified, STEME sets the ``--max-sites`` option to the number of sequences in the input set
with a cap of 10,000. The cap can be over-ridden on the command line.

If you know the consensus sequence of the motif you want to find, you can search just for this motif by using the ``--use-seed``
option. In this case, STEME will not perform a search for other seeds.





.. index:: EM algorithm

STEME as an implementation of the EM algorithm
----------------------------------------------

In the original STEME paper_ we described how STEME could be used as a replacement for the EM algorithm. We have some example
code_ showing this usage in ``python/scripts/steme-em``.  Try running this with the ``--help`` option to see a
list of possible options. If you want to run the pure EM algorithm without post-processing try the ``--dont-discretize``
option. By default, a discretization step is used that chooses the number of sites used in order to maximise the significance
of the motif.

For example, you could try::

  steme-em --help

or::

  steme-em $HOME/local/src/STEME/python/test/fasta/T00759-small.fa AACCTTGG 32

.. _code: ../../../python/scripts/steme-em




