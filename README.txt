PhyPatClassProb is a C++ implementation of a dynamic programming algorithm for
calculating the probability of classes of patterns of a given length
parsimony and observed state combinations under a phylogenetic tree 
and model pair. 

Waddell, Ota, and Penny, 2009 develop a model adequacy assessment that uses
these probabilities, but they approximated them using simulation.  The algorithm
used here is to be described (Koch and Holder, in prep).


P. J. Waddell, R. Ota, and D. Penny, “Measuring Fit of Sequence Data to 
Phylogenetic Model: Gain of Power Using Marginal Tests,” Journal of Molecular
Evolution, vol. 69, no. 4, pp. 289–299, Oct. 2009.


Building
========
Building the software requires NCL, beagle, and pytbeaglehon. The easiest way
to build the software (on Mac or *nix) is to download the snapshot at
http://phylo.bio.ku.edu/software/pattern_class_prob_and_deps.tar.gz
which contains the dependencies and a do_build.sh script.

If you want to obtain the dependencies via version control, you can get NCL v2.1
by:

    svn co https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1

Beagle with:

    svn co http://beagle-lib.googlecode.com/svn/trunk/ beagle-lib-read-only
     
and pytbeaglehon using:
    
    git clone git://github.com/mtholder/pytbeaglehon.git
    

Running
=======
The -h flag passed to src/PhyPatClassProb or src/wrapPhyPatClassProb.py
will result in the list of command line arguments. The first two lines of
test/waddellEtAl/rag1wMLTreeOutput.txt show the invocation used to produce
that analysis of the dataset from Waddell et al's paper.




##############################################################################
## PhyPatClassProb
##
## Copyright 2011 Mark T. Holder.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
