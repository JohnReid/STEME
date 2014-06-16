#!/usr/bin/env python

"""
Code for importance sampling in EM part of MEME algorithm.
"""

import logging
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

from itertools import imap
import pandas as pd
import pandas.rpy.common as com
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seqan.traverse
import numpy as npy
import numpy.random as rdm
from numpy import log, exp
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
from rpy2.robjects.packages import importr
grdevices = importr('grDevices')


ALLBASES = (seqan.DNA('A'), seqan.DNA('C'), seqan.DNA('G'), seqan.DNA('T'))
SIGMA = len(ALLBASES)
LOGQUARTER = log(.25)
UNIFORM0ORDER = npy.ones(SIGMA) / SIGMA


def logo(dist, tag, make_png=False, make_eps=True, write_title=True):
    "Generate a logo with the given tag in the given directory."
    import weblogolib as W
    import corebio.seq
    data = W.LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist)
    scale = 5.4 * 4
    options = W.LogoOptions(
        logo_title=write_title and tag or None,
        stack_width=scale,
        stack_aspect_ratio=5,
        color_scheme=W.colorscheme.nucleotide,
        show_xaxis=False,
        show_yaxis=False,
        show_fineprint=False,
    )
    format_ = W.LogoFormat(data, options)
    filename = 'logo-%s' % tag
    if make_eps:
        W.eps_formatter(data, format_, open('%s.eps' % filename, 'w'))
    if make_png:
        W.png_formatter(data, format_, open('%s.png' % filename, 'w'))


def uniform0orderloglikelihood(X):
    """A likelihood that assigns 1/4 probability to each base at each position."""
    return len(X) * LOGQUARTER


def loglikelihoodforpwm(pwm):
    """Return a function that computes the log likelihood for the pwm."""
    logpwm = log(pwm)
    def loglikelihood(X):
        return sum(logpwm[w,base.ordValue] for w, base in enumerate(X))
    return loglikelihood


class ZnSumVisitor(object):
    def __init__(self, W, Zncalculator):
        self.W = W
        self.Zncalculator = Zncalculator
        self.sums = npy.zeros((W, 4))

    def __call__(self, it):
        if it.repLength >= self.W:
            # Get the word
            Xn = it.representative[:W]
            # Calculate Zn
            Zn = self.Zncalculator(Xn)
            # Update sums
            for w, base in enumerate(Xn):
                self.sums[w,base.ordValue] += Zn * it.numOccurrences
            # Have gone deep enough in index, truncate traversal
            return False
        else:
            # Keep descending
            return True


def createZncalculator(pwm, lambda_):
    bsloglikelihoodfn = loglikelihoodforpwm(pwm)
    lambdaratio = lambda_ / (1. - lambda_)
    def calculateZn(X):
        logodds = bsloglikelihoodfn(X) - uniform0orderloglikelihood(X)
        return 1./(1 + 1/(lambdaratio * exp(logodds)))
    return calculateZn


def addpseudocounts(pwm, numsites, pseudocount):
    return (pwm * numsites + pseudocount) / (numsites + pseudocount)


def normalisepwm(pwm):
    return (pwm.T / pwm.sum(axis=1)).T


lambda_ = .01
numsites = 50
pseudocount = 1.
runx1pwm = npy.array((
    (0.384615,  0.076923,  0.115385,  0.423077),
    (0.461538,  0.076923,  0.038462,  0.423077),
    (0.153846,  0.269231,  0.038462,  0.538462),
    (0.038462,  0.038462,  0.000000,  0.923077),
    (0.076923,  0.000000,  0.884615,  0.038462),
    (0.076923,  0.307692,  0.000000,  0.615385),
    (0.000000,  0.000000,  1.000000,  0.000000),
    (0.000000,  0.000000,  1.000000,  0.000000),
    (0.000000,  0.038462,  0.000000,  0.961538),
    (0.307692,  0.076923,  0.000000,  0.615385),
    (0.500000,  0.076923,  0.153846,  0.269231),
))
runx1withpc = addpseudocounts(runx1pwm, numsites, pseudocount)
W = len(runx1withpc)

#logo(runx1pwm, 'runx1')
#logo(runx1withpc, 'runx1-pc')

#seqs = seqan.StringDNASet(('AAAAAAAA', 'ACGTACGT', 'TATATATA'))
numbases, seqs, ids = seqan.readFastaDNA('T00759-small.fa')
index = seqan.IndexStringDNASetESA(seqs)

calculateZn = createZncalculator(runx1withpc, lambda_)
sumvisitor = ZnSumVisitor(W, calculateZn)
seqan.traverse.topdownhistorytraversal(index.topdownhistory(), sumvisitor)
print sumvisitor.sums
trueZnsum = sumvisitor.sums[0].sum()
print trueZnsum
#logo(normalisepwm(sumvisitor.sums), 'learnt')


def countWmers(it, W, counts):
    if it.repLength >= W:
        count = it.numOccurrences
    else:
        count = 0
        if it.goDown():
            while True:
                count += countWmers(it, W, counts)
                if not it.goRight():
                    break
            it.goUp()
    counts[it.value.id] = count
    return count


Wmercounts = npy.zeros(2*len(index), dtype=npy.uint32)
numWmers = countWmers(index.topdownhistory(), W, Wmercounts)


def createpwmlikelihoodfn(pwm):
    def baselikelihoodfn(w):
        return pwm[w]
    return baselikelihoodfn


def bglikelihoodfn(w):
    return UNIFORM0ORDER


class ImportanceSampler(object):
    def __init__(self, baselikelihoodfn, W, Wmercounts):
        self.W = W
        self.Wmercounts = Wmercounts
        self.baselikelihoodfn = baselikelihoodfn

    def __call__(self, it, likelihoodratio=1.):
        w = it.repLength
        if w >= self.W:
            return it, likelihoodratio
        else:
            representative = it.representative
            # Descend each child to get Wmer counts
            counts = npy.zeros(SIGMA)
            for base in ALLBASES:
                if it.goDown(base):
                    counts[base.ordValue] = Wmercounts[it.value.id]
                    #logger.info('%s %d %f', base, Wmercount, likelihoodratio)
                    it.goUp()
            freqs = counts.astype(float) / counts.sum()
            # Get the importance weights and combine with Wmer frequencies
            # to create sampling distribution
            impweights = self.baselikelihoodfn(w)
            weights = freqs * impweights
            samplingdist = weights / weights.sum()
            #logger.info('%-10s: %s', representative, samplingdist)
            # Sample one of the bases
            sample = rdm.choice(ALLBASES, p=samplingdist)
            # Descend the sample
            wentDown = it.goDown(sample)
            assert wentDown
            # Calculate the updated likelihood ratio
            likelihoodratioupdate = samplingdist[sample.ordValue] / freqs[sample.ordValue]
            # recurse
            return self(it, likelihoodratio * likelihoodratioupdate)



def sample(sampler, numsamples, **kwargs):
    sampled = [
        sampler(index.topdownhistory())
        for _ in xrange(numsamples)]
    sampledit = map(lambda x: x[0], sampled)
    lr = map(lambda x: x[1], sampled)
    X = map(lambda it: it.representative[:W], sampledit)
    Z = map(calculateZn, X)
    kwargs.update({
        'X': X,
        'Z': Z,
        'lr': lr,
    })
    df = pd.DataFrame(kwargs)
    return df


numsamples = 1000
rdm.seed(1)
samplerbs = ImportanceSampler(createpwmlikelihoodfn(runx1withpc), W, Wmercounts)
samplerbg = ImportanceSampler(bglikelihoodfn, W, Wmercounts)
samplebs = sample(samplerbs, numsamples, model='BS')
samplebg = sample(samplerbg, numsamples, model='BG')
samples = pd.concat((samplebs, samplebg))
samples.index = npy.arange(len(samples))  # re-index to avoid duplicate row.names in Rdf
rsamples = com.convert_to_r_dataframe(samples)

#grdevices.png(file="sampled-Z.png", width=4, height=3, units="in", res=300)
pp = ggplot2.ggplot(rsamples) + \
    ggplot2.aes_string(x='Z', color='factor(model)') + \
    ggplot2.scale_colour_discrete(name="model") + \
    ggplot2.geom_density() + \
    ggplot2.scale_x_log10()
    #ggplot2.scale_x_continuous(limits=FloatVector((0, 1)))
pp.plot()
#grdevices.dev_off()


def estimatesum(samples):
    return sum(samples['Z']/samples['lr'])*numWmers/len(samples)


def makeestimate(sampler, numsamples, **kwargs):
    samples = sample(sampler, numsamples, **kwargs)
    return estimatesum(samples)


def makeestimates(sampler, numsamples, numestimates, **kwargs):
    estimates = [makeestimate(sampler, numsamples, **kwargs) for _ in xrange(numestimates)]
    kwargs.update({
        'estimate': estimates,
        'numsamples': numsamples,
    })
    return pd.DataFrame(kwargs)


def estimatebothsamplers(numsamples, numestimates):
    estimatesbs = makeestimates(
        samplerbs, numsamples=numsamples, numestimates=numestimates, model='BS')
    estimatesbg = makeestimates(
        samplerbg, numsamples=numsamples, numestimates=numestimates, model='BG')
    return pd.concat((estimatesbs, estimatesbg))

estimates = pd.concat(
    (estimatebothsamplers(numsamples, 10) for numsamples in [10, 50, 100, 150, 300, 1000]))
bp = estimates.boxplot(column=['estimate'], by=['numsamples', 'model'])
ax = plt.gca()
xmin, xmax = ax.get_xbound()
plt.hlines(trueZnsum, xmin, xmax, linestyles='dashdot')
ax.set_ylim(top=2*trueZnsum)
plt.savefig('estimates.png')

grouped = estimates.groupby(['numsamples', 'model'])
grouped.aggregate(npy.var)
grouped.aggregate(npy.mean)

varbs = npy.var(estimatesbs['estimate'])
varbg = npy.var(estimatesbg['estimate'])
print varbs
print varbg
print varbg / varbs

mysamples = sample(samplerbs, 300*20*2, model='BS')
myestimate = estimatesum(samples)
print myestimate
print trueZnsum


