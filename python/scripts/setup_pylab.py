#
# Copyright John Reid 2009, 2010
#

"""
Check STEME EM against MEME EM.
"""

import pylab
from math import sqrt

inches_per_pt = 1.0 / 72.27               # Convert pt to inch
golden_mean = (sqrt(5) - 1.0) / 2.0         # Aesthetic ratio


def setup_pylab_for_tex(fig_width_pt=246, normal_fontsize=10, small_fontsize=8):
    "Setup pylab to use tex backend with given fontsizes."
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})

    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean      # height in inches
    fig_size = [fig_width, fig_height]

    params = {
        'backend': 'ps',
        'axes.labelsize': normal_fontsize,
        'text.fontsize': normal_fontsize,
        'legend.fontsize': normal_fontsize,
        'xtick.labelsize': small_fontsize,
        'ytick.labelsize': small_fontsize,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    pylab.rcParams.update(params)


def correct_axes():
    "Correct axes to make sure we can see labels."
    pylab.axes([0.15, 0.2, 0.95 - 0.15, 0.95 - 0.2])
