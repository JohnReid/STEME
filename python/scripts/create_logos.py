import os
import numpy
import glob


def is_file_newer(newer, older):
    if not os.path.exists(older):
        return True
    else:
        return os.path.getmtime(newer) > os.path.getmtime(older)


def logo_from_dist(dist, eps_target, png_target):
    import hmm.weblogo
    hmm.weblogo.weblogo_from_dist(dist, eps_target)
    hmm.weblogo.convert_format(eps_target, png_target)


def nucleo_dist_freqs(nucleo_dist):
    return map(nucleo_dist.get_freq, xrange(4))


def pssm_freqs(pssm):
    return numpy.array(map(nucleo_dist_freqs, pssm))


def logo_for_pssm(pssm, basename):
    logo_from_dist(pssm_freqs(pssm.pssm), basename + '.eps', basename + '.png')


def logo_for_pssm_name(pssm_name):
    import biopsy
    logo_for_pssm(biopsy.get_pssm(pssm_name), pssm_name)

for freq_file in glob.glob('*.freqs'):
    print freq_file
    basename = os.path.splitext(os.path.basename(freq_file))[0]
    eps_target = 'eps/%s.eps' % basename
    png_target = 'png/%s.png' % basename
    if is_file_newer(freq_file, eps_target):
        print eps_target
        dist = numpy.array([map(float, l.split()) for l in open(freq_file)])
        logo_from_dist(dist, eps_target, png_target)
