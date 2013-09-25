#
# Copyright John Reid 2012, 2013
#


"""
Code to parse minimal MEME motif format (http://meme.sdsc.edu/meme/doc/meme-format.html):

    #
    # Parse and extract the data from the sample file
    #
    from pkg_resources import resource_string
    extracted = do_parse_and_extract(resource_string('stempy', 'test/sample-dna-motif.meme-io'))
    print extracted
    

Example of format:

    MEME version 4
    
    ALPHABET= ACGT
    
    strands: + -
    
    Background letter frequencies
    A 0.303 C 0.183 G 0.209 T 0.306 
    
    MOTIF crp
    letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009 
     0.000000  0.176471  0.000000  0.823529 
     0.000000  0.058824  0.647059  0.294118 
     0.000000  0.058824  0.000000  0.941176 
     0.176471  0.000000  0.764706  0.058824 
     0.823529  0.058824  0.000000  0.117647 
     0.294118  0.176471  0.176471  0.352941 
     0.294118  0.352941  0.235294  0.117647 
     0.117647  0.235294  0.352941  0.294118 
     0.529412  0.000000  0.176471  0.294118 
     0.058824  0.235294  0.588235  0.117647 
     0.176471  0.235294  0.294118  0.294118 
     0.000000  0.058824  0.117647  0.823529 
     0.058824  0.882353  0.000000  0.058824 
     0.764706  0.000000  0.176471  0.058824 
     0.058824  0.882353  0.000000  0.058824 
     0.823529  0.058824  0.058824  0.058824 
     0.176471  0.411765  0.058824  0.352941 
     0.411765  0.000000  0.000000  0.588235 
     0.352941  0.058824  0.000000  0.588235 
    
    MOTIF lexA
    letter-probability matrix: alength= 4 w= 18 nsites= 14 E= 3.2e-035 
     0.214286  0.000000  0.000000  0.785714 
     0.857143  0.000000  0.071429  0.071429 
     0.000000  1.000000  0.000000  0.000000 
     0.000000  0.000000  0.000000  1.000000 
     0.000000  0.000000  1.000000  0.000000 
     0.000000  0.000000  0.000000  1.000000 
     0.857143  0.000000  0.071429  0.071429 
     0.000000  0.071429  0.000000  0.928571 
     0.857143  0.000000  0.071429  0.071429 
     0.142857  0.000000  0.000000  0.857143 
     0.571429  0.071429  0.214286  0.142857 
     0.285714  0.285714  0.000000  0.428571 
     1.000000  0.000000  0.000000  0.000000 
     0.285714  0.214286  0.000000  0.500000 
     0.428571  0.500000  0.000000  0.071429 
     0.000000  1.000000  0.000000  0.000000 
     1.000000  0.000000  0.000000  0.000000 
     0.000000  0.000000  0.785714  0.214286 
"""


from pyparsing import Suppress, CaselessLiteral, ZeroOrMore, OneOrMore, Group, \
    Literal, Keyword, \
    Word, Optional, Combine, Empty, White, \
    alphas, printables, nums, restOfLine
from numpy import asarray
from cookbook import namedtuple


def K(literal):
    return Suppress(Keyword(literal))

comment = ZeroOrMore(Suppress(Literal('#') + restOfLine))
point = Literal(".")
e = CaselessLiteral("E")
fnumber = Combine(
    Word("+-" + nums, nums) +
    Optional(point + Optional(Word(nums))) +
    Optional(e + Word("+-" + nums, nums))
) | Literal('inf')
version = K('MEME version') + Word(printables)('VERSION')
alphabet = Group(K('ALPHABET=') + Word(alphas))('ALPHABET')
strands = Group(K('strands:')
                + Optional(Keyword('+')) + Optional(Keyword('-')))('STRANDS')
background_freqs = Group(K('Background letter frequencies')
                         + Suppress(restOfLine)
                         + K('A') + fnumber + K('C') + fnumber + K('G')
                         + fnumber + K('T')
                         + fnumber)('BACK_FREQS')
matrix_row = Group(
    fnumber + fnumber + fnumber + fnumber
)
prob_matrix = Group(
    Optional(K('alength=') + Word(nums)('ALENGTH'))
    + Optional(K('w=') + Word(nums)('W'))
    + Optional(K('nsites=') + Word(nums)('NSITES'))
    + Optional(K('E=') + fnumber('E'))
    + Group(OneOrMore(matrix_row))('ROWS')
)
letter_probs = K('letter-probability matrix:') + prob_matrix('LETTER_PROBS')
log_odds = K('log-odds matrix:') + prob_matrix('LOG_ODDS')
prob_matrix = letter_probs | log_odds
url = Optional(K('URL') + Word(printables)('URL'))
motif = Group(
    K('MOTIF') + Word(printables)('NAME') + Optional(restOfLine('ALTNAME'))
    + comment
    + Group(OneOrMore(prob_matrix))('MATRICES')
    + comment
    + url
)
meme_format = \
    comment \
    + version \
    + comment \
    + Optional(alphabet) \
    + comment \
    + Optional(strands) \
    + comment \
    + Optional(background_freqs) \
    + comment \
    + Group(OneOrMore(motif))('MOTIFS') \


Parsed = namedtuple('Parsed', 'version alphabet strands back_freqs motifs')
Motif = namedtuple('Motif', 'name altname letter_probs log_odds url')
Matrix = namedtuple('Matrix', 'alength w nsites E values')


def extract_or_none(thing, name, extractor=lambda x: x):
    """Try to extract named value out of thing but return None if it is not
    there.
    """
    if name in thing:
        return extractor(thing[name])
    return None


def extract_motif(motif):
    return Motif(
        motif['NAME'],
        altname=extract_or_none(motif, 'ALTNAME', str.strip),
        letter_probs=extract_or_none(
            motif['MATRICES'],
            'LETTER_PROBS',
            extract_matrix),
        log_odds=extract_or_none(
            motif['MATRICES'], 'LOG_ODDS', extract_matrix),
        url=extract_or_none(motif, 'URL')
    )


def extract_matrix(matrix):
    return Matrix(
        alength=extract_or_none(matrix, 'ALENGTH', int),
        w=extract_or_none(matrix, 'W', int),
        nsites=extract_or_none(matrix, 'NSITES', int),
        E=extract_or_none(matrix, 'E', float),
        values=asarray([map(float, row) for row in matrix['ROWS']]),
    )


def extract_parsed(parsed):
    return Parsed(
        parsed['VERSION'],
        alphabet=extract_or_none(parsed, 'ALPHABET', lambda x: x[0]),
        strands=extract_or_none(parsed, 'STRANDS', lambda x: set(x)),
        back_freqs=extract_or_none(
            parsed, 'BACK_FREQS', lambda x: map(float, x)),
        motifs=extract_or_none(
            parsed, 'MOTIFS', lambda x: map(extract_motif, x))
    )


def do_parse_and_extract(string):
    """Parse and extract info from string."""
    parsed = meme_format.parseString(string)
    return extract_parsed(parsed)


def write_pwm(out, motif):
    """Write one PWM out. For example:
    
        MOTIF crp
        letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009 
         0.000000  0.176471  0.000000  0.823529 
         0.000000  0.058824  0.647059  0.294118 
         0.000000  0.058824  0.000000  0.941176 
         0.176471  0.000000  0.764706  0.058824 
         0.823529  0.058824  0.000000  0.117647 
         0.294118  0.176471  0.176471  0.352941 
         0.294118  0.352941  0.235294  0.117647 
         0.117647  0.235294  0.352941  0.294118 
         0.529412  0.000000  0.176471  0.294118 
         0.058824  0.235294  0.588235  0.117647 
         0.176471  0.235294  0.294118  0.294118 
         0.000000  0.058824  0.117647  0.823529 
         0.058824  0.882353  0.000000  0.058824 
         0.764706  0.000000  0.176471  0.058824 
         0.058824  0.882353  0.000000  0.058824 
         0.823529  0.058824  0.058824  0.058824 
         0.176471  0.411765  0.058824  0.352941 
         0.411765  0.000000  0.000000  0.588235 
         0.352941  0.058824  0.000000  0.588235 
    """
    print >>out, \
        """
MOTIF %s %s
letter-probability matrix: alength= %d w= %d nsites= %d %s
%s
""" % (
        motif.name, motif.altname,
        motif.letter_probs.alength,
        motif.letter_probs.w,
        motif.letter_probs.nsites,
        motif.letter_probs.E and 'E= %e' % motif.letter_probs.E or '',
        '\n'.join(
            ' '.join(('%.6f' % p) for p in row)
            for row in motif.letter_probs.values
        )
    )


def write_pwms(out, meme_info):
    """Write the motifs out in the same format.
    """
    print >>out, \
        """
MEME version %s

ALPHABET= %s

strands: %s

Background letter frequencies
A %.3f C %.3f G %.3f T %.3f

""" % (
        meme_info.version,
        meme_info.alphabet,
        ' '.join(meme_info.strands),
        meme_info.back_freqs[0], meme_info.back_freqs[
            1], meme_info.back_freqs[2], meme_info.back_freqs[3],
    )
    for motif in meme_info.motifs:
        write_pwm(out, motif)


if '__main__' == __name__:
    #
    # Parse the sample file
    #
    from pkg_resources import resource_string
    extracted = do_parse_and_extract(
        resource_string('stempy', 'test/sample-dna-motif.meme-io'))
    print extracted
