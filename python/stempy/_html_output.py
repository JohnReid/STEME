# -*- coding: utf-8 -*-
#
# Copyright John Reid 2011, 2012
#

"""
HTML output.
"""

import os
import shutil
import logging
from . import ensure_dir_exists, OptionParser, add_options, get_default_options, __release__
from jinja2 import Environment, PackageLoader
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)
env = Environment(loader=PackageLoader('stempy', 'templates'))


class HTMLOutput(object):

    """
    Produces HTML output.
    """

    def __init__(self, output_dir):
        self.output_dir = output_dir
        "Where to write the files."

    def _static_dir(self):
        """
        The directory where static web files are kept.
        """
        return os.path.join(self.output_dir, 'static')

    def initialise(self, algorithm):
        """
        Create directories and copy static files.
        """
        from pkg_resources import resource_filename
        stylesheet_filename = resource_filename(__name__, 'static/style.css')
        ensure_dir_exists(self._static_dir())
        shutil.copy(stylesheet_filename, self._static_dir())

    @staticmethod
    def _image_url(idx, format_='png'):
        """
        A URL for the motif's image.
        """
        return u'logo-STEME-motif-%02d.%s' % (idx, format_)

    def found_motif(self, algorithm, motif, seconds_taken):
        """
        Write the motif
        """
        template = env.get_template('motif.html')
        f = open(
            os.path.join(self.output_dir, 'STEME-motif-%02d.html' % motif.idx), 'w')
        f.write(
            template.render(
                motif=motif,
                image_url=HTMLOutput._image_url(motif.idx)
            )
        )

    def _write_options(self, options):
        """
        Write the HTML options document.
        """
        template = env.get_template('options.html')
        filename = os.path.join(self.output_dir, 'options.html')
        logger.info('Writing STEME options as HTML to %s', filename)
        f = open(filename, 'w')
        option_parser = OptionParser()
        add_options(option_parser)
        variables = {
            'options': options,
            'option_parser': option_parser,
            'getattr': getattr,
            'str': str,
            'default_options': get_default_options(),
        }
        f.write(template.render(**variables))

    def _write_job(self, algorithm):
        """
        Write the HTML.
        """
        template = env.get_template('job.html')
        filename = os.path.join(self.output_dir, 'STEME.html')
        logger.info('Writing output as HTML to %s', filename)
        f = open(filename, 'w')
        tomtom_url = algorithm.options.tomtom and 'TOMTOM/tomtom.html' or None
        num_bases, num_seqs, freqs, unknown = algorithm.input_sequences.composition(
        )
        bg_num_bases, bg_num_seqs, bg_freqs, bg_unknown = algorithm.background_composition(
        )
        variables = {
            'logfile': u'STEME.log',
            'meme_like_output': algorithm.options.meme_like_output,
            'tomtom_url': tomtom_url,
            'num_bases': num_bases,
            'num_sequences': num_seqs,
            'freqs': freqs,
            'unknown': unknown,
            'bg_num_bases': bg_num_bases,
            'bg_num_sequences': bg_num_seqs,
            'bg_freqs': bg_freqs,
            'bg_unknown': bg_unknown,
            'calc_rel_enrichment': calc_rel_enrichment,
            'timings': algorithm.timings,
            'total_duration': algorithm.total_duration,
            'version': __release__,
            'completion_date_time': datetime.now().strftime('%c'),
            'total_duration_string': time_delta_readable(algorithm.total_duration)
        }
        f.write(
            template.render(
                motifs=[(motif, HTMLOutput._image_url(idx), HTMLOutput._image_url(idx, 'eps'))
                        for idx, motif in enumerate(algorithm.motifs)],
                **variables
            )
        )

    def finalise(self, algorithm):
        """
        Write the HTML.
        """
        self._write_job(algorithm)
        self._write_options(algorithm.options)


def calc_rel_enrichment(num_bases, bg_num_bases, motif):
    return motif.input_stats.per_base / motif.bg_stats.per_base


def time_delta_readable(total):
    """Return a readable string representing the total number of seconds.
    """
    total = int(total)
    days = total / 60 / 60 / 24
    hours = total / 60 / 60 % 24
    minutes = total / 60 % 60
    seconds = total % 60
    result = ''
    if days:
        result = '%d days, ' % days
    if days or hours:
        result += '%d hours ' % hours
    if days or hours or minutes:
        result += '%d minutes ' % minutes
    return result + '%d seconds' % seconds
