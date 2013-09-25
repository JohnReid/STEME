# -*- coding: utf-8 -*-
#
# Copyright John Reid 2011
#

"""
Forms for STEME web app.
"""


from flaskext.wtf.html5 import EmailField
from flaskext.wtf import Form, TextField, FileField, IntegerField, validators, ValidationError

import logging
logger = logging.getLogger(__name__)


class NewJobForm(Form):
    #email = EmailField(u'Please supply an email address<br>to refer back to this analysis later:', [validators.required()])
    name = TextField(u'Name to identify this analysis:<br>(optional)', [validators.Length(max=120)])
    sequences = FileField(u'Input sequences:<br>(FASTA format):', [validators.file_required(message=u'No file selected')])        
    num_motifs = IntegerField(u'Maximum number of motifs to find:', [validators.required()], default=1)
    bg_fasta_file = FileField(u'Background sequences:<br>(optional but recommended; FASTA format)')
    bg_model_order = IntegerField(u'Order of background Markov model:', [validators.required()], default=3)
    min_sites = IntegerField(u'Minimum number of sites:<br>(use 0 for default settings)', default=0)
    max_sites = IntegerField(u'Maximum number of sites:<br>(use 0 for default settings)', [validators.NumberRange(max=20000)], default=0)
    min_w = IntegerField(u'Minimum motif width:', default=6)
    max_w = IntegerField(u'Maximum motif width:', default=14)
    max_start_finding_time = IntegerField(
        u'Maximum start finding time:<br>(seconds; 0 for unlimited)',
        default=65000
    )
    

