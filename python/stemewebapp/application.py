# -*- coding: utf-8 -*-
#
# Copyright John Reid 2011
#

"""
STEME web application.
"""

import os
from flask import Flask


#
# Create and configure our application
#
app = Flask(__name__.split('.')[0])
"The web application."

app.popens = []
"The list of Popen objects for subprocesses."

app.config.from_object('stemewebapp.default_settings')
if 'STEMEWEBAPP_SETTINGS' in os.environ:
    app.config.from_envvar('STEMEWEBAPP_SETTINGS')

