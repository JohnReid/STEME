# -*- coding: utf-8 -*-
#
# Copyright John Reid 2011
#

"""
Default settings for STEME web app.
"""

# configuration
DEBUG = True
SECRET_KEY = 'Vv\xd1\xe8X\xc4d\x8d\x8c\x1f\x14\x1a\xed\x95\x8br\xbbx\xe6\x0b\xd9\x1f\\\xee'
USERNAME = 'admin'
PASSWORD = 'default'
ADMINS = ['john.reid@mrc-bsu.cam.ac.uk']
SQLALCHEMY_DATABASE_URI = 'sqlite:////tmp/STEME.db'
JOB_FOLDER = '/home/john/Dev/MyProjects/STEME/www/jobs'
STEME_SCRIPT = '/home/john/Dev/MyProjects/STEME/python/scripts/steme'
PYTHON_EXE = '/usr/bin/python2'
MOTIF_DBS = [
    '/home/john/Data/Jaspar/JASPAR_core.meme',
    '/home/john/Data/UNIPROBE/All_PWMs/UNIPROBE-all.meme',
]
