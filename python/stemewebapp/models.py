# -*- coding: utf-8 -*-
#
# Copyright John Reid 2011, 2012
#

"""
Models for STEME web app.
"""


from uuid import uuid4 
from datetime import datetime

try:
    from flask.ext.sqlalchemy import SQLAlchemy
except NameError:
    from flaskext.sqlalchemy import SQLAlchemy

from .application import app
db = SQLAlchemy(app)



class Job(db.Model):
    """
    Database model of a job.
    """
    id = db.Column(db.Integer, primary_key=True)
    creation_date = db.Column(db.DateTime)
    name = db.Column(db.String(120))
    uuid = db.Column(db.Integer, unique=True)
    pid = db.Column(db.Integer)
    completed = db.Column(db.Boolean)

    def __init__(self, name=None):
        self.uuid = str(uuid4())
        self.pid = 0
        self.name = name
        self.completed = False
        self.creation_date = datetime.utcnow()

    def __repr__(self):
        return '<Job %r>' % (self.id)

