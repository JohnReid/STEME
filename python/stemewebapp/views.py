# -*- coding: utf-8 -*-
#
# Copyright John Reid 2011
#

"""
Views for STEME web app.
"""

import os, logging, subprocess, psutil
from flask import render_template, flash, url_for, redirect, Markup, request
from werkzeug import secure_filename
from sqlalchemy.sql import desc

from .application import app
from .models import Job, db
from .forms import NewJobForm

logger = logging.getLogger(__name__)

def upload_file(f, d):
    """
    Upload a file to given directory (or default directory if none given).
    """
    filename = secure_filename(f.filename)
    path = os.path.join(d, filename)
    logger.info('Saving uploaded file to %s', path)
    f.save(path)
    return filename, path


def job_directory(uuid):
    """
    Return the directory for the job.
    """
    d = os.path.join(app.config['JOB_FOLDER'], str(uuid))
    return d


def job_file(uuid, *args):
    """
    Return a filename in the job directory.
    """
    return os.path.join(job_directory(uuid), *args)


def job_url(uuid):
    """
    The URL for a job.
    """
    return 'result/%s/STEME.html' % uuid


def job_log_url(uuid):
    """
    The URL for a job.
    """
    return 'result/%s/steme.log' % uuid



def named_job(job, link=True):
    """
    Return some markup that refers to a job.
    """
    html = job_file(job.uuid, 'STEME.html')
    if link and os.path.exists(html) and os.stat(html).st_size > 0:
        url = job_url(job.uuid)
        if job.name:
            return Markup('Job <a href="%s">%s</a>' % (url, job.name))
        else:
            return Markup('Unnamed <a href="%s">job</a>' % url)
    else:
        if job.name:
            return Markup('Job "%s"' % job.name)
        else:
            return Markup('Unnamed job')


def check_jobs_finished():
    """
    Check which jobs have finished and updates database.
    """
    #
    # check each popen we have in our list
    #
    for popen in app.popens:
        retval = popen.poll()
        if None != retval:
            logger.info('Job with pid %d has finished', popen.pid)
            for job in Job.query.filter_by(pid=popen.pid).all():
                job.completed = True # set completed flag on job to true
            app.popens.remove(popen) # remove from list when finished
            flash('Job "%s" has finished' % job.name)
    db.session.commit()
    
    #
    # Check each job that still has not completed
    #
    for job in Job.query.filter_by(completed=False).all():
        try:
            psutil.Process(job.pid)
        except psutil.NoSuchProcess:
            logger.warning('Could not find process for job %s with pid %d', job.uuid, job.pid)
            job.completed = True
            flash('Job "%s" has finished' % job.name)
    db.session.commit()




def check_job_successful(job):
    """
    Check that the job finished successfully.
    """
    job_dir = job_file(job.uuid)
    html = job_file(job.uuid, 'STEME.html')
    meme = job_file(job.uuid, 'meme.txt')
    successful = (
        os.path.exists(job_dir)
        and os.path.exists(html)
        and os.stat(html).st_size # non-empty HTML
        and os.path.exists(meme)
    )
    return successful



def run_steme(form, job):
    """
    Run the STEME algorithm.
    """
    # set up arguments
    seqs_filename, _path = upload_file(form.sequences.file, job_directory(job.uuid))
    args = [
        app.config['PYTHON_EXE'],
        app.config['STEME_SCRIPT'],
        '--output-dir', '.',
        '--bg-model-order', form.bg_model_order.data,
        '--num-motifs', form.num_motifs.data,
        '--minw', form.min_w.data,
        '--maxw', form.max_w.data,
        '--max-start-finding-time', form.max_start_finding_time.data,
    ]
    
    # add min/max # of sites if given
    if form.min_sites.data:
        args.append('--min-sites')
        args.append(form.min_sites.data)
    if form.max_sites.data:
        args.append('--max-sites')
        args.append(form.max_sites.data)
    
    # Add arguments for background FASTA file if given
    if form.bg_fasta_file.file:
        bg_filename, _path = upload_file(form.bg_fasta_file.file, job_directory(job.uuid))
        args.append('--bg-fasta-file')
        args.append(bg_filename)
    
    # Add TOMTOM arguments
    for motif_db in app.config['MOTIF_DBS']:
        args.append('--tomtom')
        args.append(motif_db)
    
    # Add sequence file
    args.append(seqs_filename)
        
    args = map(str, args) # convert to strings
    logger.info('Starting STEME job %s with args: %s', job.uuid, ' '.join(('"%s"' % a) for a in args))
    stdout = open(job_file(job.uuid, 'STEME.out'), 'w')
    if 'PYTHONPATH' in app.config:
        os.environ['PYTHONPATH'] = app.config['PYTHONPATH']
    popen = subprocess.Popen(
        args, 
        bufsize=-1, 
        stdout=stdout, 
        stderr=subprocess.STDOUT, 
        close_fds=True,
        cwd=job_directory(job.uuid)
    )
    app.popens.append(popen) # remember popen to determine when completed
    return popen


@app.route('/')
def home():
    """
    The home view.
    """
    return render_template("home.html")


@app.route('/new', methods=['GET', 'POST'])
def new_job():
    """
    The new job view.
    """
    form = NewJobForm()
    if form.validate_on_submit():
        job = Job(form.name.data)
        job_dir = job_directory(job.uuid)
        if not os.path.exists(job_dir):
            os.mkdir(job_dir)
        logger.info('Adding new job uuid=%s submitted from %s to database', job.uuid, request.remote_addr)
        popen = run_steme(form, job)
        job.pid = popen.pid
        db.session.add(job)
        db.session.commit()
        flash(u'Started new job')
        return redirect(url_for('list_jobs'))
    return render_template("newjob.html", form=form)


@app.route('/jobs')
def list_jobs():
    """
    The list jobs view.
    """
    check_jobs_finished()
    return render_template(
        "listjobs.html", 
        jobs=Job.query.order_by(desc(Job.creation_date)).all(), 
        named_job=named_job,
        check_job_successful=check_job_successful,
        job_log_url=job_log_url,
    )



@app.route('/result/<jobid>/')
def job_result(jobid):
    """
    The view of the results of a job.
    """
    return open(os.path.join(app.config['JOB_FOLDER'], jobid, 'STEME.html')).read()



