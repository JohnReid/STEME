#!/usr/bin/env python2

import logging
logging.basicConfig(level=logging.INFO)

from stemewebapp import app

#
# Set up error emailing and logging if not debug app.
#
if not app.debug:
    from logging.handlers import SysLogHandler

    syslog_handler = SysLogHandler()
    syslog_handler.setLevel(logging.WARNING)
    app.logger.addHandler(syslog_handler)

#    mail_handler = SMTPHandler(
#        '127.0.0.1',
#        'server-error@example.com',
#        ADMINS,
#        'YourApplication Failed'
#    )
#    mail_handler.setLevel(logging.ERROR)
#    app.logger.addHandler(mail_handler)

from logging.handlers import RotatingFileHandler
file_handler = RotatingFileHandler(
    'STEME_webapp.log', maxBytes=1000000, backupCount=10)
file_handler.setLevel(logging.INFO)
app.logger.addHandler(file_handler)

app.run()
