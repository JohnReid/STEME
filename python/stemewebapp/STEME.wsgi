#
# Copyright John Reid 2012
#

#
# Let STEME web application know where configuration is
#
import os
os.environ['STEMEWEBAPP_SETTINGS'] = '/etc/STEME/settings.py'

#
# Import the application
#
from stemewebapp import app as application
