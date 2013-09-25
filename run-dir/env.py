import os

for var, value in os.environ.iteritems():
    print '%20s : %s' % (var, value)
