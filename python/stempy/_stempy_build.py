#
# Copyright John Reid 2009, 2010, 2011, 2012
#

"""
Code to choose whether to import debug or release build of stempy module.
"""

import sys
import logging
logger = logging.getLogger(__name__)

#
# decide whether to import debug or release stempy C++ extension
#
# only available in python debug build
_python_debug_build = hasattr(sys, "gettotalrefcount")
if _python_debug_build:

    from ._debug_build import *
    from ._debug_build import _dummy_fn, _debug, _has_google_profiler
    from . import _debug_build as S
    if _has_google_profiler:
        from _debug_build import __google_profiler_start, __google_profiler_stop

else:

    from ._release_build import *
    from ._release_build import _dummy_fn, _debug, _has_google_profiler
    from . import _release_build as S
    if _has_google_profiler:
        from _release_build import __google_profiler_start, __google_profiler_stop

logger.info('Loaded STEME C++-python interface from %s', S.__name__)
