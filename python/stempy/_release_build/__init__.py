#
# Copyright John Reid 2011, 2012
#

from cookbook.script_basics import boost_python_dlopen_flags
import os, logging
logger = logging.getLogger(__name__)


#
# Do the importing with the correct dlopen flags.
#
with boost_python_dlopen_flags():
#    import sys
#    logger.info('dlopen flags: %s', sys.getdlopenflags())
    USE_GENOME_ENV_VAR='STEME_USE_GENOME_INDEX'
    if USE_GENOME_ENV_VAR in os.environ:
        # Use genome index
        from ._index_genome_wrapper import *
    else:
        # Use standard index
        from ._index_standard_wrapper import *

    from ._stempy import *
    from ._stempy import _dummy_fn, _debug, _has_google_profiler
    if _has_google_profiler:
        from ._stempy import __google_profiler_start, __google_profiler_stop



