#!/usr/bin/env python
# -*- coding: latin-1 -*-
#
# Copyright John Reid 2010, 2011, 2012, 2013
#

"""
setuptools setup script for STEME.
Adapted from http://git.tiker.net/pyublas.git/tree.
"""

import os, sys

# is it a python debug build?
_python_debug_build = hasattr(sys, "gettotalrefcount")
_build = _python_debug_build and 'debug' or 'release'


def read(*fnames):
    """
    Utility function to read the README file.
    Used for the long_description.  It's nice, because now 1) we have a top level
    README file and 2) it's easier to type in the README file than to put a raw
    string in below ...
    """
    return open(os.path.join(os.path.dirname(__file__), *fnames)).read()


def get_config_schema():
    from aksetup_helper import ConfigSchema, Option, \
            IncludeDir, LibraryDir, Libraries, BoostLibraries, \
            Switch, StringListOption, make_boost_base_options

    import sys
    if 'linux' in sys.platform:
        default_libs = ['rt']

    return ConfigSchema(
        make_boost_base_options() + [
            BoostLibraries("python"),
            Option("SEQAN_DIR", default=''),
            StringListOption(
                "CXXFLAGS", [], help="Any extra C++ compiler options to include"),
            StringListOption(
                "LDFLAGS", [], help="Any extra linker options to include"),
        ]
    )


def main():
    from aksetup_helper import hack_distutils, get_config, setup, NumpyExtension
    from setuptools import find_packages
    import sys

    hack_distutils()
    conf = get_config(get_config_schema())

    INCLUDE_DIRS = [
        'c++',
        os.path.join('c++', 'include'),
        os.path.join('c++', 'hmm', 'myrrh'),
        os.path.join('c++', 'hmm'),
        os.path.join('c++', 'indexing_suite_v2'),
        os.path.join('c++', 'FAST'),
        os.path.join(conf['SEQAN_DIR'], 'include'),
    ] + conf['BOOST_INC_DIR']
    LIBRARY_DIRS = conf['BOOST_LIB_DIR']
    LIBRARIES = conf['BOOST_PYTHON_LIBNAME'] + ['fftw3']
    CXXFLAGS = conf['CXXFLAGS']
    LDFLAGS = conf['LDFLAGS']
    EXTRA_DEFINES = [
        ("SEQAN_ENABLE_TESTING"  , "0"),
        ("SEQAN_ENABLE_DEBUG"    , "0"),
        ("BOOST_DISABLE_ASSERTS" , None),
        ("MYRRH_DISABLE_ASSERTS" , None),
        ("PVAL_LOOKUP"           , None),
        ("NDEBUG"                , None),
    ]

    # need to link against librt otherwise get an undefined aio_cancel
    # error or somesuch.
    if 'linux' in sys.platform:
        LIBRARIES.append('rt')

    # want to add these options to all GCC builds
    if 'linux' in sys.platform or 'darwin' in sys.platform:
        CXXFLAGS += [
            "-Wno-unused",
            "-Wno-deprecated",
            "-Wno-long-long",
            "-Wno-variadic-macros",
            "-Wno-sequence-point",
            "-Wno-uninitialized",
            "-finline-functions",
            "-ftemplate-depth-128",
        ]

    try:
        from distutils.command.build_py import build_py_2to3 as build_py
    except ImportError:
        # 2.x
        from distutils.command.build_py import build_py

    #
    # C++ extensions
    #
    lib_srcs = [
        'c++/hmm/myrrh/src/python/multi_array_to_numpy.cpp',
        'c++/FAST/HS.cpp',
        'c++/FAST/motif_evaluator.cpp',
        'c++/FAST/theta_fns.cpp',
        'c++/FAST/minimize.cpp',
        'c++/FAST/convolution.cpp',
        'c++/FAST/fft.cpp',
        'c++/FAST/utils.cpp',
        'c++/FAST/llr_score.cpp',
        'c++/pvalues/parse.cpp',
        'c++/pvalues/bejerano.cpp',
        'c++/pvalues/fast.cpp',
        'c++/pvalues/shifted_hirji.cpp',
        'c++/pvalues/pval.cpp',
        'c++/pvalues/pvalue_test_defs.cpp',
        'c++/python/steme.cpp',
    ]
    module_stempy_srcs = [
        'c++/python/module_stempy.cpp',
        'c++/python/python_bs.cpp',
        'c++/python/python_llr_pvalues.cpp',
        'c++/python/python_utility.cpp',
    ]
    module_index_srcs = [
        'c++/python/module_index.cpp',
        'c++/python/python_descender.cpp',
        'c++/python/python_bg.cpp',
        'c++/python/python_data.cpp',
        'c++/python/python_em.cpp',
        'c++/python/python_find_best_w_mers.cpp',
        'c++/python/python_find_instances.cpp',
        'c++/python/python_markov.cpp',
        'c++/python/python_model.cpp',
        'c++/python/python_seqan.cpp',
        'c++/python/python_significance.cpp',
        'c++/python/python_start_finder.cpp',
    ]
    cStempy = NumpyExtension(
        'stempy._%s_build._stempy' % _build,
        lib_srcs + module_stempy_srcs,
        include_dirs         = INCLUDE_DIRS,
        library_dirs         = LIBRARY_DIRS,
        libraries            = LIBRARIES,
        define_macros        = EXTRA_DEFINES,
        extra_compile_args   = CXXFLAGS,
        extra_link_args      = LDFLAGS,
    )
    cIndex = NumpyExtension(
        'stempy._%s_build._index' % _build,
        module_index_srcs,
        include_dirs         = INCLUDE_DIRS,
        library_dirs         = LIBRARY_DIRS,
        libraries            = LIBRARIES,
        define_macros        = EXTRA_DEFINES + [('STEME_INDEX_MODULE_NAME', '_index')],
        extra_compile_args   = CXXFLAGS,
        extra_link_args      = LDFLAGS,
    )
    cIndexGenome = NumpyExtension(
        'stempy._%s_build._index_genome' % _build,
        module_index_srcs,
        include_dirs         = INCLUDE_DIRS,
        library_dirs         = LIBRARY_DIRS,
        libraries            = LIBRARIES,
        define_macros        = EXTRA_DEFINES + [('STEME_INDEX_MODULE_NAME', '_index_genome'), ('STEME_INDEX_MODULE_TYPE', 1)],
        extra_compile_args   = CXXFLAGS,
        extra_link_args      = LDFLAGS,
    )

    setup(
        name                  = 'STEME',
        version               = read('python', 'stempy', 'VERSION').strip().split('-')[0],
        description           = 'STEME: an accurate efficient motif finder for large data sets.',
        long_description      = read('python', 'stempy', 'README'),
        author                ='John Reid',
        author_email          ='johnbaronreid@netscape.net',
        license               = 'BSD',
        url                   ='http://sysbio.mrc-bsu.cam.ac.uk/johns/STEME/',
        classifiers           = [
                                    'Development Status :: 5 - Production/Stable',
                                    'Environment :: Console',
                                    'Intended Audience :: Developers',
                                    'Intended Audience :: Science/Research',
                                    'License :: OSI Approved :: BSD License',
                                    'Operating System :: MacOS :: MacOS X',
                                    'Operating System :: POSIX',
                                    'Operating System :: Microsoft :: Windows',
                                    'Programming Language :: Python',
                                    'Programming Language :: C++',
                                    'Topic :: Scientific/Engineering',
                                    'Topic :: Scientific/Engineering :: Mathematics',
                                    'Topic :: Scientific/Engineering :: Bio-Informatics',
                                    'Topic :: Utilities',
                                ],

        packages              = find_packages(where='python'),
        package_dir           = { '' : 'python' },
        package_data          = {
                                    'stempy': [
                                        'README',
                                        'LICENSE',
                                        'VERSION',
                                        'static/style.css',
                                        'templates/job.html',
                                        'templates/layout.html',
                                        'templates/motif.html',
                                        'templates/options.html',
                                        'templates/scan-stats.html',
                                        'templates/wide-layout.html',
                                    ],
                                    'stemewebapp': [
                                        'STEME.wsgi',
                                        'INSTALL',
                                        'VERSION',
                                        'static/style.css',
                                        'static/logo.png',
                                        'static/logo.svg',
                                        'templates/_formhelpers.html',
                                        'templates/home.html',
                                        'templates/layout.html',
                                        'templates/listjobs.html',
                                        'templates/newjob.html',
                                    ],
                                },
        include_package_data  = True,
        install_requires      = ['cookbook', 'numpy', 'matplotlib', 'weblogo'],
        scripts               = [
                                    'python/scripts/steme',
                                    'python/scripts/steme-em',
                                    'python/scripts/steme-find-spacings',
                                    'python/scripts/steme-motif-logos',
                                    'python/scripts/steme-pwm-scan',
                                    'python/scripts/steme-scan-stats',
                                    'python/scripts/steme-seqs-2-motif',
                                    'python/scripts/steme-spacing-analysis',
                                    'python/scripts/stemewebapp-create-db',
                                    'python/scripts/stemewebapp-remove-job',
                                ],
        ext_modules           = [cStempy, cIndex, cIndexGenome],

        # 2to3 invocation
        cmdclass              = {'build_py': build_py},
    )





if __name__ == '__main__':
    main()

