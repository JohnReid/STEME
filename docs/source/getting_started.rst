
.. index:: installation
.. _installation:

Installation
============
These installation instructions are specifically for an Ubuntu Linux system but should work with minor 
changes on other Linux systems and with some changes on other operating systems (Windows, MacOS, etc..). 
For other Linux systems you will need to replace any invocations of ``apt-get`` with your system's package
manager. For other operating systems you will need to work out the corresponding method of install.
The instructions install packages using ``apt-get`` where possible, otherwise they are installed in ``$HOME/local/``.
There are also some `Darwin/MacOS specific installation instructions`_.

I have written these instructions to install STEME into a virtual Python environment using the virtualenv_ tool.
This has the benefits of isolating the STEME install from the rest of your system and can be done without
administrator privileges. 

.. _virtualenv: http://www.virtualenv.org/en/latest/index.html


.. index:: prerequisites

Prerequisites
-------------
The following software should be downloaded and unpacked/installed somewhere on your system:

- GCC_: I have tested STEME with g++ versions 4.5.1 and 4.6.3 but other recent versions should work. I have had
  reports that gcc 4.2.1 works on a Mac but that gcc 4.0.1 does not. Every release of
  g++ seems to give it a stricter interpretation of the C++ standard so
  perhaps it is not a good idea to try a version too much newer. Other modern 
  compilers such as Microsoft C++ should work but I have not tried them. Please let me know
  if you try other compilers and run into any issues. 
  
  Make sure you have gcc (with g++) installed::
    
    sudo apt-get install gcc build-essential
  
  If gcc is not already installed and you do not have administrator privileges you could install it locally from source.

  
- `Python 2.7`_: Other versions of python 2 will probably work such as 2.5 or 2.6. Install Python 
  with virtualenv_::

	sudo apt-get install python-dev python-virtualenv

  STEME requires the following Python packages:
  
  * numpy
  * matplotlib
  * Bio
  * corebio
  * weblogolib
  * cookbook

  Recent versions of weblogolib contain the corebio package which has been removed from PyPI so you probably won't
  need to install corebio separately.
  We will install some of these packages into the operating system as they are fairly standard and widely used. In fact you
  may well have them installed already::
    
    sudo apt-get install python-numpy python-matplotlib python-biopython
  
  Now create a virtual environment to install STEME into and update `PATH` so it is used::
  
	mkdir -p $HOME/local
	virtualenv --system-site-packages $HOME/local/steme-virtualenv
	export PATH=$HOME/local/steme-virtualenv/bin:$PATH
	
  Now install the other requirements into the virtual environment. You will not need to install corebio if you use
  a recent version of weblogo::
  
    pip install corebio
    pip install weblogo
    pip install cookbook
    pip install jinja2
    
  The jinja2 package is not absolutely required but with it STEME will produce easy to read HTML output.
  
    
- `Boost C++ libraries`_: I have tested STEME using version 1.47, 1.48 and 1.49 of boost, 
  although any recent version should work. If you do not already have the libraries installed,
  first download them::
  
    mkdir -p $HOME/local/src
    cd $HOME/local/src
    wget \
      http://sourceforge.net/projects/boost/files/boost/1.49.0/boost_1_49_0.tar.bz2/download \
      -O boost_1_49_0.tar.bz2
    tar jxf boost_1_49_0.tar.bz2
     

  Install boost following the commands_ given at the Boost website and update your 
  ``LD_LIBRARY_PATH`` environment variable so that the shared objects can be found
  at runtime::
  
    cd $HOME/local/src/boost_1_49_0
    ./bootstrap.sh --help
    ./bootstrap.sh --prefix=$HOME/local --with-libraries=python
    ./b2 --ignore-config install release
    export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
  



- `SeqAn sequence analysis library`_: The SeqAn C++ library provides the suffix array implementation that
  STEME uses. I have used the latest SVN version (r12013) of SeqAn, although any recent version should work.
  You will need subversion to retrieve it and might need cmake to install it::

    sudo apt-get install subversion
    cd $HOME/local/src
    svn co http://svn.mi.fu-berlin.de/seqan/trunk/seqan
    cd seqan
    python util/bin/build_forwards.py core/include/seqan
    python util/bin/build_forwards.py extras/include/seqan



- `FFTW3`_: A library for computing discrete Fourier transforms. Install it::

    sudo apt-get install libfftw3-dev
    
  If you do not already have FFTW3 installed and are not running Ubuntu or you do not have 
  administrative privileges you can install `FFTW3 from source`_::
    
    cd $HOME/local/src
    wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.2.tar.gz
    tar zxf fftw-3.3.2.tar.gz
    cd fftw-3.3.2
    CFLAGS=-fPIC ./configure --prefix=$HOME/local
    make
    make install
    
    
.. _GCC: http://gcc.gnu.org/
.. _Python 2.7: http://www.python.org/
.. _Boost C++ libraries: http://www.boost.org/
.. _commands: http://www.boost.org/doc/libs/1_49_0/more/getting_started/unix-variants.html#easy-build-and-install
.. _SeqAn sequence analysis library: http://www.seqan.de/
.. _FFTW3: _http://www.fftw.org/
.. _FFTW3 from source: _http://www.fftw.org/download.html






.. index:: download

Download STEME
--------------

If you have not already done so, download STEME and unpack it locally. You will need to replace
``*.*.*`` with the version of STEME you want to install (check at PyPI_)::

    cd $HOME/local/src
    wget http://pypi.python.org/packages/source/S/STEME/STEME-*.*.*.tar.gz
    tar zxf STEME-*.*.*.tar.gz
    cd $HOME/local/src/STEME-*.*.*

.. _PyPI: http://pypi.python.org/pypi/STEME/





.. index:: build environment

Configure, build, install
-------------------------

We are ready to configure STEME. STEME uses aksetup for installation, which means that
installation should be easy and quick. Try::
  
    python configure.py --help

to examine the possible options. By the way, if a configuration option says ``several ok``,
then you may specify several values, separated by commas. We need to tell STEME
where the boost and seqan C++ libraries are::

    python configure.py \
      --seqan-dir=$HOME/local/src/seqan/ \
      --boost-inc-dir=$HOME/local/include \
      --boost-lib-dir=$HOME/local/lib

Configuration is obtained from files in this order::

    /etc/aksetup-defaults.py
    $HOME/.aksetup-defaults.py
    $PACKAGEDIR/siteconf.py

Once you've run configure, you can copy options from your ``siteconf.py`` file to
one of these files, and you won't ever have to configure them again manually.
In fact, you may pass the options ``--update-user`` and ``--update-global`` to
configure, and it will automatically update these files for you. This is particularly 
handy if you want to perform an unattended or automatic installation via pip_ or easy_install_.

Now install STEME::
    
    python setup.py install

This can take some time to compile. When it has finished, check that STEME has been successfully installed::

    steme --help

You should see a list of STEME's runtime options. 

.. _easy_install: http://packages.python.org/distribute/easy_install.html
.. _pip: http://pypi.python.org/pypi/pip





.. index:: Darwin/MacOS specific installation

Darwin/MacOS specific installation instructions
-----------------------------------------------

I have installed STEME successfully on MacOS 10.6.8 using Darwin gcc 4.2.1. Here are some notes that
might help you if you run into any problems.

* On my MacOS I do not have ``wget``, you should replace it in the above instructions with ``curl -O``. This
  may not work for the boost download link which you can download manually.


* The environment variable ``LD_LIBRARY_PATH`` is called ``DYLD_LIBRARY_PATH`` on the Mac 
  so this must be changed in the above.


* When installing packages using pip, you may see errors such as::

    IndentationError: unindent does not match any outer indentation level

  You need to edit the file mentioned and remove any blank lines at the end.


* If you see `errors <http://superuser.com/questions/242190/how-to-install-matplotlib-on-os-x>`_
  related to freetype font header files when installing matplotlib,
  you might be able to work around them with::

    # see: http://superuser.com/questions/242190/how-to-install-matplotlib-on-os-x
    export LDFLAGS="-L/usr/X11/lib"
    export CFLAGS="-I/usr/X11/include -I/usr/X11/include/freetype2 -I/usr/X11/include/libpng12"


* You may need to add ``-arch i386 -arch ppc -arch x86_64`` to CFLAGS when compiling boost and FFTW3.
