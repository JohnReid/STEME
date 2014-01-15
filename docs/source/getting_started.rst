
.. index:: installation
.. _installation:

Installation
============

These installation instructions are specifically for an Ubuntu Linux system 
with a bash shell but should work with minor changes on other Linux systems and
with some changes on other operating systems (Windows, MacOS, etc..). 
For other Linux systems you will need to replace any invocations of
``apt-get`` with your system's package manager. For other operating systems you
will need to work out the corresponding method of install. The instructions
install packages using ``apt-get`` where possible, otherwise they are
installed in ``$HOME/local/``. There are also some `Darwin/MacOS specific
installation instructions`_. All of the prerequisites can be installed from
source if necessary but it is far easier to use your system's package
manager if possible.

I have written these instructions to install STEME into a virtual Python
environment using the virtualenv_ tool. This has the benefits of isolating the
STEME install from the rest of your system and can be done without
administrator privileges. 

.. _virtualenv: http://www.virtualenv.org/en/latest/index.html


.. index:: prerequisites

Prerequisites
-------------
The following software should be downloaded and unpacked/installed somewhere
on your system:

- GCC_: I have tested STEME with g++ versions 4.5.1 and 4.6.3 but other recent versions should work. I have had
  reports that gcc 4.2.1 works on a Mac but that gcc 4.0.1 does not. Every release of
  g++ seems to give it a stricter interpretation of the C++ standard so
  perhaps it is not a good idea to try a version too much newer. Other modern 
  compilers such as Microsoft C++ should work but I have not tried them. Please let me know
  if you try other compilers and run into any issues. 

  Make sure you have gcc (with g++) installed::

    sudo apt-get install gcc build-essential

  If gcc is not already installed and you do not have administrator privileges
  you could install it locally from source.


- `Boost C++ libraries`_: Versions more recent than 1.47 should work. The
  easiest way to install them is through your package manager::

    sudo apt-get install libboost-dev

  Otherwise you can download them and install the libraries following the
  commands_ given at the Boost website.


- `SeqAn sequence analysis library`_: The SeqAn C++ library provides the suffix
  array implementation that STEME uses. On Ubuntu you can install the library
  with::

    sudo apt-get install seqan-dev


- `FFTW3`_: A library for computing discrete Fourier transforms. Install it::

    sudo apt-get install libfftw3-dev


- `Python 2.7`_: Other versions of python 2 will probably work such as 2.5 or
  2.6. Install Python with virtualenv_::

    sudo apt-get install python-dev python-virtualenv

  STEME requires the following Python packages:

  * numpy
  * matplotlib
  * Bio
  * corebio
  * weblogolib
  * cookbook

  Recent versions of weblogolib contain the corebio package which has been
  removed from PyPI so you probably won't need to install corebio separately.
  We will install some of these packages into the operating system as they are
  fairly standard and widely used. In fact you may well have them installed
  already::

    sudo apt-get install python-numpy python-matplotlib python-biopython python-jinja2

  The jinja2 package is not absolutely required but with it STEME will produce
  easy to read HTML output.


.. _GCC: http://gcc.gnu.org/
.. _Python 2.7: http://www.python.org/
.. _Boost C++ libraries: http://www.boost.org/
.. _commands: http://www.boost.org/doc/libs/1_49_0/more/getting_started/unix-variants.html#easy-build-and-install
.. _SeqAn sequence analysis library: http://www.seqan.de/
.. _FFTW3: _http://www.fftw.org/
.. _FFTW3 from source: _http://www.fftw.org/download.html



Install STEME
-------------

Create a virtual environment to install STEME into and activate it::

    mkdir -p $HOME/local
    virtualenv --system-site-packages $HOME/local/steme-virtualenv
    source $HOME/local/steme-virtualenv/bin/activate

Install STEME::

    pip install steme

You will probably see a message about configuration not being run. Everything
should be fine unless you have built some of the libraries from source and
they are not in standard locations. In that case, you download the STEME
source package from PyPI_, unpack it and see the instructions in the
configuration_ section.

The install can take some time to compile. When it has finished, check that
STEME has been successfully installed::

    steme --help

You should see a list of STEME's runtime options. 

.. _PyPI: https://pypi.python.org/pypi/STEME/




.. index:: _configuration

Configuration
-------------

If you are installing STEME using headers in non-standard locations then you
will need to configure it first. STEME uses aksetup for installation, which
means that this should be easy and quick. Try::

    python configure.py --help

to examine the possible options. By the way, if a configuration option says
``several ok``, then you may specify several values, separated by commas.
For example, we might need to tell STEME where the boost and SeqAn C++
libraries are::

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

Once STEME is configured, you can install it using the normal method::

    python setup.py install


.. _easy_install: http://packages.python.org/distribute/easy_install.html
.. _pip: http://pypi.python.org/pypi/pip




.. index:: Darwin/MacOS specific installation

Darwin/MacOS specific installation instructions
-----------------------------------------------

I have installed STEME successfully on MacOS 10.6.8 using Darwin gcc 4.2.1. Here are some notes that
might help you if you run into any problems.

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
