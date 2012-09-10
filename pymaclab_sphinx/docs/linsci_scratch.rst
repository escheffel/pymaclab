.. index:: cloud; installation

=======================================
Building a Linux Scientific Environment
=======================================

Requirements
============
* Python >= 2.5 or Python 3
* `Sphinx <http://sphinx.pocoo.org/>`_ 1.1 or newer.

Installing
==========
* To install from pypi using pip::

   pip install cloud_sptheme

* To install from pypi using easy_install::

   easy_install cloud_sptheme

* To install from source using ``setup.py``::

    python setup.py build
    sudo python setup.py install

.. index:: readthedocs.org; installation on

ReadTheDocs
===========
To use this theme on `<http://readthedocs.org>`_:

1. If it doesn't already exist, add a pip ``requirments.txt`` file to your documentation (e.g. alongside ``conf.py``).
   It should contain a minimum of the following lines::

       sphinx
       cloud_sptheme

   ... as well as any other build requirements for your project's documentation.

2. When setting up your project on ReadTheDocs, enter the path to ``requirements.txt``
   in the *requirements file* field on the project configuration page.

3. ReadTheDocs will now automatically download the latest version of :mod:`!cloud_sptheme`
   when building your documentation.

Documentation
=============
The latest copy of this documentation should always be available at:
    `<http://packages.python.org/cloud_sptheme>`_

If you wish to generate your own copy of the documentation,
you will need to:

1. Install `Sphinx <http://sphinx.pocoo.org/>`_ (1.1 or better)
2. Download the :mod:`!cloud_sptheme` source.
3. From the source directory, run ``python setup.py build_sphinx -E``.
4. Once Sphinx is finished, point a web browser to the file :samp:`{SOURCE}/build/sphinx/html/index.html`.
