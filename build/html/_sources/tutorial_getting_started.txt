Getting started
===============

Installing
----------

You can either use git to clone the latest version of the repository or download the repository as a .zip.

>>> git clone git://github.com/tfgg/magres-format.git
>>> cd magres-format

or

>>> wget https://github.com/tfgg/magres-format/archive/master.zip
>>> unzip master.zip 
>>> cd magres-format-master

Now run setup.py to install the module and associated scripts. If you want to install it system wide, run

>>> sudo python setup.py install

If you want to install it locally (e.g. if you don't have sudo), run

>>> python setup.py install --user

and it should now be available.

Python package
--------------

The package comes with a number of Python modules, documented separately.

.. toctree::
   :maxdepth: 2

   atoms
   format

Import the modules into your code with the standard import commands

>>> import magres.atoms
>>> import magres.format

or just import individual components

>>> from magres.atoms import MagresAtoms


