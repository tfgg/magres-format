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


Helper scripts
--------------

The package comes with a number of helper scripts that can be run from the command line.

Dump JSON: magresjson.py
~~~~~~~~~~~~~~~~~~~~~~~~

This will give you a JSON dump of a given .magres text file.

>>> magresjson.py sample.magres > sample.magres.json

Merge .magres files: magresmerge.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will merge together several given .magres (e.g. different perturbing atoms in an indirect spin-spin coupling calculation) and output the result.

>>> magresmerge.py */*.magres > ethanol-jc.magres

Quickly analyse J-couplings: extract-all-jc-compare.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> extract-all-jc-compare.py [directory] [species] [index]

This will find all two-way indirect spin-spin coupling calculations below [directory] between atom [species] [index] and other atoms and print them out for analysis. For example, extracting all couplings to the C1 atom in ethanol:

.. code::

  >>> extract-all-jc-compare.py . C 1
  #atm1 atm2  isc isc_fc  isc_spin  isc_orbital_p isc_orbital_d   r
  13C1  13C2  41.06 40.66 1.21  -0.96 0.16  1.51993253781
  13C2  13C1  41.12 40.68 1.22  -0.93 0.15  1.51993253781
  1H1 13C1  34.49 33.80 0.14  0.39  0.15  1.09435820235
  13C1  1H1 34.62 33.87 0.07  0.59  0.09  1.09435820235
  1H2 13C1  35.33 34.68 0.14  0.39  0.12  1.09508258166
  13C1  1H2 35.50 34.75 0.07  0.58  0.10  1.09508258166
  13C1  1H3 35.96 35.23 0.07  0.56  0.10  1.09218503188
  1H3 13C1  35.80 35.16 0.14  0.37  0.14  1.09218503188
  13C1  1H4 1.05  1.06  0.01  0.06  -0.09 2.1653422802
  1H4 13C1  1.04  1.09  0.02  -0.01 -0.05 2.1653422802
  1H5 13C1  0.18  0.25  0.02  -0.03 -0.05 2.16077342428
  13C1  1H5 0.23  0.23  0.02  0.07  -0.09 2.16077342428
  1H6 13C1  0.22  0.24  -0.01 0.02  -0.04 2.62546755589
  13C1  1H6 0.30  0.24  -0.00 0.13  -0.07 2.62546755589
  17O1  13C1  0.95  2.51  -0.06 -1.45 -0.06 2.44843140426
  13C1  17O1  0.96  2.54  -0.03 -1.51 -0.04 2.44843140426

Samples and example scripts
---------------------------

There are some sample .magres files in the samples/ directory and some example Python processing scripts in the examples/ directory. Run the example scripts from the package root, like

>>> python examples/atoms_ms.py

or they will not find the correct .magres files.
