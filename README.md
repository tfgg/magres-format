magres-format
=============

Code for parsing CASTEP ab-initio magnetic resonance file format. See more at http://www.ccpnc.ac.uk/pmwiki.php/CCPNC/Fileformat

Installing
----------

Clone the repository or download and extract the .zip file somewhere. From the command line, run:

    python setup.py install

and the Python module and associated scripts should now be installed.

Conversion script usage
-----------------------

The convertoldmagres.py script installed by the above command will convert an old-style CASTEP magres file to
the new-style format for use with the new tools. You use it from the command line like:

    convertoldmagres.py sample.magres > sample.new.magres

and optionally with the associated job's .castep file, to capture the lattice information,

    convertoldmagres.py sample.magres sample.castep > sample.new.magres

Python module usage
-------------------

From inside Python you can import the magres.format module and use the MagresFile class.


     >>> from magres.format import MagresFile
     >>> magres_file = MagresFile(open('samples/T1Si0.magres'))
     >>> magres_file
     >>> <magres.format.MagresFile object at 0x26ad690>

The 'data_dict' member variable contains all the information parsed from the file. You can also request the data
serialized as JSON by the .as_json() method and load JSON data with the .load_json() method.
