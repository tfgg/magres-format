Processing electric field gradients
========================

Gehlenite: A disordered material
--------------------------------

Gehlenite is a disordered aluminosilicate consisting of two types of aluminium site and one type of silicon site. We want to look at aluminium EFGs to quantify the amount of disorder in the material, in particular the degree of violation of the "Loewenstein" rule (avoidance of Al-O-Al linkages) and its structural effects. This section is based on analysis performed in `Elucidation of the Al/Si Ordering in Gehlenite Ca2Al2SiO7 by Combined 29Si and 27Al NMR Spectroscopy/Quantum Chemical Calculations <http://pubs.acs.org/doi/abs/10.1021/cm3016935>`_.

Two sets of EFG calculations have been performed. Both have 25 relaxed structures, representing five different central aluminium sites with five random occupations of other sites. The first set, "low", has no Al-O-Al or Si-O-Si bonds present, the other "random" has some present. We want to aggregate aluminium EFGs according to their site type (indexed by number of silicon second-neighbours) and compute statistics.

`Download the files here <http://tfgg.github.io/magres-format/workshop/tutorial_disordered/gehlenite_efgs.zip>`_ and unzip them somewhere. `Download the example script <http://tfgg.github.io/workshop/tutorial_disordered/Cq_dist.py>`_ and put it in the same directory.

.. code::

  curl -L "http://tfgg.github.io/magres-format/workshop/tutorial_disordered/gehlenite_efgs.zip" > gehlenite_efgs.zip
  unzip gehlenite_efgs.zip
  curl -L "http://tfgg.github.io/magres-format/workshop/tutorial_disordered/Cq_dist.py" > Cq_dist.py

Inspect Cq_dist.py with your favourite text editor and try running it on the two sets of structures:

.. code::

  python Cq_dist.py files/orig
  python Cq_dist.py files/low

Note the means and standard deviations for the sites. What does this tell you about the presence of Loewenstein rule violating sites and disorder?

Some key techniques used in this script:

* We used magres.utils.find_all to quickly locate all the .magres files in the given directory and make a MagresAtoms object for each.

* We loop over each atoms object and loop over every aluminium atom in each by using the atoms.species('Al') selector.

* In order to classify the aluminium site as T1 or T2 (see paper), we count the number of aluminium and silicon neighbours within 3.5 Angstroms (chosen for this specific system) using the atoms.species(['Al', 'Si']).within(Al_atom, 3.5) selector.

* We further refine this list of neighbours to just silicon atoms in order to count them for n_si using the neighbours.species('Si') selector.

* We access the C_q value for the aluminium atom by access the efg.Cq property on it. This is available on all atoms with EFG data.

Plotting with matplotlib
~~~~~~~~~~~~~~~~~~~~~~~~

If the matplotlib Python module is installed on your computer, you could try plotting a histogram of the EFGs with the following code at the end of the script:

.. code:: python

  import matplotlib.pyplot as plt
  plt.hist(Al_Cqs['T1'][0], bins=20)
  plt.hist(Al_Cqs['T1'][1], bins=20)
  plt.hist(Al_Cqs['T1'][2], bins=20)
  plt.hist(Al_Cqs['T1'][3], bins=20)
  plt.hist(Al_Cqs['T1'][4], bins=20)
  plt.show()

Calculating distortion parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modify the script to calculate Ghose's tetrahedral distortion parameter (Ghose, S.; Tsang, T. Am. Mineral. 1973, 58, 748): the mean tangent of the absolute deviation of the O-Al-O bond angles from the tetrahedral angle (109.47°). Does this correlate with the aluminium Cqs?

A working script which does this `is available at here <http://tfgg.github.io/magres-format/workshop/tutorial_disordered/Cq_ghose.py>`_. Run it like

.. code:: bash

  python Cq_ghose files/low > ghose.dat

and use gnuplot to plot a scatterplot of tetrahedral distortion vs. aluminium Cq.
