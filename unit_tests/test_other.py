# -*- coding: utf-8 -*-

# =============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
# christian.fernandez@ensicaen.fr
#
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
#
# See Licence.txt in the main directory
# =============================================================================

"""Short description.

A longer description.

"""

from __future__ import print_function

from magres.atoms import MagresAtoms, MagresAtomsView
ethanol_atoms = MagresAtoms.load_magres('../samples/ethanol-all.magres')

ethanol_atoms.calculate_bonds()
