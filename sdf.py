#-------------------------------------------------------------------------------
# Name:        sdf
# Purpose:     operations with sdf files (reading)
#
# Author:      Pavel Polishchuk
#
# Created:     11.01.2013
# Copyright:   (c) Pavel Polishchuk 2013
# Licence:     GPLv3
#-------------------------------------------------------------------------------

import os, sys

from mols import Mol3 as Mol

formal_charges_table = {'0': 0,
                        '1': 3,
                        '2': 2,
                        '3': 1,
                        '4': 0,
                        '5': -1,
                        '6': -2,
                        '7': -3}

def ReadSDF(fname):
    """
    INPUT: sdf-filename
    OUTPUT: dict of molecules, where key is the title of the moleculae taken from the first line of mol-record
    """

    def MolstrToMol(molstr):
        mol = Mol()
        mol.title = molstr[0]
        natoms = int(molstr[3][0:3])
        nbonds = int(molstr[3][3:6])
        # read atoms
        id = 0
        for line in molstr[4:4+natoms]:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            label = line[30:33].strip()
            formal_charge = line[36:39].strip()
            id += 1
            mol.AddAtom(id, label, x, y, z, formal_charges_table[formal_charge])
        # read bonds
        for line in molstr[4+natoms:4+natoms+nbonds]:
            id1 = int(line[0:3])
            id2 = int(line[3:6])
            bond_type = int(line[6:9])
            mol.AddBond(id1, id2, bond_type)
        return mol

    mols = {}
    molstr = []
    f = open(fname)
    for line in f:
        if line.find("$$$$") < 0:
            molstr.append(line.rstrip())
        else:
            m = MolstrToMol(molstr)
            mols[m.title] = m
            molstr = []
    return mols
