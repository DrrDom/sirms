#-------------------------------------------------------------------------------
# Name:        mols
# Purpose:     Mol class to operate with molecules
#
# Author:      Pavel Polishchuk
#
# Created:     11.01.2013
# Copyright:   (c) Pavel Polishchuk 2013
# Licence:     GPLv3
#-------------------------------------------------------------------------------

from collections import defaultdict

class Mol3:
    def __init__(self):
        self.atoms = {}
        self.bonds = dict()
        self.title = ""
        self.component_ratios = {}  # {id1: value1, id2, value2}
    def AddAtom(self, id, label, x, y, z):
        self.atoms[id] = {'label': label, 'x': x, 'y': y, 'z': z, 'property': {}, 'component_id': None}
        self.bonds[id] = dict()
    def AddBond(self, id1, id2, bond_type):
        if id1 not in self.bonds.keys():
            self.bonds[id1] = dict()
        if id2 not in self.bonds.keys():
            self.bonds[id2] = dict()
        self.bonds[id1][id2] = self.bonds[id2][id1] = bond_type
    def GetBondType(self, id1, id2):
        try:
            return(self.bonds[id1][id2])
        except KeyError:
            return(0)
